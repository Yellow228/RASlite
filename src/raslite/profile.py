"""标准步进法 (Standard Step Method) 一维恒定流水面线推算。

理论: HEC-RAS Hydraulic Reference Manual v6.0, Ch.2
  能量方程:  Z2 + Y2 + a2*V2^2/2g = Z1 + Y1 + a1*V1^2/2g + he   (式 2-1)
  损失:      he = L*Sf_bar + C*|a2*V2^2/2g - a1*V1^2/2g|         (式 2-2)
  长度:      L = sum(Li*Qi_avg)/sum(Qi_avg)  (左滩/主槽/右滩流量加权, 式 2-3)
  摩阻坡度:  默认平均输水能力法 (Q1+Q2)^2/(K1+K2)^2, 另有 3 种可选 (式 2-13~2-16)

迭代 (HEC-RAS 计算流程):
  试值1: 上一断面水深投影;
  试值2: 假定值 + 0.70*误差;
  试值3+: 割线法外推 (式 2-18);
  收敛判据: |计算值 - 假定值| <= 容差 (默认 0.003 m = 0.01 ft)。
  不收敛时 HEC-RAS 默认取临界水深并报警 —— 此处复刻该行为。

系数约定 (对齐 HEC-RAS):
  - reach_lengths / Exp-Cntr 系数存于上游断面, 应用于该断面与其下游
    断面之间的河段 (即成对断面取上游断面的系数);
  - 收缩/扩散判定: 下游流速水头 > 上游 -> 收缩 (用 c_contract), 反之为扩散。

方向约定:
  断面 chainage 沿流向递增 (下游大)。亚临界流从最下游已知边界向上游推,
  超临界流从最上游已知边界向下游推。

v1 范围: 天然河道段 (无桥涵堰), 亚临界 / 超临界; 混合流 (比力判别) 二期。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import G
from .geometry import CrossSection
from .hydraulics import critical_ws

_TOL_DEFAULT = 0.003    # 0.01 ft
_MAX_ITER_DEFAULT = 40

SBAR_METHODS = ("avg_conveyance", "avg_friction", "geo_mean", "harmonic")
_PARTS = ("lob", "ch", "rob")


@dataclass
class SectionResult:
    """单断面计算结果 (字段名对齐 HEC-RAS Detailed Output 便于核对)。"""
    section: CrossSection
    q: float
    ws: float                 # 水面高程 WSE
    crit_ws: float = float("nan")   # 临界水面高程 (不可用时保持 NaN)
    area: float = 0.0
    velocity: float = 0.0
    hv: float = 0.0           # alpha*V^2/2g
    total_head: float = 0.0   # H = WS + hv
    conveyance: float = 0.0
    alpha: float = 1.0
    top_width: float = 0.0
    froude: float = 0.0
    warnings: list = field(default_factory=list)

    @property
    def depth(self) -> float:
        return self.ws - self.section.z_min

    @property
    def critical_depth(self) -> float:
        return self.crit_ws - self.section.z_min


def _sbar(method: str, q_us: float, q_ds: float, k1: float, k2: float,
          sf1: float, sf2: float) -> float:
    """代表性摩阻坡度 (式 2-13 ~ 2-16)。q1/q2 为两断面各自流量。"""
    if method == "avg_conveyance":
        return ((q_us + q_ds) / (k1 + k2)) ** 2
    if method == "avg_friction":
        return 0.5 * (sf1 + sf2)
    if method == "geo_mean":
        return math.sqrt(max(sf1 * sf2, 0.0))
    if method == "harmonic":
        if sf1 + sf2 <= 0:
            return 0.0
        return 2.0 * sf1 * sf2 / (sf1 + sf2)
    raise ValueError(f"未知摩阻坡度算法: {method}, 可选 {SBAR_METHODS}")


class SteadyProfile:
    """一维恒定流水面线计算器。

    Parameters
    ----------
    sections : [CrossSection, ...]  顺序任意, 内部按 chainage 升序排序;
        每个断面需带 chainage 属性 (沿流向里程, 下游大)
    q : float  恒定流量 (沿程不变; 沿程汇流二期)
    regime : 'subcritical' | 'supercritical'
    boundary : tuple
        ('ws', elev)        已知水面高程
        ('critical',)       临界水深
        ('normal', S0)      正常水深 (需给河床比降)
        ('rating', fn)      水位流量关系 fn(q)->ws
    sbar_method : 摩阻坡度算法, 默认 'avg_conveyance'
    tol / max_iter : 迭代容差 (m) 与最大迭代次数 (默认同 HEC-RAS)
    """

    def __init__(self, sections, q, regime="subcritical", boundary=("ws",),
                 sbar_method="avg_conveyance", tol=_TOL_DEFAULT,
                 max_iter=_MAX_ITER_DEFAULT):
        if regime not in ("subcritical", "supercritical"):
            raise ValueError("regime 必须为 subcritical / supercritical")
        if sbar_method not in SBAR_METHODS:
            raise ValueError(f"sbar_method 必须为 {SBAR_METHODS}")
        for s in sections:
            if not hasattr(s, "chainage"):
                raise TypeError(f"断面 {s.name} 缺少 chainage 属性 (沿流向里程)")
        self.sections = sorted(sections, key=lambda s: s.chainage)
        if isinstance(q, (list, tuple)):
            if len(q) != len(self.sections):
                raise ValueError("q 序列长度必须与断面数一致")
            self.qs = [float(v) for v in q]
        else:
            self.qs = [float(q)] * len(self.sections)
        self.q = self.qs[0]  # 兼容旧接口 (标量情形)
        self.regime = regime
        self.boundary = boundary
        self.sbar_method = sbar_method
        self.tol = tol
        self.max_iter = max_iter

    # ---------------- 边界条件 ----------------

    def _boundary_ws(self, sec: CrossSection, q: float) -> float:
        kind = self.boundary[0]
        if kind == "ws":
            return float(self.boundary[1])
        if kind == "critical":
            ws, _ = critical_ws(sec, q)
            return ws
        if kind == "normal":
            from .hydraulics import normal_ws
            return normal_ws(sec, q, float(self.boundary[1]))
        if kind == "rating":
            return float(self.boundary[1](q))
        raise ValueError(f"未知边界条件类型: {kind}")

    # ---------------- 河段长度与损失 ----------------

    def _reach_length(self, sec_us: CrossSection, p_us, p_ds,
                      q_us: float, q_ds: float) -> float:
        """流量加权断面间距 (式 2-3)。sec_us 为成对断面的上游断面,
        其 reach_lengths 指向下游断面; Q_i 取上下断面分区流量的算术平均。"""
        k_us = p_us.k_by_part()
        k_ds = p_ds.k_by_part()
        num = den = 0.0
        for part, L in zip(_PARTS, sec_us.reach_lengths):
            qu = q_us * k_us[part] / p_us.conveyance
            qd = q_ds * k_ds[part] / p_ds.conveyance
            q_avg = 0.5 * (qu + qd)
            num += L * q_avg
            den += q_avg
        if den <= 1e-9:
            return float(sec_us.reach_lengths[1])
        return num / den

    def _iterate(self, sec_known, p_known, sec_unknown,
                 q_known: float, q_unknown: float):
        """标准步进迭代。返回 (SectionResult, SectionProps)。"""
        # 成对断面的上游断面及其系数 (HEC-RAS 约定)
        if self.regime == "subcritical":
            sec_us = sec_unknown          # known=下游, unknown=上游
        else:
            sec_us = sec_known            # known=上游, unknown=下游

        def ws_calc_of(ws_guess: float):
            p = sec_unknown.props(ws_guess)
            # 上/下游 props 与各自断面流量
            if self.regime == "subcritical":
                p_us, p_ds = p, p_known
                q_us, q_ds = q_unknown, q_known
                sec_L = sec_unknown
            else:
                p_us, p_ds = p_known, p
                q_us, q_ds = q_known, q_unknown
                sec_L = sec_known
            L = self._reach_length(sec_L, p_us, p_ds, q_us, q_ds)
            sf_us = (q_us / p_us.conveyance) ** 2
            sf_ds = (q_ds / p_ds.conveyance) ** 2
            sbar = _sbar(self.sbar_method, q_us, q_ds,
                         p_us.conveyance, p_ds.conveyance, sf_us, sf_ds)
            hv_us = p_us.velocity_head_of(q_us)
            hv_ds = p_ds.velocity_head_of(q_ds)
            # 下游流速水头 > 上游 -> 收缩 (HEC-RAS 规则)
            if hv_ds > hv_us:
                c = sec_us.c_contract
            else:
                c = sec_us.c_expand
            he = L * sbar + c * abs(hv_ds - hv_us)
            if self.regime == "subcritical":
                h_new = p_known.total_head_of(q_known) + he
            else:
                h_new = p_known.total_head_of(q_known) - he
            return h_new - p.velocity_head_of(q_unknown), p

        # ---- 试值序列 (HEC-RAS 策略) ----
        z1, z2 = sec_known.z_min, sec_unknown.z_min
        g1 = z2 + (p_known.ws - z1)                 # 试值1: 水深投影
        c1, p1 = ws_calc_of(g1)
        err1 = c1 - g1
        if abs(err1) <= self.tol:
            ws_final, p_final = c1, p1
        else:
            g2 = g1 + 0.70 * err1                   # 试值2: 0.7 倍误差修正
            c2, p2 = ws_calc_of(g2)
            err2 = c2 - g2
            if abs(err2) <= self.tol:
                ws_final, p_final = c2, p2
            else:
                ws_final, p_final = None, None
                gu, gl = g1, g2                     # 旧/新试值
                eu, el = err1, err2
                # 限幅上界: 封闭断面取顶高程; 开放断面不设上界
                if sec_unknown.z_max > sec_unknown.z_min + 1e-3:
                    g_hi = sec_unknown.z_max - 1e-6
                else:
                    g_hi = float("inf")
                for _ in range(self.max_iter - 2):
                    if abs(el - eu) < 1e-15:
                        break
                    gn = gl - el * (gl - gu) / (el - eu)   # 割线法 (式 2-18)
                    gn = min(max(gn, sec_unknown.z_min + 1e-9), g_hi)
                    cn, pn = ws_calc_of(gn)
                    en = cn - gn
                    if abs(en) <= self.tol:
                        ws_final, p_final = cn, pn
                        break
                    gu, gl = gl, gn
                    eu, el = el, en
                if ws_final is None:
                    # 不收敛 -> HEC-RAS 行为: 默认取临界水深并报警
                    warns = []
                    try:
                        cw, _ = critical_ws(sec_unknown, q_unknown)
                    except Exception as e:
                        raise RuntimeError(
                            f"{sec_unknown.name}: 迭代不收敛且临界水深"
                            f"求解失败 ({e})") from e
                    p_c = sec_unknown.props(cw)
                    warns.append(
                        f"迭代 {self.max_iter} 次未收敛 (最后误差 "
                        f"{abs(el):.4f} m), 默认取临界水深 {cw:.3f}")
                    return self._make_result(sec_unknown, cw, p_c,
                                             q_unknown, warns), p_c

        # ---- 临界水深校核与默认处理 ----
        warns = []
        try:
            cw, _ = critical_ws(sec_unknown, q_unknown)
        except Exception as e:
            cw = None
            warns.append(f"临界水深求解失败 ({e})")
        if cw is not None:
            if self.regime == "subcritical" and ws_final < cw - 1e-9:
                warns.append(
                    f"计算水位 {ws_final:.3f} 低于临界 {cw:.3f}, "
                    f"默认取临界水深 (亚临界)")
                ws_final, p_final = cw, sec_unknown.props(cw)
            elif self.regime == "supercritical" and ws_final > cw + 1e-9:
                warns.append(
                    f"计算水位 {ws_final:.3f} 高于临界 {cw:.3f}, "
                    f"默认取临界水深 (超临界)")
                ws_final, p_final = cw, sec_unknown.props(cw)

        return self._make_result(sec_unknown, ws_final, p_final,
                                 q_unknown, warns), p_final

    def _make_result(self, sec, ws, p, q, warns) -> SectionResult:
        try:
            cw, _ = critical_ws(sec, q)
        except Exception:
            cw = float("nan")
        return SectionResult(
            section=sec, q=q, ws=ws, crit_ws=cw,
            area=p.area, velocity=q / p.area,
            hv=p.velocity_head_of(q), total_head=p.total_head_of(q),
            conveyance=p.conveyance, alpha=p.alpha,
            top_width=p.top_width, froude=math.sqrt(p.froude_sq_of(q)),
            warnings=warns,
        )

    # ---------------- 主流程 ----------------

    def run(self):
        """返回 [SectionResult, ...], 按推算顺序排列 (边界断面在前)。"""
        secs = self.sections
        if len(secs) < 2:
            raise ValueError("至少需要 2 个断面")

        if self.regime == "subcritical":
            order = list(reversed(secs))     # 最下游起, 向上游推
        else:
            order = list(secs)               # 最上游起, 向下游推

        # 各断面的流量 (与排序后的 sections 对齐)
        if not hasattr(self, "qs") or len(self.qs) != len(secs):
            base = getattr(self, "qs", [self.q])
            self.qs = [base[0]] * len(secs) if len(base) == 1 else base
        q_by_sec = {id(s): q for s, q in zip(secs, self.qs)}

        sec0 = order[0]
        q0 = q_by_sec[id(sec0)]
        ws0 = self._boundary_ws(sec0, q0)
        p0 = sec0.props(ws0)
        results = [self._make_result(sec0, ws0, p0, q0, [])]

        p_known = p0
        for i in range(1, len(order)):
            res, p_known = self._iterate(order[i - 1], p_known, order[i],
                                         q_by_sec[id(order[i - 1])],
                                         q_by_sec[id(order[i])])
            results.append(res)
        return results

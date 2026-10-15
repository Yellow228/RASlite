"""临界水深、正常水深与 Froude 数。

临界水深: 解 Fr^2(ws) = alpha*Q^2*T/(g*A^3) = 1。
  该条件由总水头 H = WS + alpha*V^2/2g 对 WS 取极小导出, 对任意
  断面形状严格成立 (矩形棱柱渠道退化为 yc=(q^2/g)^(1/3))。
  Fr^2 随水位总体单调递减, 用二分法求根, 稳定可靠。

  开放断面 (断面两端低于水面可无限抬升, z_max == 最低河底) 时,
  由河底附近逐级扩张上界直至 Fr^2 < 1。

  局限 (v1): 复式断面 T 在滩肩突变可能使 Fr^2(ws) 非单调 (多重
  临界水深现象), 此处取最浅的第一个解; 与 HEC-RAS 比能最小法
  可能有差异, L2 验证时核对。

正常水深: 解 Q = K(ws)*sqrt(S0) (曼宁公式), K 随 ws 单调增, 二分求根。

HEC-RAS 行为对齐说明:
  - 断面内找不到临界水深时, 本模块抛 CriticalDepthError;
    推算过程由 profile 层决定降级行为 (报警/取边界值)。
"""

from __future__ import annotations

from .geometry import CrossSection


class CriticalDepthError(ValueError):
    pass


def _lo_bound(cs: CrossSection) -> float:
    z = cs.z_min
    # 相对高程量纲取下限, 避免极浅水深时面积数值下溢
    return z + max(1e-6, abs(z) * 1e-9)


def critical_ws(cs: CrossSection, q: float, tol: float = 1e-9):
    """返回 (临界水面高程, 临界比能 E = H - z_min)。"""
    if q <= 0:
        raise ValueError("流量必须为正")

    def fr2(ws: float) -> float:
        return cs.props(ws).froude_sq_of(q)

    lo = _lo_bound(cs)
    if fr2(lo) <= 1.0:
        # 理论上极浅处 Fr^2 -> inf, 此情形仅出现在数值异常
        raise CriticalDepthError(f"{cs.name}: 河底附近 Froude 数异常")

    # 确定包含 Fr^2 = 1 的区间上界
    span = cs.z_max - cs.z_min
    if span > 1e-9:
        hi = cs.z_max
        if fr2(hi) > 1.0:
            raise CriticalDepthError(
                f"{cs.name}: 断面内无临界水深 (抬升至断面顶 Fr^2 仍 > 1), "
                f"断面过小或流量过大")
    else:
        # 开放断面: 上界逐级扩张
        d0 = max(1e-3, 0.001 * max(1.0, abs(cs.z_min)))
        hi = cs.z_min + d0
        it = 0
        while fr2(hi) > 1.0:
            d0 *= 2.0
            hi = cs.z_min + d0
            it += 1
            if it > 200:
                raise CriticalDepthError(
                    f"{cs.name}: 开放断面内未找到临界水深 (扩张超过上限)")

    a, b = lo, hi
    for _ in range(200):
        mid = 0.5 * (a + b)
        if (b - a) < tol:
            break
        if fr2(mid) > 1.0:
            a = mid
        else:
            b = mid
    ws = 0.5 * (a + b)
    p = cs.props(ws)
    return ws, p.total_head_of(q) - cs.z_min


def normal_ws(cs: CrossSection, q: float, s0: float, tol: float = 1e-9,
              max_iter: int = 200):
    """正常水深对应的水面高程 (二分法解 Q = K(ws)*sqrt(S0))。"""
    if q <= 0 or s0 <= 0:
        raise ValueError("流量与坡度必须为正")

    lo = _lo_bound(cs)
    hi = cs.z_max
    if hi <= lo:
        # 开放断面: 上界扩张至 K 足够
        d0 = max(1e-3, 0.001 * max(1.0, abs(cs.z_min)))
        hi = cs.z_min + d0
        it = 0
        while cs.props(hi).conveyance * s0 ** 0.5 < q:
            d0 *= 2.0
            hi = cs.z_min + d0
            it += 1
            if it > 200:
                raise CriticalDepthError(
                    f"{cs.name}: 断面内无正常水深 (扩张超过上限)")

    def f(ws: float) -> float:
        return cs.props(ws).conveyance * s0 ** 0.5 - q

    f_lo, f_hi = f(lo), f(hi)
    if f_lo > 0 or f_hi < 0:
        raise CriticalDepthError(
            f"{cs.name}: 断面内无正常水深 (f(lo)={f_lo:.3f}, f(hi)={f_hi:.3f})")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if (hi - lo) < tol:
            break
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def froude_squared(props, q: float) -> float:
    """Fr^2 = alpha*Q^2*T/(g*A^3)"""
    return props.froude_sq_of(q)

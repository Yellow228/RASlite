"""HEC-RAS 几何文件 (.g01) 解析器。

读取工程几何: 河段/断面 (站距-高程、糙率分区、滩肩、河段长度、
收缩/扩散系数) 与汇合口 (入流出流河段、汇合口长度),
转换为 raslite 的 Reach / Junction / CrossSection 对象。

单位: .g01 本身无单位标记, 单位记录在 .prj (English Units /
Metric Units)。parse_g01 的 unit 参数需调用方指定 (或用
read_prj_units 辅助)。

格式要点 (v1 支持范围):
  - Type RM Length L Ch R = 1 ,<RS> ,<LOB>,<CH>,<ROB>   (type 1 = 断面)
  - #Sta/Elev= N   随后 N 对站距/高程 (可跨行)
  - #Mann= N ,0 ,0 随后 N 组 (站距, n, 标志) 三元组; 过滤 n<=0
    或超大的占位条目
  - Bank Sta=<左滩肩>,<右滩肩>
  - Exp/Cntr=<扩散>,<收缩>
  - Junct Name= / Up River,Reach= / Dn River,Reach= / Junc L&A=<长>,<角>
    (Junc L&A 按出现顺序对应各入流河段)
  - 桥涵等其他节点 (type != 1) v1 跳过并告警
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from .geometry import CrossSection
from .network import Junction, Reach

_FT = 0.3048
_MI = 1609.344


def read_prj_units(prj_path: str) -> str:
    """从 .prj 读单位制: 返回 'ft' 或 'm'。"""
    with open(prj_path, encoding="latin-1") as f:
        for line in f:
            if line.startswith("English Units"):
                return "ft"
            if line.startswith("Metric Units") or line.startswith("SI Units"):
                return "m"
    return "ft"   # HEC-RAS 历史默认


@dataclass
class G01Section:
    rs: str
    points: list            # [(station, elev)] 已换算
    bank_left: float
    bank_right: float
    n_main: float
    n_left: list
    n_right: list
    reach_lengths: tuple    # (LOB, CH, ROB) 已换算
    c_contract: float
    c_expand: float
    description: str = ""
    levees: list = field(default_factory=list)   # [(side, sta, elev)]

    def to_cross_section(self, name: str) -> CrossSection:
        return CrossSection(
            name=name, points=self.points,
            bank_left=self.bank_left, bank_right=self.bank_right,
            n_main=self.n_main, n_left=self.n_left, n_right=self.n_right,
            reach_lengths=self.reach_lengths,
            c_contract=self.c_contract, c_expand=self.c_expand,
            levees=self.levees,
        )


@dataclass
class G01Reach:
    river: str
    reach: str
    sections: list = field(default_factory=list)   # [G01Section] US->DS

    @property
    def full_name(self):
        return f"{self.river}/{self.reach}"


@dataclass
class G01Junction:
    name: str
    inflows: list          # [(river, reach)] 按文件顺序
    outflow: tuple         # (river, reach)
    lengths: list          # [float] 与 inflows 顺序对应, 已换算


def parse_g01(path: str, unit: str = "ft"):
    """解析 .g01, 返回 (reaches, junctions)。

    reaches   : [G01Reach]
    junctions : [G01Junction]
    """
    L = _FT if unit == "ft" else 1.0
    fname = os.path.basename(path)

    reaches: dict = {}       # (river, reach) -> G01Reach (按出现顺序)
    reach_order = []
    junctions = []

    cur_reach = None
    cur_node = None          # 正在读取的 G01Section
    se_count = se_read = 0
    mann_count = mann_read = 0
    mann_buf = []
    se_buf = []
    junc_in_up = []          # 当前 junction 的入流 (river, reach)
    junc_len_left = 0        # 待配对的 Junc L&A 数量

    def _nums(line):
        out = []
        for tok in line.split():
            try:
                out.append(float(tok))
            except ValueError:
                pass
        return out

    def _nums8(line):
        """HEC-RAS 数组行为 8 字符定宽字段; 相邻数字可能粘连
        (如 '43.5611004.059'), 按定宽切分, 失败时回退空白切分。"""
        chunks = [line[i:i + 8] for i in range(0, len(line), 8)]
        try:
            return [float(c) for c in chunks if c.strip()]
        except ValueError:
            return _nums(line)

    def _flush_node():
        nonlocal cur_node, se_count, se_read, mann_count, mann_read
        nonlocal mann_buf, se_buf
        if cur_node is None:
            return
        pts = [(se_buf[i] * L, se_buf[i + 1] * L)
               for i in range(0, len(se_buf) - 1, 2)]
        triples = [(mann_buf[i], mann_buf[i + 1])
                   for i in range(0, len(mann_buf) - 2, 3)]
        nzones = [(s * L, n) for s, n in triples if 0.0 < n < 1e30]
        bl, br = cur_node.bank_left, cur_node.bank_right
        n_main_c = [(s, n) for s, n in nzones if bl <= s < br]
        n_left = [(s, n) for s, n in nzones if s < bl]
        n_right = [(s, n) for s, n in nzones if s >= br]
        if not n_main_c:
            raise ValueError(f"{fname} {cur_reach.full_name}@{cur_node.rs}: "
                             "主槽无 n 值")
        cur_node.points = pts
        cur_node.n_main = n_main_c[0][1]
        cur_node.n_left = n_left
        cur_node.n_right = n_right
        cur_reach.sections.append(cur_node)
        cur_node = None
        se_count = se_read = mann_count = mann_read = 0
        se_buf, mann_buf = [], []

    with open(path, encoding="latin-1") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if line.startswith("River Reach="):
                _flush_node()
                parts = [p.strip() for p in line.split("=", 1)[1].split(",")]
                key = (parts[0], parts[1])
                cur_reach = G01Reach(parts[0], parts[1])
                reaches[key] = cur_reach
                reach_order.append(key)
                continue

            if line.startswith("Type RM Length L Ch R ="):
                _flush_node()
                parts = [p.strip() for p in line.split("=", 1)[1].split(",")]
                ntype = int(float(parts[0]))
                if ntype != 1:
                    import warnings
                    warnings.warn(f"{fname}: 跳过非断面节点 (type={ntype}, "
                                  f"RS={parts[1]}) — v1 不支持桥涵等建筑物")
                    cur_node = None
                    continue
                lengths = tuple(
                    (float(v) if v.strip() else 0.0) * L for v in parts[2:5])
                cur_node = G01Section(
                    rs=parts[1], points=[], bank_left=0.0, bank_right=0.0,
                    n_main=0.0, n_left=[], n_right=[],
                    reach_lengths=lengths, c_contract=0.1, c_expand=0.3)
                continue

            if cur_node is not None and line.startswith("#Sta/Elev="):
                se_count = int(float(line.split("=")[1].split(",")[0]))
                continue
            if cur_node is not None and line.startswith("#Mann="):
                mann_count = int(float(line.split("=")[1].split(",")[0]))
                continue
            if cur_node is not None and line.startswith("Bank Sta="):
                p = line.split("=", 1)[1].split(",")
                cur_node.bank_left = float(p[0]) * L
                cur_node.bank_right = float(p[1]) * L
                continue
            if cur_node is not None and line.startswith("Levee="):
                # 格式: LF,Sta,Elev,RF,Sta,Elev (flag 0=无堤; -1/1=有堤)
                p = [t.strip() for t in line.split("=")[1].split(",")]
                for flag_i, sta_i, elev_i, side in ((0, 1, 2, -1),
                                                    (3, 4, 5, 1)):
                    if flag_i < len(p) and p[flag_i] not in ("", "0"):
                        try:
                            cur_node.levees.append(
                                (side, float(p[sta_i]) * L,
                                 float(p[elev_i]) * L))
                        except (ValueError, IndexError):
                            pass
                continue
            if cur_node is not None and line.startswith("Exp/Cntr="):
                p = line.split("=", 1)[1].split(",")
                cur_node.c_expand = float(p[0])
                cur_node.c_contract = float(p[1])
                continue

            if se_count and se_read < se_count * 2:
                vals = _nums8(line)
                take = min(len(vals), se_count * 2 - se_read)
                se_buf.extend(vals[:take])
                se_read += take
                if se_read >= se_count * 2:
                    se_count = 0
                continue
            if mann_count and mann_read < mann_count * 3:
                vals = _nums8(line)
                take = min(len(vals), mann_count * 3 - mann_read)
                mann_buf.extend(vals[:take])
                mann_read += take
                if mann_read >= mann_count * 3:
                    mann_count = 0
                continue

            if line.startswith("Junct Name="):
                _flush_node()
                junctions.append(G01Junction(
                    name=line.split("=", 1)[1].strip(),
                    inflows=[], outflow=None, lengths=[]))
                continue
            if junctions and line.startswith("Up River,Reach="):
                p = [x.strip() for x in line.split("=", 1)[1].split(",")]
                junctions[-1].inflows.append((p[0], p[1]))
                continue
            if junctions and line.startswith("Dn River,Reach="):
                p = [x.strip() for x in line.split("=", 1)[1].split(",")]
                junctions[-1].outflow = (p[0], p[1])
                continue
            if junctions and line.startswith("Junc L&A="):
                p = [x.strip() for x in line.split("=", 1)[1].split(",")]
                # 按文件顺序对应各入流河段; 值 = 汇合口长度 (每区同值)
                junctions[-1].lengths.append(float(p[0]) * L)
                continue

        _flush_node()

    # 汇合口长度与入流配对
    for j in junctions:
        if len(j.lengths) != len(j.inflows):
            raise ValueError(f"Junction {j.name}: Junc L&A 数 ({len(j.lengths)})"
                             f" 与入流河段数 ({len(j.inflows)}) 不一致")

    return [reaches[k] for k in reach_order], junctions


def build_network(g01_reaches, g01_junctions, outlet, boundary,
                  sbar_method="avg_conveyance", tol=0.003, max_iter=40):
    """把 parse_g01 的结果转成 SteadyNetwork。

    outlet      : (river, reach) 最下游河段
    boundary    : 出流河段下游边界元组
    各河段 chainage 自上游累计 (主槽河段长度)。
    流量: Reach.q 默认 0 — 由调用方按 .f01 或实际资料设置。
    """
    from .network import SteadyNetwork as _Net

    reaches = {}
    for gr in g01_reaches:
        cs_list = []
        c = 0.0
        for s in gr.sections:
            cs = s.to_cross_section(name=f"{gr.full_name}@{s.rs}")
            cs.chainage = c
            c += s.reach_lengths[1]
            cs_list.append(cs)
        reaches[gr.full_name] = Reach(gr.full_name, cs_list, q=0.0)

    junctions = []
    for gj in g01_junctions:
        junctions.append(Junction(
            name=gj.name,
            inflows={f"{r}/{rch}": L for (r, rch), L
                     in zip(gj.inflows, gj.lengths)},
            outflow=f"{gj.outflow[0]}/{gj.outflow[1]}"))

    outlet_name = f"{outlet[0]}/{outlet[1]}"
    return _Net(reaches, junctions, outlet_name, boundary,
                sbar_method=sbar_method, tol=tol, max_iter=max_iter)

"""HEC-RAS 计划结果 HDF (.p01.hdf) 读取器。

读取几何 (断面属性、站距-高程、n 值分区) 与稳态计算结果
(Water Surface / Energy Grade / 临界水位 / Alpha / 分区输水能力等),
转换为 raslite 的 CrossSection (SI 单位) 与对比数据。

单位换算: HEC-RAS 美制 (ft, cfs) -> SI (m, m3/s)
  长度 1 ft = 0.3048 m;  1 cfs = 0.0283168466 m3/s
  输水能力 K 与流量同量纲: K_si = K_eng * 0.0283168466
  面积: *0.3048^2

已知局限: 仅支持 'Basic Mann n' 摩阻模式; 忽略无效流动区/障碍物
(遇到含这些特性的断面会给出警告)。
"""

from __future__ import annotations

from dataclasses import dataclass, field

FT = 0.3048                # m / ft
CFS = 0.0283168466         # m3/s / cfs
FT2 = FT * FT              # m2 / ft2


def _s(b) -> str:
    return b.decode("latin-1").strip() if isinstance(b, bytes) else str(b).strip()


@dataclass
class HecrasSection:
    """单个断面的几何与结果数据 (已换算 SI)。"""
    river: str
    reach: str
    rs: str                    # 河流站距标识 (如 '9.79')
    points: list               # [(station_m, elev_m)]
    bank_left: float
    bank_right: float
    n_main: float
    n_left: list               # [(station_m, n), ...] 或 float
    n_right: list
    reach_lengths: tuple       # (L_lob, L_ch, L_rob) m, 指向下一断面
    c_contract: float
    c_expand: float
    # 结果 (按剖面名索引, None 表示无)
    results: dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        return f"{self.river}/{self.reach}@{self.rs}"


def read_plan_hdf(path: str, to_si: bool = True):
    """读取 .p01.hdf, 返回 (sections, profiles, meta)。

    sections : [HecrasSection] 按文件内顺序 (US->DS)
    profiles : [str] 剖面名列表 (按列顺序)
    meta     : dict: units, version 等原始属性
    """
    import h5py
    import numpy as np

    f = h5py.File(path, "r")
    meta = {k: _s(v) if isinstance(v, bytes) else v for k, v in f.attrs.items()}
    us = _s(f.attrs.get("Units System", b"")) if "Units System" in f.attrs else ""
    if "US Customary" not in us:
        # 公制工程文件: 不换算
        to_si_eff = False
    else:
        to_si_eff = to_si
    L = FT if to_si_eff else 1.0
    Q = CFS if to_si_eff else 1.0
    A = FT2 if to_si_eff else 1.0

    cs = f["Geometry/Cross Sections"]
    attrs = cs["Attributes"][:]
    se_info = cs["Station Elevation Info"][:]
    se_vals = cs["Station Elevation Values"][:]
    n_info = cs["Manning's n Info"][:]
    n_vals = cs["Manning's n Values"][:]

    sections = []
    for i, row in enumerate(attrs):
        river, reach, rs = _s(row[0]), _s(row[1]), _s(row[2])
        off, cnt = int(se_info[i][0]), int(se_info[i][1])
        pts = [(float(x) * L, float(z) * L) for x, z in se_vals[off:off + cnt]]
        noff, ncnt = int(n_info[i][0]), int(n_info[i][1])
        # 过滤哨兵值: HEC-RAS 用 float32 max 表示未定义 n (滩地无定义时),
        # n<=0 亦为占位垃圾 (如 Bank Sta 与断面端重合时的 0 值条目)
        nzones = [(float(s) * L, float(n)) for s, n in n_vals[noff:noff + ncnt]
                  if 0.0 < float(n) < 1e30]
        bl, br = float(row[8]) * L, float(row[9]) * L

        # n 值按滩肩拆分: <bl -> 左滩; bl<=x<br -> 主槽; >=br -> 右滩
        n_main_cands = [(s, n) for s, n in nzones if bl <= s < br]
        n_left = [(s, n) for s, n in nzones if s < bl]
        n_right = [(s, n) for s, n in nzones if s >= br]
        if not n_main_cands:
            raise ValueError(f"{river}/{reach}@{rs}: 主槽无 n 值定义")
        if len(n_main_cands) > 1:
            import warnings
            warnings.warn(f"{river}/{reach}@{rs}: 主槽有多个 n 值, 取第一个")
        n_main = n_main_cands[0][1]
        # HEC-RAS 语义: 滩地列表首元素给出首区 n (起点即断面端)
        # 保留原列表 (geometry._norm_n_zones 按 HEC-RAS 语义解析)

        friction_mode = _s(row[10])
        if friction_mode not in ("Basic Mann n", "Horiz Mann n"):
            raise NotImplementedError(f"暂不支持摩阻模式: {friction_mode}")

        sections.append(HecrasSection(
            river=river, reach=reach, rs=rs, points=pts,
            bank_left=bl, bank_right=br,
            n_main=n_main, n_left=n_left or 0.0, n_right=n_right or 0.0,
            # Attributes 字段序: 5=Len Left, 6=Len Channel, 7=Len Right,
            #                    11=Contr, 12=Expan (13=HP Count 勿用!)
            reach_lengths=(_nz(float(row[5])) * L, _nz(float(row[6])) * L,
                           _nz(float(row[7])) * L),
            c_contract=float(row[11]), c_expand=float(row[12]),
            results={},
        ))

    # ---- 结果 ----
    g = f["Results/Steady/Output/Output Blocks/Base Output/Steady Profiles"]
    xsg = g["Cross Sections"]
    ws = np.array(xsg["Water Surface"])          # (nprof, nxs)
    eg = np.array(xsg["Energy Grade"])
    q = np.array(xsg["Flow"])
    add = xsg["Additional Variables"]

    prof_attrs = {}
    for key in g.attrs:
        prof_attrs[key] = g.attrs[key]
    profiles = None
    # 剖面名通常在 Steady Profiles 组属性或 Cross Sections 组属性
    for src_ in (g, xsg, xsg["Water Surface"]):
        for k in src_.attrs:
            if "profile" in k.lower() and profiles is None:
                profiles = [_s(v) for v in list(src_.attrs[k])]
    if profiles is None:
        profiles = [f"P{i+1}" for i in range(ws.shape[0])]

    extra = {
        "alpha": np.array(add["Alpha"]),
        "K_total": np.array(add["Conveyance Total"]),
        "K_ch": np.array(add["Conveyance Channel"]),
        "K_lob": np.array(add["Conveyance Left OB"]),
        "K_rob": np.array(add["Conveyance Right OB"]),
        "crit_ws": np.array(add["Critical Water Surface"]),
        "area_total": np.array(add["Area Flow Total"]),
        "area_ch": np.array(add["Area Flow Channel"]),
        "topw_ch": np.array(add["Top Width Channel"]),
        "sf": np.array(add["Friction Slope"]),
        "vel_ch": np.array(add["Vel Channel"]) if "Vel Channel" in add else None,
    }

    for j in range(ws.shape[0]):
        pname = profiles[j]
        for i, sec in enumerate(sections):
            sec.results[pname] = {
                "ws": float(ws[j, i]) * L,
                "eg": float(eg[j, i]) * L,
                "q": float(q[j, i]) * Q,
                "alpha": float(extra["alpha"][j, i]),
                "K_total": float(extra["K_total"][j, i]) * Q,
                "K_ch": float(extra["K_ch"][j, i]) * Q,
                "K_lob": float(extra["K_lob"][j, i]) * Q,
                "K_rob": float(extra["K_rob"][j, i]) * Q,
                "crit_ws": float(extra["crit_ws"][j, i]) * L,
                "area_total": float(extra["area_total"][j, i]) * A,
                "area_ch": float(extra["area_ch"][j, i]) * A,
                "sf": float(extra["sf"][j, i]),
            }

    # ---- 河段中心线 (junction 距离, 单位: 英里) ----
    try:
        cl = []
        for row in f["Geometry/River Centerlines/Attributes"][:]:
            cl.append({
                "river": _s(row[0]), "reach": _s(row[1]),
                "us_type": _s(row[2]), "us_name": _s(row[3]),
                "ds_type": _s(row[4]), "ds_name": _s(row[5]),
                # DS 断面到 junction / junction 到 US 断面 (英里)
                "ds_xs_to_junction_mi": float(row[6]),
                "junc_to_us_xs_mi": float(row[7]),
            })
        meta["centerlines"] = cl
    except KeyError:
        pass

    f.close()
    return sections, profiles, meta


def _nz(v: float) -> float:
    """NaN (无滩地时的河段长度) -> 0"""
    return 0.0 if v != v else v


def read_flow_file(path: str, to_si: bool = True):
    """解析 HEC-RAS 稳定流流量文件 (.f01)。

    返回 (profiles, flows, boundaries):
      profiles   : [str] 剖面名
      flows      : {(river, reach, rs): [q, ...]} 各位置逐剖面流量 (SI)
      boundaries : {(river, reach): {'up_type', 'dn_type', 'dn_slope', ...}}
    """
    import re
    Q = CFS if to_si else 1.0
    profiles = []
    flows = {}
    boundaries = {}
    cur_flow_loc = None
    cur_bnd = None
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Profile Names="):
                profiles = [p.strip() for p in line.split("=", 1)[1].split(",")]
            elif line.startswith("River Rch & RM="):
                parts = [p.strip() for p in line.split("=", 1)[1].split(",")]
                cur_flow_loc = (parts[0], parts[1], parts[2])
                flows[cur_flow_loc] = []
            elif cur_flow_loc is not None and re.match(r"^\s*[\d.eE+-]+(\s+[\d.eE+-]+)*\s*$", line):
                flows[cur_flow_loc].extend(
                    float(v) * Q for v in line.split())
            elif line.startswith("Boundary for River Rch & Prof#="):
                parts = [p.strip() for p in line.split("=", 1)[1].split(",")]
                key = (parts[0], parts[1])
                cur_bnd = boundaries.setdefault(key, {})
            elif line.startswith("Up Type=") and cur_bnd is not None:
                cur_bnd["up_type"] = int(line.split("=")[1])
            elif line.startswith("Up Slope=") and cur_bnd is not None:
                cur_bnd["up_slope"] = float(line.split("=")[1])
            elif line.startswith("Dn Type=") and cur_bnd is not None:
                cur_bnd["dn_type"] = int(line.split("=")[1])
            elif line.startswith("Dn Slope=") and cur_bnd is not None:
                cur_bnd["dn_slope"] = float(line.split("=")[1])
    # 去掉非流量的杂散行 (如 Profile Names 行残留)
    return profiles, flows, boundaries

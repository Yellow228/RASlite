"""test04 (YXG) 算例: g01 导入 -> raslite 计算 -> 与 HEC-RAS 6.2 结果对比。

真实河道算例: River 1/Reach 1, SI 单位, 6 个设计频率
(1%~50%, Q=95.5~442 m3/s), 逐剖面下游已知水位。
注意: g01 含桥涵节点 (RS 1800 / RS 70), v1 跳过, 相应壅水缺失。
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import SteadyNetwork
from raslite.g01 import parse_g01, build_network, read_prj_units
from raslite.hecras_hdf import read_plan_hdf, read_flow_file

BASE = os.path.join(HERE, "..", "examples", "test04")
G01 = os.path.join(BASE, "YXG.g01")
PRJ = os.path.join(BASE, "YXG.prj")
HDF = os.path.join(BASE, "YXG.p01.hdf")
F01 = os.path.join(BASE, "YXG.f01")

RNAME = "River 1/Reach 1"


def main():
    # 1. 几何: g01 导入 (单位自动识别)
    unit = read_prj_units(PRJ)
    g01_reaches, g01_junctions = parse_g01(G01, unit=unit)
    net = build_network(g01_reaches, g01_junctions,
                        outlet=("River 1", "Reach 1"), boundary=None)

    # 2. 流量与边界: 来自 .f01 (逐剖面)
    f_profiles, flows, bnds = read_flow_file(F01, to_si=(unit == "ft"))
    loc = next(iter(flows))
    q_list = flows[loc]
    bnd_list = bnds[(loc[0], loc[1])]

    # 3. HEC-RAS 结果 (对比基准)
    sections, hdf_profiles, meta = read_plan_hdf(HDF)
    ref = {f"{s.river}/{s.reach}@{s.rs}": s.results for s in sections}

    lines = []
    lines.append("# test04 (YXG): raslite vs HEC-RAS 6.2\n")
    n_sec = len(net.reaches[RNAME].sections)
    lines.append(f"- 几何: YXG.g01 导入 (单位 {unit}), 单河段 {n_sec} 断面")
    lines.append("- 6 个设计频率 (1%~50%), 流量与下游已知水位来自 YXG.f01")
    lines.append("- 注意: g01 含桥涵节点 (RS 1800 / RS 70), v1 跳过, "
                 "其上游壅水缺失\n")
    lines.append("| 频率 | RS | Q(cms) | HEC WS (m) | raslite WS (m) | Δ (cm) | 流态(HEC/ras) |")
    lines.append("|---|---|---|---|---|---|---|")

    pmax = 0.0
    for ip, pname in enumerate(f_profiles):
        net.reaches[RNAME].q = q_list[ip]
        net.boundary = ("ws", bnd_list[ip]["dn_known_ws"])
        net.results = {}
        net._run_reach(RNAME, net.boundary)
        for r in net.results[RNAME]:
            refd = ref.get(r.section.name, {}).get(hdf_profiles[ip])
            if refd is None:
                continue
            dws = (r.ws - refd["ws"]) * 100.0
            pmax = max(pmax, abs(dws))
            hec_zone = "临界" if abs(refd["ws"] - refd["crit_ws"]) < 0.005 \
                else ("超" if refd["ws"] < refd["crit_ws"] else "亚")
            ras_zone = "临界" if r.defaulted_critical else (
                "超" if r.ws < r.crit_ws else "亚")
            lines.append(f"| {pname} | {r.section.name.split('@')[1]} | "
                         f"{r.q:.1f} | {refd['ws']:.4f} | {r.ws:.4f} | "
                         f"{dws:+.2f} | {hec_zone}/{ras_zone} |")
    lines.append("")
    # 统计: 排除 HEC 默认临界断面后的最大误差
    pmax2 = 0.0
    for ip, pname in enumerate(f_profiles):
        net.reaches[RNAME].q = q_list[ip]
        net.boundary = ("ws", bnd_list[ip]["dn_known_ws"])
        net.results = {}
        net._run_reach(RNAME, net.boundary)
        for r in net.results[RNAME]:
            refd = ref.get(r.section.name, {}).get(hdf_profiles[ip])
            if refd is None:
                continue
            if abs(refd["ws"] - refd["crit_ws"]) < 0.005:
                continue          # HEC 默认临界的断面不计
            pmax2 = max(pmax2, abs((r.ws - refd["ws"]) * 100.0))
    lines.append(f"**最大 |ΔWS| = {pmax:.2f} cm**")
    lines.append(f"**排除 HEC 默认临界断面后 = {pmax2:.2f} cm**\n")
    lines.append("各断面警告:")
    any_warn = False
    for r in net.results[RNAME]:
        for w in r.warnings:
            any_warn = True
            lines.append(f"- {r.section.name}: {w}")
    if not any_warn:
        lines.append("- 无")

    text = "\n".join(lines)
    print(text)
    out = os.path.join(HERE, "..", "L2_验证报告_test04.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\n已保存: {out}")
    return pmax


if __name__ == "__main__":
    main()

"""test05 (XHYXC) 算例: g01 导入 -> raslite 计算 -> 与 HEC-RAS 6.2 结果对比。

真实河道算例: River 1/Reach 1, 10 断面, Q=333 m3/s (RS 897.0463 入流),
下游边界已知水位 1269.75 m, SI 单位制。
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import SteadyNetwork
from raslite.g01 import parse_g01, build_network, read_prj_units
from raslite.hecras_hdf import read_plan_hdf, read_flow_file

BASE = os.path.join(HERE, "..", "examples", "test05")
G01 = os.path.join(BASE, "XHYXC.g01")
PRJ = os.path.join(BASE, "XHYXC.prj")
HDF = os.path.join(BASE, "XHYXC.p01.hdf")
F01 = os.path.join(BASE, "XHYXC.f01")
RNAME = "River 1/Reach 1"


def main():
    # 1. 几何: g01 导入 (单位自动识别 SI)
    unit = read_prj_units(PRJ)
    g01_reaches, g01_junctions = parse_g01(G01, unit=unit)
    net = build_network(g01_reaches, g01_junctions,
                        outlet=("River 1", "Reach 1"),
                        boundary=("ws", 1269.75))

    # 2. 流量: .f01 (RS 897.0463 = 333 m3/s, 沿程不变)
    f_profiles, flows, bnds = read_flow_file(F01, to_si=(unit == "ft"))
    loc = next(iter(flows))
    q = flows[loc]
    net.reaches[RNAME].q = q[0]  # 沿程不变

    # 3. 计算
    results = net.run_subcritical()[RNAME]

    # 4. 对比 HEC-RAS 结果
    sections, hdf_profiles, meta = read_plan_hdf(HDF)
    ref = {f"{s.river}/{s.reach}@{s.rs}": s.results[hdf_profiles[0]]
           for s in sections}

    lines = []
    lines.append("# test05 (XHYXC): raslite vs HEC-RAS 6.2\n")
    lines.append(f"- 几何: XHYXC.g01 导入 (单位 {unit}), 单河段 10 断面")
    lines.append(f"- 流量: {q[0]:.2f} m3/s (RS {loc[2]});  "
                 f"下游边界: 已知水位 1269.75 m")
    lines.append("- 计划: 亚临界, Std Step Tol=0.003, Friction Slope Method=1\n")
    lines.append("| RS | Q(cms) | HEC WS (m) | raslite WS (m) | Δ (cm) | 流态(HEC/ras) |")
    lines.append("|---|---|---|---|---|---|")
    pmax = 0.0
    for r in results:            # DS->US 推算顺序
        key = r.section.name
        refd = ref[key]
        dws = (r.ws - refd["ws"]) * 100.0
        pmax = max(pmax, abs(dws))
        hec_zone = "临界" if abs(refd["ws"] - refd["crit_ws"]) < 0.005 \
            else ("超" if refd["ws"] < refd["crit_ws"] else "亚")
        ras_zone = "临界" if r.defaulted_critical else (
            "超" if r.ws < r.crit_ws else "亚")
        lines.append(f"| {key.split('@')[1]} | {r.q:.2f} | "
                     f"{refd['ws']:.4f} | {r.ws:.4f} | {dws:+.2f} | "
                     f"{hec_zone}/{ras_zone} |")
    lines.append("")
    lines.append(f"**最大 |ΔWS| = {pmax:.2f} cm**\n")
    lines.append("各断面警告:")
    any_warn = False
    for r in results:
        for w in r.warnings:
            any_warn = True
            lines.append(f"- {r.section.name}: {w}")
    if not any_warn:
        lines.append("- 无")

    text = "\n".join(lines)
    print(text)
    out = os.path.join(HERE, "..", "L2_验证报告_test05.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\n已保存: {out}")
    return pmax


if __name__ == "__main__":
    main()

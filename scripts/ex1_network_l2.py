"""L2 验证: 河网 + Junction — raslite vs HEC-RAS, EX1 全模型。

河网结构 (来自 EX1):
  Butte Cr. (Tributary, 3 断面, Q=100cfs)  ──┐
  Fall River (Upper Reach, 3 断面, Q=500) ───┤ Junction "Sutter"
                                             ↓
  Fall River (Lower Reach, 4 断面, Q=600→650)

汇合口长度 (来自 g01 Junc L&A / r01 Junction Information):
  Upper Reach 60ft, Tributary 50ft。
下游边界: Fall River Lower Reach 末断面正常水深 S0=0.0004。
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import SteadyNetwork
from raslite.g01 import parse_g01, build_network, read_prj_units
from raslite.hecras_hdf import read_plan_hdf, read_flow_file

BASE = os.path.join(HERE, "..", "examples", "1D Steady Flow Hydraulics",
                    "Chapter 4 Example Data")
HDF = os.path.join(BASE, "EX1.p01.hdf")
F01 = os.path.join(BASE, "EX1.f01")
G01 = os.path.join(BASE, "EX1.g01")
PRJ = os.path.join(BASE, "EX1.prj")


def build(profile_flows):
    """profile_flows: {reach_full_name: [q,...]} 已是 SI (cms)。
    几何全部来自 .g01 导入器 (含汇合口长度), 结果对比仍用 HDF。"""
    unit = read_prj_units(PRJ)
    g01_reaches, g01_junctions = parse_g01(G01, unit=unit)
    net = build_network(g01_reaches, g01_junctions,
                        outlet=("Fall River", "Lower Reach"),
                        boundary=("normal", 0.0004))
    for name, qs in profile_flows.items():
        net.reaches[name].q = qs

    sections, hdf_profiles, meta = read_plan_hdf(HDF)
    ref = {}
    for s in sections:
        ref[(s.river, s.reach, s.rs)] = s.results
    upper = net.reaches["Fall River/Upper Reach"]
    butte = net.reaches["Butte Cr./Tributary"]
    lower = net.reaches["Fall River/Lower Reach"]
    return net, ref, hdf_profiles, upper, butte, lower


def main():
    hdf_profiles = read_plan_hdf(HDF)[1]
    f_profiles, flows, bnds = read_flow_file(F01)

    report = []
    report.append("# L2 验证: 河网 + Junction — raslite vs HEC-RAS (EX1 全模型)\n")
    report.append("- 河网: Butte Cr.(3断面) + Fall River Upper(3) → Junction Sutter")
    report.append("  → Fall River Lower(4);几何与汇合口长度均来自 EX1.g01 导入器")
    report.append("- 流量: 每剖面按 EX1.f01 重建 (Upper 500/2000/5000, Butte 100/500/1500,")
    report.append("  Lower 600/2500/6500 → 650/2700/7000 cfs)")
    report.append("- 边界: Lower Reach 末断面正常水深 S0=0.0004;单位 SI\n")

    grand_max = 0.0
    for ip, pname in enumerate(f_profiles):
        hname = hdf_profiles[ip]
        qmap = {
            "Fall River/Upper Reach": flows[("Fall River", "Upper Reach", "10")][ip],
            "Butte Cr./Tributary": flows[("Butte Cr.", "Tributary", "0.2")][ip],
            "Fall River/Lower Reach": [
                flows[("Fall River", "Lower Reach", "9.79")][ip],
                flows[("Fall River", "Lower Reach", "9.79")][ip],
                flows[("Fall River", "Lower Reach", "9.6")][ip],
                flows[("Fall River", "Lower Reach", "9.6")][ip],
            ],
        }
        net, ref, _, upper, butte, lower = build(qmap)
        results = net.run_subcritical()

        report.append(f"## 剖面 {pname}\n")
        report.append("| 河段 | RS | Q(cms) | HEC WS (m) | raslite WS (m) | Δ (cm) |")
        report.append("|---|---|---|---|---|---|")
        pmax = 0.0
        for rname, secs in [("Butte Cr./Tributary", butte.sections),
                            ("Fall River/Upper Reach", upper.sections),
                            ("Fall River/Lower Reach", lower.sections)]:
            for r, cs in zip(results[rname], secs):
                key = None
                for (riv, rch, rs), resd in ref.items():
                    if f"{riv}/{rch}@{rs}" == cs.name:
                        key = (riv, rch, rs)
                        break
                refd = ref[key][hname]
                dws = (r.ws - refd["ws"]) * 100.0
                pmax = max(pmax, abs(dws))
                report.append(
                    f"| {rname.split('/')[1]} | {key[2]} | {r.q:.2f} | "
                    f"{refd['ws']:.4f} | {r.ws:.4f} | {dws:+.2f} |")
        report.append("")
        report.append(f"**{pname} 最大 |ΔWS| = {pmax:.2f} cm**\n")
        grand_max = max(grand_max, pmax)

    report.append("\n## 总结\n")
    report.append(f"- 全部剖面最大 |ΔWS| = **{grand_max:.2f} cm**")
    report.append("- 验收目标: ≤ 1.0 cm")

    text = "\n".join(report)
    print(text)
    out = os.path.join(HERE, "..", "L2_验证报告_EX1河网.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\n已保存: {out}")
    return grand_max


if __name__ == "__main__":
    main()

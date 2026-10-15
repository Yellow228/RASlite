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

from raslite import CrossSection, Reach, Junction, SteadyNetwork
from raslite.hecras_hdf import read_plan_hdf, read_flow_file, CFS

BASE = os.path.join(HERE, "..", "examples", "1D Steady Flow Hydraulics",
                    "Chapter 4 Example Data")
HDF = os.path.join(BASE, "EX1.p01.hdf")
F01 = os.path.join(BASE, "EX1.f01")

FT = 0.3048


def build(profile_flows):
    """profile_flows: (q_upper, q_butte, q_lower_at_979, q_lower_at_96) cfs"""
    qu, qb, ql1, ql2 = profile_flows  # 已是 SI (cms)
    sections, hdf_profiles, meta = read_plan_hdf(HDF)
    f_profiles, flows, bnds = read_flow_file(F01)

    def mk(sec):
        cs = CrossSection(
            name=f"{sec.river}/{sec.reach}@{sec.rs}",
            points=sec.points,
            bank_left=sec.bank_left, bank_right=sec.bank_right,
            n_main=sec.n_main, n_left=sec.n_left, n_right=sec.n_right,
            reach_lengths=sec.reach_lengths,
            c_contract=sec.c_contract, c_expand=sec.c_expand,
        )
        return cs

    def reach_of(river, reach_name, q):
        ss = sorted([s for s in sections
                     if s.river == river and s.reach == reach_name],
                    key=lambda s: -float(s.rs))
        cs_list = [mk(s) for s in ss]
        # chainage: 自上游向下游累计 (沿流向递增)
        c = 0.0
        for cs, s in zip(cs_list, ss):
            cs.chainage = c
            c += s.reach_lengths[1]
        return Reach(f"{river}/{reach_name}", cs_list, q=q)

    upper = reach_of("Fall River", "Upper Reach", qu)
    butte = reach_of("Butte Cr.", "Tributary", qb)
    lower = reach_of("Fall River", "Lower Reach", [ql1, ql1, ql2, ql2])

    junc = Junction("Sutter",
                    inflows={"Fall River/Upper Reach": 60 * FT,
                             "Butte Cr./Tributary": 50 * FT},
                    outflow="Fall River/Lower Reach")

    net = SteadyNetwork(
        reaches={"Fall River/Upper Reach": upper,
                 "Butte Cr./Tributary": butte,
                 "Fall River/Lower Reach": lower},
        junctions=[junc],
        outlet_reach="Fall River/Lower Reach",
        boundary=("normal", bnds[("Fall River", "Lower Reach")]["dn_slope"]),
    )
    # HDF 结果映射: {(river, reach, rs): results}
    ref = {}
    for s in sections:
        ref[(s.river, s.reach, s.rs)] = s.results
    return net, ref, f_profiles, hdf_profiles, upper, butte, lower


def main():
    sections0, hdf_profiles, meta = read_plan_hdf(HDF)
    f_profiles, flows, bnds = read_flow_file(F01)

    report = []
    report.append("# L2 验证: 河网 + Junction — raslite vs HEC-RAS (EX1 全模型)\n")
    report.append("- 河网: Butte Cr.(3断面) + Fall River Upper(3) → Junction Sutter")
    report.append("  → Fall River Lower(4);汇合口长度 Upper 60ft / Butte 50ft")
    report.append("- 流量: 每剖面按 EX1.f01 重建 (Upper 500/2000/5000, Butte 100/500/1500,")
    report.append("  Lower 600/2500/6500 → 650/2700/7000 cfs)")
    report.append("- 边界: Lower Reach 末断面正常水深 S0=0.0004;单位 SI\n")

    grand_max = 0.0
    for ip, pname in enumerate(f_profiles):
        hname = hdf_profiles[ip]
        pf = (flows[("Fall River", "Upper Reach", "10")][ip],
              flows[("Butte Cr.", "Tributary", "0.2")][ip],
              flows[("Fall River", "Lower Reach", "9.79")][ip],
              flows[("Fall River", "Lower Reach", "9.6")][ip])
        net, ref, _, _, upper, butte, lower = build(pf)
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

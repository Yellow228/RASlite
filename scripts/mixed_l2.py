"""L2 验证: 混合流 (Mixed Flow Regime) — raslite vs HEC-RAS, MIXED 官方算例。

模型: 单一矩形渠道 (20ft 宽, n=0.015), 19 断面, 陡坡 (0.01) 接缓坡
(下游边界正常水深 0.000586), 两剖面 (500 / 1000 cfs)。
HEC-RAS 计划开启 Mixed Flow 模式; raslite 用
regime="mixed" (亚临界遍 + 超临界遍 + 比力合并) 复现同款算法。
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import CrossSection, SteadyProfile
from raslite.hecras_hdf import read_plan_hdf, read_flow_file, CFS

BASE = os.path.join(HERE, "..", "examples", "1D Steady Flow Hydraulics",
                    "Mixed Flow Regime Channel")
HDF = os.path.join(BASE, "MIXED.p01.hdf")
F01 = os.path.join(BASE, "MIXED.f01")


def main():
    sections, hdf_profiles, meta = read_plan_hdf(HDF)
    f_profiles, flows, bnds = read_flow_file(F01)

    # 单河段: RS 0.5682 (US) -> 0 (DS)
    secs = sorted(sections, key=lambda s: -float(s.rs))
    xs_list = []
    chain = 0.0
    for s in secs:
        cs = CrossSection(
            name=f"XS@{s.rs}",
            points=s.points,
            bank_left=s.bank_left, bank_right=s.bank_right,
            n_main=s.n_main, n_left=s.n_left, n_right=s.n_right,
            reach_lengths=s.reach_lengths,
            c_contract=s.c_contract, c_expand=s.c_expand,
        )
        cs.chainage = chain
        chain += s.reach_lengths[1] if s.reach_lengths[1] > 0 else 30.48
        xs_list.append((cs, s))

    # 边界 (来自 .f01): 上游正常水深 S0=0.01, 下游正常水深 S0=0.000586
    bnds_fall = bnds[("Mixed Reach", "Mixed Reach")]
    ds_bc = ("normal", bnds_fall["dn_slope"])
    us_bc = ("normal", bnds_fall["up_slope"])

    report = []
    report.append("# L2 验证: 混合流 — raslite vs HEC-RAS (MIXED 官方算例)\n")
    report.append("- 模型: 矩形渠道 20ft 宽, n=0.015, 19 断面, 陡坡(0.01)接缓坡")
    report.append("- 算法: regime='mixed' — 亚临界遍(自下游) + 超临界遍(自上游)")
    report.append("  + 比力 SF 合并 (复刻 HEC-RAS Mixed Flow 流程)")
    report.append(f"- 边界: 上游正常水深 S0={bnds_fall['up_slope']}, "
                  f"下游正常水深 S0={bnds_fall['dn_slope']}\n")

    grand_max = 0.0
    for ip, pname in enumerate(f_profiles):
        hname = hdf_profiles[ip]
        q = flows[("Mixed Reach", "Mixed Reach", "0.5682")][ip]
        prof = SteadyProfile([cs for cs, _ in xs_list], q,
                             regime="mixed", boundary=(ds_bc, us_bc),
                             tol=1e-6)
        results = prof.run()          # US->DS
        by_name = {r.section.name: r for r in results}

        report.append(f"## 剖面 {pname} (Q = {q:.2f} cms)\n")
        report.append("| RS | HEC WS (m) | raslite WS (m) | Δ (cm) | "
                      "HEC Fr 区| ras 流态 | 备注 |")
        report.append("|---|---|---|---|---|---|---|")
        pmax = 0.0
        for cs, s in xs_list:
            r = by_name[cs.name]
            ref = s.results[hname]
            dws = (r.ws - ref["ws"]) * 100.0
            pmax = max(pmax, abs(dws))
            hec_zone = "超" if ref["ws"] < ref["crit_ws"] else "亚"
            ras_zone = "超" if r.ws < r.crit_ws else "亚"
            notes = "; ".join(r.warnings) if r.warnings else ""
            notes = notes[:40]
            report.append(
                f"| {s.rs} | {ref['ws']:.4f} | {r.ws:.4f} | {dws:+.2f} | "
                f"{hec_zone} | {ras_zone} | {notes} |")
        report.append("")
        report.append(f"**{pname} 最大 |ΔWS| = {pmax:.2f} cm**\n")
        grand_max = max(grand_max, pmax)

    report.append("\n## 总结\n")
    report.append(f"- 全部剖面最大 |ΔWS| = **{grand_max:.2f} cm**")
    report.append("- 验收目标: ≤ 1.0 cm")

    text = "\n".join(report)
    print(text)
    out = os.path.join(HERE, "..", "L2_验证报告_MIXED.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\n已保存: {out}")
    return grand_max


if __name__ == "__main__":
    main()

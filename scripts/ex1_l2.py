"""L2 验证: raslite vs HEC-RAS, EX1 官方算例 Fall River 单河段。

模型说明:
- 取 Fall River (Upper Reach 3 断面 + Lower Reach 4 断面, 共 7 断面),
  在支流汇合口 (Butte Cr.) 下游直接衔接 (junction 摩阻近似);
- 流量沿程变化: 上游 500 cfs, 汇合口后 600 cfs, 9.6 断面后 650 cfs
  (010 yr 剖面), 其余剖面同理, 按 HEC-RAS 输入逐断面给定;
- 下游边界: 正常水深 S0=0.0004 (来自 EX1 计划设置);
- 几何/流量由 .p01.hdf 读取并换算为 SI。
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import CrossSection, SteadyProfile
from raslite.hecras_hdf import read_plan_hdf, read_flow_file

BASE = os.path.join(HERE, "..", "examples", "1D Steady Flow Hydraulics",
                    "Chapter 4 Example Data")
HDF = os.path.join(BASE, "EX1.p01.hdf")
F01 = os.path.join(BASE, "EX1.f01")


def build_model():
    sections, profiles, meta = read_plan_hdf(HDF)
    # 只取 Fall River (上下游 reach 顺序拼接)
    fall = [s for s in sections if s.river == "Fall River"]
    fall.sort(key=lambda s: -float(s.rs))  # US -> DS (RS 大在上游)

    # chainage: 自上游向下游累计 (沿流向递增, m)。
    # 汇合口两侧断面 (RS 9.8 / 9.79) 的 HDF 河段长度为 0 或含义特殊,
    # 用中心线的 junction 距离 (英里 -> m)。
    MI = 1609.344
    cl = {c["reach"]: c for c in meta.get("centerlines", [])}
    step_of = {}
    for s in fall:
        step_of[s.rs] = s.reach_lengths[1]  # 默认: 本断面 -> 下一断面
    step_of["9.8"] = (cl["Upper Reach"]["ds_xs_to_junction_mi"]
                      + cl["Lower Reach"]["junc_to_us_xs_mi"]) * MI
    chains = [0.0]
    for s in fall[:-1]:
        chains.append(chains[-1] + step_of[s.rs])

    xs_list = []
    for s, c in zip(fall, chains):
        rl = s.reach_lengths
        if s.rs == "9.8":
            # 汇合口上游侧断面: HDF 存 0。HEC-RAS junction 能量法实际采用
            # Junction Information 中的汇合口长度 (r01: 60ft, 三区同值),
            # 而非中心线里程差 (0.07mi 会高估摩擦损失)。
            d = 60.0 * 0.3048
            rl = (d, d, d)
        cs = CrossSection(
            name=f"{s.river}@{s.rs}",
            points=s.points,
            bank_left=s.bank_left, bank_right=s.bank_right,
            n_main=s.n_main, n_left=s.n_left, n_right=s.n_right,
            reach_lengths=rl,
            c_contract=s.c_contract, c_expand=s.c_expand,
        )
        cs.chainage = c
        xs_list.append((cs, s))
    return xs_list, profiles, meta


def flow_profile(fall_sections, pname):
    """按 EX1 流量数据构造逐断面流量 (m3/s):
    010yr: 500@US..汇合口, 600@汇合口下, 650@RS9.6 以下
    050yr: 2000 / 2500 / 2700 ; 0100yr: 5000 / 6500 / 7000
    """
    base = {"010 yr": (500, 600, 650, 7),
            "050 yr": (2000, 2500, 2700, 7),
            "0100 yr": (5000, 6500, 7000, 9)}
    # 注: 节点序号 (1..10) 中 Fall River 为 4..10
    #   010yr: 500 从节点 4 起; 600 从节点 7 起 (汇合口); 650 从节点 9 起 (RS 9.6)
    key = None
    for k in base:
        if k.replace(" ", "") in pname.replace(" ", ""):
            key = k
    if key is None:
        raise ValueError(f"未知剖面: {pname}")
    q_us, q_junc, q_low, n_start_low = base[key]
    from raslite.hecras_hdf import CFS
    qs = []
    for s in fall_sections:
        rs = float(s.rs)
        if rs >= 9.8:          # Upper Reach (含 9.8, 汇合口上游侧)
            q = q_us
        elif rs >= 9.6:        # 9.79, 9.7 (汇合口下游)
            q = q_junc
        else:                  # 9.6 及以下 (9.6 处有沿程入流? 节点 9 起 650)
            q = q_low if rs < 9.6 else q_junc
        qs.append(q * CFS)
    # 修正: 650 从节点 9 (RS 9.6) 起 -> RS 9.6 本身也是 650
    for i, s in enumerate(fall_sections):
        if abs(float(s.rs) - 9.6) < 1e-9:
            qs[i] = q_low * CFS
    return qs


def main():
    xs_list, profiles, meta = build_model()
    fall_secs = [s for _, s in xs_list]

    # 从 .f01 读剖面名 / 逐位置流量 / 边界
    f_profiles, flows, bnds = read_flow_file(F01)
    # Fall River 各断面流量: 用 (Fall River, reach, rs) 位置的流量
    # (该位置的流量适用于该位置至下一流量位置之间的断面)
    q_by_rs = {}
    for (riv, rch, rs), qlist in flows.items():
        if riv == "Fall River":
            q_by_rs[float(rs)] = qlist

    report = []
    report.append("# L2 验证: raslite vs HEC-RAS (EX1 Fall River, 单河段)\n")
    report.append(f"- HDF: {os.path.basename(HDF)} ({meta.get('File Version', '')})")
    report.append("- 模型: Fall River 7 断面 (RS 10 → 9.5), 支流汇合口按普通断面衔接 (近似)")
    report.append("- 流量: 来自 EX1.f01 (Upper 500, 汇合口下 600, RS9.6 下 650 cfs 等, 已换算 SI)")
    report.append("- 边界: 下游 RS 9.5 正常水深 S0=0.0004 (EX1 计划设置); 单位制 SI\n")

    grand_max = 0.0
    assert len(f_profiles) == len(profiles), "剖面数不一致"
    for ip, pname in enumerate(f_profiles):
        hname = profiles[ip]          # HDF 内部名 (P1/P2/...)
        qs = []
        for s in fall_secs:
            rs = float(s.rs)
            # 河流站距向下游递减: 断面的流量 = 数值上 >= rs 的最小流量位置
            cands = [r for r in q_by_rs if r >= rs - 1e-9]
            rs_key = min(cands) if cands else max(q_by_rs)
            qs.append(q_by_rs[rs_key][ip])
        sections = [cs for cs, _ in xs_list]
        prof = SteadyProfile(sections, qs, regime="subcritical",
                             boundary=("normal", 0.0004))
        results = prof.run()
        by_name = {r.section.name: r for r in results}

        report.append(f"## 剖面 {pname}\n")
        report.append("| RS | Q(cms) | HEC-RAS WS (m) | raslite WS (m) | Δ (cm) | "
                      "HEC K_tot | ras K_tot | ΔK% | α(HEC) | α(ras) |")
        report.append("|---|---|---|---|---|---|---|---|---|---|")
        pmax = 0.0
        for (cs, s) in xs_list:  # US -> DS 显示
            r = by_name[cs.name]
            ref = s.results[hname]
            dws = (r.ws - ref["ws"]) * 100.0
            pmax = max(pmax, abs(dws))
            dk = (r.conveyance - ref["K_total"]) / max(ref["K_total"], 1e-9) * 100
            report.append(
                f"| {s.rs} | {r.q:.2f} | {ref['ws']:.4f} | {r.ws:.4f} | {dws:+.2f} | "
                f"{ref['K_total']:.1f} | {r.conveyance:.1f} | {dk:+.2f} | "
                f"{ref['alpha']:.3f} | {r.alpha:.3f} |")
        report.append("")
        report.append(f"**{pname} 最大 |ΔWS| = {pmax:.2f} cm**\n")
        grand_max = max(grand_max, pmax)

    report.append(f"\n## 总结\n")
    report.append(f"- 全部剖面最大 |ΔWS| = **{grand_max:.2f} cm**")
    report.append(f"- 验收目标: ≤ 1.0 cm")

    text = "\n".join(report)
    print(text)
    out = os.path.join(HERE, "..", "L2_验证报告_EX1.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\n已保存: {out}")
    return grand_max


if __name__ == "__main__":
    main()

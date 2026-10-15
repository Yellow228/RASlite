"""临界水深对标: 收集 test01~04 中 HEC 默认临界的断面, 逐一对比。

输出两类诊断:
  A. 几何一致性: 在 HEC 临界水位处, raslite 的 K/alpha/A 与 HEC 报告值
     是否一致 (不一致 = 几何解析问题)
  B. 判据一致性: K 一致但临界水位不同 = 临界判据/求解差异
"""

import sys
import os
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from raslite import CrossSection, critical_ws
from raslite.hecras_hdf import read_plan_hdf

def collect():
    """收集所有 (算例, 断面, Q, HEC临界水位) 样本"""
    samples = []
    roots = {
        "test01": os.path.join(HERE, "..", "examples", "test01"),
        "test02": os.path.join(HERE, "..", "examples", "test02"),
        "test03": os.path.join(HERE, "..", "examples", "test03"),
        "test04": os.path.join(HERE, "..", "examples", "test04"),
    }
    plans = {
        "test01": ["QYG.p01.hdf"],
        "test02": ["xccez.p01.hdf"],
        "test03": ["QYG.p02.hdf", "QYG.p03.hdf", "QYG.p04.hdf",
                   "QYG.p05.hdf", "QYG.p06.hdf"],
        "test04": ["YXG.p01.hdf"],
    }
    for case, root in roots.items():
        for phdf in plans[case]:
            path = os.path.join(root, phdf)
            if not os.path.exists(path):
                continue
            sections, profiles, meta = read_plan_hdf(path)
            for s in sections:
                for p in profiles:
                    r = s.results[p]
                    if abs(r["ws"] - r["crit_ws"]) < 0.005 and r["q"] > 0:
                        samples.append({
                            "case": case, "plan": phdf,
                            "name": f"{s.river}/{s.reach}@{s.rs}",
                            "sec": s, "q": r["q"], "crit_ws": r["ws"],
                            "hec_K": r["K_total"], "hec_alpha": r["alpha"],
                            "hec_area": r["area_total"],
                        })
    return samples


def main():
    os.chdir(HERE)
    samples = collect()
    print(f"HEC 默认临界样本总数: {len(samples)}\n")

    n_match_K = n_mismatch_K = 0
    crit_diffs = []
    worst = []
    for sm in samples:
        s = sm["sec"]
        try:
            cs = CrossSection(sm["name"], s.points,
                              s.bank_left, s.bank_right, s.n_main,
                              s.n_left, s.n_right, s.reach_lengths,
                              s.c_contract, s.c_expand)
        except Exception as e:
            worst.append((sm, f"构建失败: {e}"))
            continue

        # A. 几何一致性: 在 HEC 临界水位处对比 K / alpha
        try:
            p = cs.props(sm["crit_ws"])
            k_our, a_our = p.conveyance, p.alpha
        except Exception as e:
            worst.append((sm, f"props失败: {e}"))
            continue
        k_rel = (k_our - sm["hec_K"]) / max(sm["hec_K"], 1e-9) * 100
        if abs(k_rel) > 0.5:
            n_mismatch_K += 1
            worst.append((sm, f"K差{k_rel:+.2f}% (our={k_our:.1f} "
                              f"hec={sm['hec_K']:.1f}) "
                              f"alpha our={a_our:.3f} hec={sm['hec_alpha']:.3f}"))
        else:
            n_match_K += 1

        # B. 判据一致性: 我们的临界水位
        try:
            ws_our, _ = critical_ws(cs, sm["q"])
            d = (ws_our - sm["crit_ws"]) * 100
            crit_diffs.append(d)
            if abs(d) > 5:
                worst.append((sm, f"临界差{d:+.1f}cm"))
        except Exception as e:
            worst.append((sm, f"critical失败: {e}"))

    import statistics
    print(f"几何一致 (|ΔK%|<=0.5): {n_match_K}")
    print(f"几何不一致 (>0.5%):   {n_mismatch_K}")
    if crit_diffs:
        crit_diffs_abs = [abs(d) for d in crit_diffs]
        print(f"\n临界水位差异 (cm): n={len(crit_diffs)}")
        print(f"  中位数: {statistics.median(crit_diffs_abs):.2f}")
        print(f"  平均:   {sum(crit_diffs_abs)/len(crit_diffs_abs):.2f}")
        print(f"  最大:   {max(crit_diffs_abs):.2f}")
        print(f"  <=1cm:  {len([d for d in crit_diffs_abs if d <= 1])}")
        print(f"  <=5cm:  {len([d for d in crit_diffs_abs if d <= 5])}")

    print(f"\n--- 问题样本明细 (前 40) ---")
    for sm, msg in worst[:40]:
        print(f"[{sm['case']} {sm['plan']} {sm['name']} q={sm['q']:.2f}] {msg}")


if __name__ == "__main__":
    main()

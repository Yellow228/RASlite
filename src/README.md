# raslite — 一维恒定流水面线推算程序

对标 HEC-RAS 1D Steady Flow 的自研水面线计算核心（第一期：恒定渐变流，无建筑物）。
纯 Python 标准库实现，无第三方依赖。单位制：公制 SI (m, m³/s)。

## 运行验证

```bash
C:/Users/Mechrevo/.workbuddy/binaries/python/versions/3.13.12/python.exe tests/run_all.py
```

当前 25 项测试全部通过（L0 闭式解 + L1 解析解验证）。

## 最小用例

```python
import sys; sys.path.insert(0, "src")
from raslite import CrossSection, SteadyProfile

# 1. 定义横断面（站距-高程点列、滩肩、糙率、到下一断面的距离）
cs = CrossSection(
    "XS1",                       # 断面名
    [(0, 8), (12, 5), (15, 2), (35, 2), (38, 5), (50, 8)],  # (站距, 高程)
    bank_left=15, bank_right=35, # 滩肩 (主槽范围)
    n_main=0.03,                 # 主槽糙率
    n_left=0.07, n_right=0.07,   # 滩地糙率 (也可传 [(station, n), ...] 突变点)
    reach_lengths=(100, 80, 100),# (左滩, 主槽, 右滩) 到下一断面距离
    c_contract=0.1, c_expand=0.3 # 收缩/扩散系数 (HEC-RAS 约定存上游断面)
)
cs.chainage = 1000.0             # 沿流向里程 (下游大) —— 必须手动赋值

# 2. 推算水面线（亚临界流, 从下游边界向上游推）
prof = SteadyProfile(
    sections=[cs1, cs2, ...],    # 顺序任意, 内部按 chainage 排序
    q=250.0,                     # 恒定流量
    regime="subcritical",        # 或 "supercritical"
    boundary=("ws", 3.5),        # 边界: ("ws",值) / ("critical",) / ("normal",S0) / ("rating",fn)
)
results = prof.run()             # 按推算顺序返回 SectionResult 列表

for r in results:
    print(r.section.name, f"WSE={r.ws:.3f}", f"yc={r.critical_depth:.3f}",
          f"V={r.velocity:.2f}", f"Fr={r.froude:.3f}")
    for w in r.warnings:
        print("  [警告]", w)
```

## 模块结构

| 文件 | 内容 |
|------|------|
| `src/raslite/geometry.py` | 横断面几何：A/P/T、分区输水能力 K（HEC-RAS 默认划分法）、α（式 2-11）、形心压强头 |
| `src/raslite/hydraulics.py` | 临界水深（Fr²=1 二分）、正常水深（曼宁二分） |
| `src/raslite/profile.py` | 标准步进法：能量方程迭代（投影→0.7 误差→割线法）、4 种摩阻坡度平均、流量加权长度、收缩/扩散损失、不收敛默认临界水深（复刻 HEC-RAS 行为） |
| `tests/run_all.py` | L0/L1 验证套件 |

## 与 HEC-RAS 的约定对齐

- reach_lengths 与收缩/扩散系数存于**上游断面**，用于它与下游断面之间的河段
- 收缩/扩散判定：下游流速水头 > 上游 → 收缩
- 滩地按 n 值突变点划分输水分区，主槽整体一个分区；分区竖直虚拟分界不计入湿周
- 断面两端竖直端墙计入湿周；水面切边（水深 0）不计
- 迭代容差默认 0.003 m（= HEC-RAS 0.01 ft），最大迭代 40 次
- 不收敛时默认取临界水深并报警（HEC-RAS 同款行为）

## 已验证精度（vs 渐变流 ODE 的 RK4 参考解）

- M1 壅水曲线（亚临界，100 m 间距）：全线最大误差 **0.25 mm**
- 陡坡均匀流（超临界）：全线**精确**（0.0000 mm）
- S2 落水曲线（超临界，12.5 m 间距）：2.46 mm，已验证随断面间距二阶收敛
  （该残差为标准步进法固有的摩阻坡度离散误差，HEC-RAS 同样具有）

## L2 验证：vs HEC-RAS 官方算例 EX1（已通过）

用 `raslite/hecras_hdf.py` 读取 HEC-RAS 6.2 计算结果（EX1.p01.hdf，需用户先用
HEC-RAS 打开算例 Compute 一次），构建 Fall River 7 断面模型（含汇流口、沿程
变流量 500→600→650 cfs、下游正常水深边界），三个重现期剖面对比：

| 剖面 | 最大 \|ΔWS\| |
|------|------|
| 10 yr | 0.26 cm |
| 50 yr | 0.30 cm |
| 100 yr | 0.30 cm |

**全部剖面最大误差 0.30 cm，验收目标 ≤1 cm 达成**（报告：`L2_验证报告_EX1.md`，
脚本：`scripts/ex1_l2.py`，需 `.venv` 里的 h5py：`.venv/Scripts/python.exe scripts/ex1_l2.py`）。

汇合口处理：HEC-RAS junction 能量法的汇合口摩擦长度取自 Junction Information
（EX1 为 60ft），脚本中按此设置；通用 junction 模块在二期实现。

## v1 已知局限（二期计划）

- 无桥/涵/堰等建筑物（二期）
- 无混合流比力判别（二期）
- 无无效流动区（ineffective flow areas）（二期）
- 无 HEC-RAS 几何文件（.g01）导入器（二期，格式已勘察）
- 复式断面多重临界水深取最浅解，与 HEC-RAS 比能最小法可能略有差异

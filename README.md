# raslite — 一维恒定流水面线推算程序

对标 HEC-RAS 一维恒定流模块的 Python 实现：能量方程标准步进法推算水面线，
**零第三方依赖**（核心计算纯标准库），并已在与 HEC-RAS 6.2 逐断面对比中
达到**厘米级一致**。

## 功能

- **标准步进法**：亚临界 / 超临界推算，HEC-RAS 试值策略（投影 → 0.7 误差
  修正 → 割线法）、4 种摩阻坡度平均（默认平均输水能力法）、流量加权断面间距、
  收缩/扩散损失、不收敛默认临界水深并告警
- **混合流**：亚临界遍 + 超临界遍 + 比力 SF 合并，自动定位水跃（复刻
  HEC-RAS Mixed Flow 算法）
- **河网汇流**：`Reach` / `Junction` / `SteadyNetwork`，汇合口按
  reach-by-reach 能量平衡处理
- **变流量**：逐断面流量序列（沿程汇入/分出）
- **临界水深**：HEC-RAS secant 法（30 片切片、找全部局部极小、取比能
  最小者），正确处理复式断面 α 随水位变化与多重临界问题
- **堤防（Levee）**：水位低于堤顶时堤后面积不参与过流，漫堤后恢复全断面
- **HEC-RAS 文件生态**：`.g01` 几何导入（断面/糙率分区/滩肩/汇合口长度）、
  `.f01` 流量与边界解析、`.p01.hdf` 结果读取（对比验证用，需 h5py）

## 验证结果

### 与 HEC-RAS 6.2 逐断面对比（L2）

| 算例 | 内容 | 最大 \|ΔWS\| |
|------|------|------|
| EX1（官方 Fall River，单河段） | 汇流 + 变流量，3 个重现期 | **0.30 cm** |
| MIXED（官方混合流） | 陡坡接缓坡 + 水跃，2 剖面 | **0.22 cm**（水跃位置完全一致） |
| EX1 全河网 | 3 河段 + Junction，10 断面 × 3 剖面 | **0.42 cm** |
| test01 / test02（真实河道） | 临界控制山区河道 | **0.53 / 0.39 cm** |
| test03（真实设计项目） | 5 方案 × 6 频率，130 工况 | 无建筑物方案 **1.00~1.15 cm** |
| test05（真实河道） | 大流量 333 m³/s + 堤防 | **1.30 cm** |

### 与解析解对比（L1，对渐变流 ODE 的 RK4 参考解）

- M1 壅水曲线（100m 断面间距）：**0.25 mm**
- 陡坡均匀流：**精确（0.0000 mm）**
- L0 闭式解（矩形/三角/复式断面几何、临界/正常水深）：全部精确

## 快速上手

```python
import sys; sys.path.insert(0, "src")
from raslite import CrossSection, SteadyProfile

# 断面（站距-高程点列、滩肩、糙率、到下一断面的距离）
cs1 = CrossSection(
    "XS1",
    [(0, 105.0), (50, 98.0), (60, 95.0), (200, 95.0), (210, 98.0), (260, 105.0)],
    bank_left=60, bank_right=210,
    n_main=0.030, n_left=0.055, n_right=0.055,
    reach_lengths=(0, 100, 0),
)
cs1.chainage = 0        # 里程沿流向递增（下游大），必须赋值
# ... cs2, cs3 ...

# 推算水面线
prof = SteadyProfile(
    [cs1, cs2, cs3],
    q=350.0,                      # 流量 m³/s（也可传逐断面列表）
    regime="subcritical",         # / "supercritical" / "mixed"
    boundary=("normal", 0.001),   # 或 ("ws", 高程) / ("critical",)
)
for r in prof.run():
    print(f"{r.section.name}  WS={r.ws:.3f}  V={r.velocity:.2f}  Fr={r.froude:.3f}")
```

直接吃 HEC-RAS 工程：

```python
from raslite.g01 import parse_g01, build_network, read_prj_units
from raslite.hecras_hdf import read_flow_file

unit = read_prj_units("工程.prj")
reaches, junctions = parse_g01("工程.g01", unit=unit)
net = build_network(reaches, junctions,
                    outlet=("River 1", "Reach 1"), boundary=("ws", 652.99))
net.reaches["River 1/Reach 1"].q = 15.4
for r in net.run_subcritical()["River 1/Reach 1"]:
    print(r.section.name, f"{r.ws:.4f}")
```

更多用例见 [`src/README.md`](src/README.md)、
[`scripts/`](scripts/) 下的验证脚本（test01~test05、EX1 官方算例），
验证报告由各对比脚本运行时本地生成（`L2_验证报告_*.md`，不入库）。

## 项目结构

```
src/raslite/
├── geometry.py      # 横断面几何: A/P/T、分区 K（HEC-RAS 划分法）、α、堤防截断
├── hydraulics.py    # 临界水深（secant 法）、正常水深、比力
├── profile.py       # 标准步进法主循环（单河段 + 混合流）
├── network.py       # 河网: Reach / Junction / SteadyNetwork
├── g01.py           # HEC-RAS .g01 几何文件导入器
└── hecras_hdf.py    # .p01.hdf 结果读取 / .f01 流量文件解析
tests/run_all.py     # 验证套件（L0 闭式解 + L1 解析解, 25 项）
scripts/             # 与 HEC-RAS 的对比验证脚本
examples/            # 算例（不入库, 见下方说明）
```

## 运行验证

```bash
# 核心验证（纯标准库, 任意 Python ≥3.10）
python tests/run_all.py

# 与 HEC-RAS 对比（需算例数据 + h5py）: .venv/Scripts/python.exe scripts/test01_run.py
```

## 与 HEC-RAS 对齐的关键约定

- 复式断面输水能力按滩肩划分、滩地按糙率突变点再分（HEC-RAS 默认法）
- 断面间距按左滩/主槽/右滩流量加权
- 收缩/扩散系数存于上游断面，用于本断面与下游断面之间
- 临界水深 = 总水头最小（secant 法，30 片切片取比能最小极小）
- 迭代不收敛默认取临界水深并告警；边界水位低于/高于临界同样默认临界
- 河流站距向下游递减；汇合口能量步只计摩擦

## 已知局限

- 桥、涵、堰等建筑物节点未实现（g01 中 type≠1 自动跳过并告警），
  其上游壅水缺失——二期首要项
- 无效流动区未实现
- 河网仅支持亚临界过汇合口（超临界需动量法）
- `.hdf` 读取需 `h5py`（仅对比验证用，核心计算无依赖）

## 参考资料

- HEC-RAS Reference Manual / 1D Technical Reference（HEC 官网免费公开）
- 相关算例下载：HEC-RAS 7.0 Example Projects（HEC 官网）

## 许可证

[MIT](LICENSE) — 可自由使用、修改与分发，需保留版权声明；软件按"现状"提供，不含任何担保。

## 免责声明

本程序用于工程学习与研究。虽已通过与 HEC-RAS 6.2 的逐断面对比验证，
正式工程设计仍应以经审定的商业软件成果为准。

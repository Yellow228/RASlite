# HEC-RAS 一维水面线推算计算原理（资料汇编）

> 整理日期：2026-10-15
> 主要来源：HEC 官方《Hydraulic Reference Manual》（v5.0 / v6.0）、HEC Confluence 技术文档、国内河道管理范围划定报告中的原理章节等（文末附全部链接）。

---

## 1. 总体定位

HEC-RAS（River Analysis System）由美国陆军工程兵团水文工程中心（HEC）开发。其一维恒定流模块用于计算**天然或人工河道中恒定渐变流（steadily gradually varied flow）的水面线**，可处理**亚临界、超临界及混合流态**，并能考虑桥梁、涵洞、堰、堤防等涉水建筑物的影响。

核心思路：从已知水位的断面出发，**逐断面求解能量方程**（标准步进法 Standard Step Method），迭代求出相邻断面的未知水面高程。

---

## 2. 基本方程：能量方程

相邻两断面（1 为下游、2 为上游）间的能量守恒：

$$Z_2 + Y_2 + \frac{\alpha_2 V_2^2}{2g} = Z_1 + Y_1 + \frac{\alpha_1 V_1^2}{2g} + h_e$$

| 符号 | 含义 |
|------|------|
| Z₁, Z₂ | 下/上游断面主槽河底高程 |
| Y₁, Y₂ | 下/上游断面水深 |
| V₁, V₂ | 断面平均流速（总流量 / 总过水面积） |
| α₁, α₂ | 动能修正系数（流速权重系数） |
| g | 重力加速度 |
| hₑ | 两断面间能量水头损失 |

> 说明：该方程是伯努利方程在明渠中的简化——忽略大气压变化，保留位置水头、压强水头（用水深表示）与流速水头。

---

## 3. 能量损失 hₑ 的构成

$$h_e = L\,\bar{S}_f + C\left|\frac{\alpha_2 V_2^2}{2g} - \frac{\alpha_1 V_1^2}{2g}\right|$$

即 **沿程摩阻损失（L·S̄f）＋ 收缩/扩散局部损失** 两部分。

### 3.1 流量加权断面间距 L

复式断面按左滩、主槽、右滩三部分分别量取距离，再按流量加权：

$$L = \frac{L_{lob}\bar{Q}_{lob} + L_{ch}\bar{Q}_{ch} + L_{rob}\bar{Q}_{rob}}{\bar{Q}_{lob} + \bar{Q}_{ch} + \bar{Q}_{rob}}$$

其中 L_lob / L_ch / L_rob 为左滩/主槽/右滩的断面间距，Q̄ 为相应分区在两断面间流量的算术平均值。

### 3.2 代表性摩阻坡度 S̄f（4 种可选，默认第 1 种）

各断面摩阻坡度由曼宁公式反算：S_f = (Q/K)²，K 为断面输水能力（conveyance）。

1. **平均输水能力法（Average Conveyance，HEC-RAS 默认）**：
   $$\bar{S}_f = \left(\frac{Q_1+Q_2}{K_1+K_2}\right)^2$$
2. **算术平均摩阻坡度**：S̄f = (S_f1 + S_f2)/2
3. **几何平均**：S̄f = √(S_f1 · S_f2)
4. **调和平均**：S̄f = 2·S_f1·S_f2 / (S_f1 + S_f2)

用户可根据流态与水面线类型（M1、S1 等）在程序中切换。

### 3.3 收缩 / 扩散损失系数 C

- 程序判定规则：下游流速水头 > 上游 → 发生**收缩**；反之 → 发生**扩散**。
- C 取值范围：收缩 0.0～0.5（默认约 0.1），扩散 0.0～1.0（典型 0.3～0.5），突变过渡取大值。

---

## 4. 断面划分与输水能力计算

### 4.1 输水能力 K

HEC-RAS 将滩地按**糙率 n 值突变点**划分为若干子区（主槽通常整体作为一个分区），每一分区用曼宁公式：

$$Q = K\sqrt{S_f}, \qquad K = \frac{1.486}{n} A R^{2/3} \ \text{(英制；公制系数为 1/n)}$$

总输水能力 = 左滩 + 主槽 + 右滩各分区 K 值之和。另有 HEC-2 风格的"逐坐标点划分"备选方法（为复现老成果保留）。

### 4.2 复合糙率

主槽由不同材料组成时，按各子区面积加权计算复合 Manning n 值。

### 4.3 动能修正系数 α

按三区输水能力与面积计算：

$$\alpha = \frac{(A_t)^2\left[\dfrac{K_{lob}^3}{A_{lob}^2} + \dfrac{K_{ch}^3}{A_{ch}^2} + \dfrac{K_{rob}^3}{A_{rob}^2}\right]}{K_t^3}$$

---

## 5. 标准步进法迭代计算流程（亚临界流为例）

1. **假定**上游断面一个水面高程（超临界流则自上游向下游推算，假定下游水位）；
2. 由假定水位计算该断面的总输水能力 K 与流速水头；
3. 由第 2 步数值计算 S̄f，代入损失公式求 hₑ；
4. 代入能量方程，解出上游水面高程 WS₂；
5. 比较计算值与假定值：不一致则重新假定，直至差值 ≤ 0.01 ft（约 0.003 m）或用户自定义容差。

**试值修正策略**：
- 第 1 次试值：将上一断面水深直接投影到当前断面；
- 第 2 次试值：假定值 + 0.70 ×（计算值 − 假定值）；
- 第 3 次及以后：**割线法（Secant method）**外推，利用前两次试值的误差变化率加速收敛。

### 边界条件（4 种）

亚临界流须在**最下游断面**给定，超临界流在最上游断面给定：
1. 已知水面高程；
2. 临界水深；
3. 正常水深（需输入河床比降，用曼宁公式求）；
4. 水位—流量关系曲线（率定曲线）。

计算河段应足够长：下游起始断面离研究区域足够远、上游延伸至水面线收敛于正常水深，以免边界条件误差影响成果。

---

## 6. 临界水深与混合流计算

- **临界水深**：断面比能最小时的水深（Froude 数 Fr = 1）。能量方程在临界流附近不适用，程序在无法收敛时默认取临界水深并做标记。
- **混合流（Mixed Flow Regime）**：用**动量方程导出的比力（Specific Force）**判别控制流态：

$$SF = \frac{Q^2\beta}{gA} + A\bar{Y}$$

计算分两步：
1. 先自下游边界向上推**亚临界**水面线，把所有"默认取临界水深"的位置做标记；
2. 再自上游边界向下推**超临界**水面线；在同时存在亚临界与超临界两解的断面，比较两者的比力——**比力大者控制**；当某断面的亚临界解比力更大时，判定其在与上一断面之间发生了**水跃**。

---

## 7. 动量方程的适用场合

水流**急变**（断面突变、流态转换）时能量方程失效，改用动量方程：

$$P_2 - P_1 + W_x - F_f = \rho Q \Delta v_x$$

（P 为静水压力，W_x 为水体重量沿流向分量，F_f 为摩阻力。）典型应用：
- 混合流态与水跃定位；
- 桥梁壅水计算（低流 Class B 必用动量法；Class A 可选能量法/动量法/Yarnell 公式/USGS 收缩开口法）；
- 河道汇流（干支流交汇）处的能量或动量平衡；
- 消能设施、闸孔出流等。

---

## 8. 模型主要假设与局限（官方提示）

- 一维恒定、渐变流假设：断面上水压力按静水压力分布，流线近似平行；
- 断面间距应使水流在每一子河段内"渐变"，断面突变处（桥、涵、陡坡）需专门处理；
- 弯道横向水面超高、急变流的非静水压力分布未在基本方程中考虑（有单独的经验修正）；
- 高速水流需考虑掺气（air entrainment）对水深与摩阻的影响；
- 恒定流模块不能模拟洪水波的传播过程（需用一维非恒定流模块，解圣维南方程组）。

---

## 9. 参考资料（原始链接）

**官方权威文档**
1. HEC-RAS 5.0 Hydraulic Reference Manual（第 2 章"Basis of 1D Steady Flow Water Surface Profiles"为核心）：
   https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS%205.0%20Reference%20Manual.pdf
2. HEC-RAS 6.0 Hydraulic Reference Manual（2021 正式版）：
   https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS_6.0_Reference_Manual.pdf
3. HEC Confluence 在线技术参考——Friction Loss Evaluation（4 种摩阻坡度平均法）：
   https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.0/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-steady-flow-water-surface-profiles/friction-loss-evaluation
4. HEC Confluence——Mixed Flow Regime Calculations（比力方程与混合流流程）：
   https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.5/overview-of-optional-capabilities/mixed-flow-regime-calculations
5. HEC 官方培训讲义 "Water Surface Profile Calculations with HEC-RAS"（幻灯片，适合快速入门）：
   https://www.hec.usace.army.mil/confluence/prospect/files/124452884/124453107/1/1669663021202/1.6+L+-+Water+Surface+Profile+Calculations.pdf
6. HEC 技术论文 TP-147（含混合流、桥梁/涵洞例程、河网汇流处理）：
   https://www.hec.usace.army.mil/publications/TechnicalPapers/TP-147.pdf

**方法学补充**
7. Standard Step Method 条目（含 E-y 曲线、水面线型 M/S 分类、Froude 数）：
   https://www.wikiwand.com/en/articles/Standard_step_method
8. 印度水资源 portals 的 HEC-RAS 水面线计算综述（公式推导完整、含算例）：
   https://prod-qt-images.s3.amazonaws.com/indiawaterportal/import/sites/default/files/iwp2/Computation_of_Water_Surface_Profiles_using_HEC_River_Analysis_System.pdf

**中文应用文献（原理章节均转述官方手册，可作中文表述参考）**
9. 三亚市抱前一溪河道管理范围线调整报告（§3.2 水面线计算原理，中文公式对照）：
   https://ty.sanya.gov.cn/tyqsite/qzfwjxx/202511/8db3ffcd424b4a5183416cd4e79d99cb/files/4da80896d1fa4c56b0ba57aab3a71d25.pdf
10. 益阳市资阳区水口山河管理范围划定方案（§4 中文能量方程表述）：
    http://www.hnziyang.gov.cn/uploadfiles/202312/2023121913281560171.pdf

**经典教材**（方法学源头）
- Chow, V.T. (1959). *Open Channel Hydraulics*——标准步进法、比能曲线、水面线分类（M1/M2/S1…）的原始出处；HEC-RAS 混合流示例即取自该书习题 9-8。
- HEC, *Hydrologic Engineering Methods for Water Resources Development, Vol.6: Water Surface Profiles*。

---

## 10. 快速记忆卡片

| 问题 | 答案 |
|------|------|
| 一维水面线的基本方程？ | 相邻断面能量方程（Z+Y+αV²/2g 守恒 + hₑ） |
| 求解方法？ | 标准步进法，逐断面迭代（割线法加速），容差默认 0.003 m |
| 损失项？ | L·S̄f（摩阻）+ C·Δ流速水头（收缩/扩散） |
| S̄f 默认算法？ | 平均输水能力法 (Q₁+Q₂)²/(K₁+K₂)² |
| 断面间距 L？ | 左滩/主槽/右滩按流量加权 |
| 亚临界从哪端算起？ | 已知下游边界，向上游推；超临界相反 |
| 流态转换怎么判？ | 动量方程 → 比力 SF 大者控制，判定水跃位置 |
| 桥梁壅水？ | 低流 Class A 可四法选一；Class B/压力流/急变处用动量方程 |

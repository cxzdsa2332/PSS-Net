# PSS-Net 创新点评价与提升方向（2026-06 文献调研后）

## 1. 当前方法的创新定位（诚实评价）

方法 = 扰动稳态(PSS) + 非参加性 ODE 基展开 + 双稀疏 ADSIHT + 微生物应用。
**每个成分均来自已有工作**，组合属增量贡献：

| 成分 | 来源 | 原创性 |
|------|------|--------|
| 扰动稳态（免数值微分） | Meister 2013、Xiao 2017 | 已有 |
| 非参加性 ODE + 基展开 | Henderson & Michailidis 2014、Wu 2014 SA-ODE | 已有 |
| 双稀疏 ADSIHT | Zhang 2023 | 调包 |
| 微生物互作应用 | Xiao 2017、Venturelli/Clark | 已有方向 |

**核心风险**：易被质疑为"Meister 2013 + Henderson 2014 + ADSIHT 的拼装"。

**几个卖点并不独占**：
- 免数值微分 → Xiao 2017、Meister 2013 已有；
- 稳态推断 → Xiao 2017 核心；
- 与最接近的 idopNetwork（Dong 2023）相比，仅多"主动扰动 + 双稀疏"，差异不够。

**唯一真正差异内核**：相对 Xiao/idopNetwork 的 gLV 线性，可估计非线性互作函数 $f_{ji}(\cdot)$；
相对 Meister/Henderson，有主动扰动设计 + 组内稀疏。但两点都未做深。

## 2. 可提升的创新组合（按推荐度）

### 🥇 A. 最优/主动扰动设计（把"推断"升级为"设计+推断"）
当前 $u^{(k)}\sim\text{Uniform}(-0.4,0.8)$ 为随机扰动。扰动是 PSS 独占杠杆——
**无人形式化"为可识别性该施加哪些扰动"**。
- 在逐节点 $\Psi_c$ 上做 D-/A-最优设计，或序贯自适应扰动选择（按 Fisher 信息选下一条件）。
- 价值：Xiao/idopNetwork/SINDy 都没有的维度；对湿实验直接省成本。**最强单点提升。**

### 🥈 B. 稳态可识别性 + 样本复杂度理论
回答"N 个扰动稳态何时能恢复网络"：给 $U$ 的秩/覆盖条件，把 Zhang 2023 minimax 率
经稳态映射传到 $N\gtrsim s\log p$ 边恢复保证。与 A 配套。

### 🥉 C. 超越两两的高阶互作（蹭 Maynard 2024 PNAS 热点，代码已有雏形）
ref/v0.1.txt 已模拟 $f_{ij}=x_ix_j$ 二阶项。扩为 functional ANOVA（主效应 + 选中的成对
交互函数）+ 层次组稀疏。结构上甩开所有 gLV 方法（Xiao/idopNetwork/MDSINE）。

### 可信度层（真实数据必备，非亮点）
- D. 组成型数据：16S 需 CLR/log-ratio 嵌入稳态回归（datasets.md 已注明）。
- E. Errors-in-variables：$x$ 同时在含噪观测与 $\psi(x)$ 中，经典 EIV 偏差；
  做测量误差校正的基回归，也解释"预平滑反而有害"的现象。

## 3. 推荐新定位

> **"主动扰动稳态下，带可识别性保证的稀疏非参（含高阶）微生物互作网络重构"**
> 骨架 A（最优扰动设计）+ B（可识别性理论）；差异化 C（高阶交互）；可信度 D/E。
> 精力有限时**优先 A**（投入产出比最高，wet-lab 买账）。

## 参考文献（详见 references.bib）
meister2013pss, henderson2014network, wu2014saode, zhang2023minimax,
xiao2017mapping, dong2023idopnetwork, chen2019omnidirectional,
brunton2016sindy, ishizawa2024beyond, kurtz2015sparse.

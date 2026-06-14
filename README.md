# PSS-Net

**稀疏加性非参 ODE 模型用于微生物互作网络重构**
(Sparse Additive Nonparametric ODE for Microbial Interaction Network Reconstruction)

## 项目简介

PSS-Net 旨在从**扰动稳态（Perturbed Steady-State, PSS）**实验数据中重构微生物群落的**稀疏非线性互作网络**。

核心思想：对 $p$ 个物种组成的群落，其动力学由耦合 ODE 系统描述

$$\frac{dx_j}{dt} = \mu_j + \sum_{i=1}^{p} f_{ji}(x_i) + u_j, \quad j = 1,\ldots,p$$

其中 $f_{ji}(\cdot)$ 为物种 $i$ 对 $j$ 的未知单变量互作函数。在多个扰动条件下系统达到稳态（$\dot{x}_j = 0$）后，方程退化为代数约束，无需数值微分即可建立线性回归问题。

技术路线：

1. **稳态方程**：扰动稳态下 $\mu_j + \sum_i f_{ji}(x_i^*) + u_j = 0$，响应变量为 $-u_j$。
2. **非参基展开**：用 B 样条（或单项式）基 $\psi(x)$ 近似每个 $f_{ji}$，约束 $f_{ji}(0)=0$（无截距基）。
3. **逐节点回归**：对每个目标物种 $j$ 独立构建设计矩阵 $\Psi_j$，列中心化消除截距。
4. **双稀疏推断**：用 **ADSIHT** 求解组间稀疏（哪些边存在）+ 组内稀疏（互作函数低阶），DSIC 自动选模型。
5. **Jacobian 提取与网络判定**：在野生型稳态处求 $J_{ji} = \psi'(x_i^{\mathrm{wt}})^\top \hat\theta_{ji}$，由 $\|\hat\theta_{ji}\|_2 \ge \tau$ 判定边，符号区分促进/抑制。

详细方法与公式见 [CLAUDE.md](CLAUDE.md)。

## 目录结构

| 目录 | 作用 |
|------|------|
| `data/` | 真实数据（稳态丰度矩阵、扰动设计等）。`datasets.md` 列出按 PSS 框架适配度排序的候选公开数据集（Clark 2021、Venturelli 2018、MDSINE2、Stein 2013、Maier 2024）。 |
| `sim_script/` | 模拟与推断脚本。`pss_net_*.R` 为 PSS 框架主线（v0 大规模欠定、v1 含非线性互作、v3 线性 GLV、`pss_net_compare.R` 为多次重复 MCC 基准）；`sindy_ss_*.R` 为 SINDy 风格对照实现。 |
| `analysis_script/` | 分析脚本（ODE 推断、网络重构、可视化）。当前待填充。 |
| `methods/` | 方法说明文档（如 `sindy_ss_method.md`）。 |
| `results/` | 输出结果。`mcc_comparison.csv/.txt` 记录 ADSIHT vs Group Lasso 在不同方案下的 TP/FP/FN、Precision/Recall/F1/MCC 等指标。 |
| `ref/` | 参考文献（PDF、`references.bib`、`ODE_solve.md` 等）。 |
| `CLAUDE.md` | 项目完整方法学说明与编码规范（建模范式、公式推导、ADSIHT 调用约定）。 |

## 技术栈

- **主语言**：R（建模、分析、绘图）；Python（辅助计算）
- **ODE 求解**：R `deSolve` / Python `scipy.integrate`
- **基展开**：`splines::bs()`（无截距 B 样条）或 Legendre 多项式
- **稀疏回归**：`ADSIHT`（首选）> `sparsegl`/`gglasso`（次选）> `glmnet`（基线）
- **绘图**：统一 `ggplot2`；网络图配色见 CLAUDE.md（促进=`tomato3`、抑制=`steelblue3`、FP=`orange2`、FN=`grey60`）

## 当前进展

- 模拟框架（线性 GLV、含非线性互作）已跑通；ADSIHT 在无平滑方案下 MCC 显著优于 Group Lasso（详见 `results/`）。
- 已确认：B 样条预平滑在多物种系统中**降低**推断质量，PSS 推断直接使用含噪观测。
- 真实数据集已调研并排序，尚未接入分析流程。

后续待办见 [CLAUDE.md](CLAUDE.md) 末尾的 **TODO List**。

## 参考文献

- Henderson & Michailidis (2014) *Network reconstruction using nonparametric additive ODE models.*
- Wu et al. (2014) *Sparse additive ODEs for dynamic gene regulatory network recovery.*
- Zhang et al. *Minimax optimal estimation via ADSIHT.*
- Barzel & Barabási (2013) *Universality in network dynamics.*

完整 BibTeX 见 [ref/references.bib](ref/references.bib)。

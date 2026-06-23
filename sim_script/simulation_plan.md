# PSS-Net 模拟研究规划

本文档用于把正式模拟脚本组织成论文级模拟路线。当前主文建议围绕三张多面板图展开：图 1 证明 PSS-Net 的可识别性与基础恢复能力，并在视觉上明确真实生成机制是 ODE 动力学、PSS 数据只是稳态切片；图 2 证明样本复杂度与扰动设计价值；图 3 证明鲁棒性、结构依赖与外部 benchmark。模型错设与敏感性分析放入 `sup/`；真实数据作为后续 case study，不混入模拟脚本。

## 总体原则

主文模拟需要回答三个审稿人会直接追问的问题：

1. 扰动稳态（perturbed steady state, PSS）数据能识别什么，不能识别什么？
2. PSS-Net 是否能在高维稀疏网络中随样本预算增加而稳定恢复有向耦合网络？
3. 与现有网络推断方法相比，PSS-Net 使用扰动输入 `u` 能带来什么额外信息？

所有正式数值模拟脚本应只生成数值结果，写入 `results/sim_results/`。当前主图汇总脚本仍分别位于 `sim_script/01_foundation_recovery/Fig1.R`、`sim_script/02_scaling_design/Fig2.R` 和 `sim_script/03_robustness_benchmarks/Fig3.R`，读取数值结果并在工作区生成绘图对象；`analysis_script/` 保留通用汇总与绘图脚本。`sim_script/manual/` 仍保留探索性脚本，可以混合绘图和检查。

## 图 1：PSS 数据来源、可识别边界与基础恢复能力

目标：图 1 作为方法主图，先把 PSS-Net 的数据对象说清楚：真实系统是 ODE 动力学，实验只观测每个扰动条件达到稳态后的 `(x*, u)`；随后展示 PSS 方程能识别稳态函数形状、但不能单独区分所有瞬态机制；最后用同一批可控模拟证明 node-wise ADSIHT 可以恢复函数形状和有向带符号网络，并在不同维度下保持合理的基础恢复能力。

对应文件夹：`sim_script/01_foundation_recovery/`

模拟数据简述：Fig1a-Fig1c 使用历史 8 节点设定和/或其高维扩展，展示 ODE 到 PSS 测量、稳态函数可识别性，以及 ADSIHT 与 group lasso 的基础比较。Fig1d-Fig1f 使用一个固定 10 节点非线性加性 ODE 系统，在 `N = 200` 个扰动稳态条件和 `SNR = 30` 的观测噪声下拟合 node-wise ADSIHT，用于展示动态效应分解、边函数形状恢复，以及真实网络与推断网络的并排比较。

视觉规则：Fig1 当前以 `sim_script/01_foundation_recovery/Fig1.R` 为准。Fig1a 使用浅色背景区分 additive ODE 与 gLV ODE，稳态采样用黑色垂直虚线标注，节点名直接标在线末端。Fig1b 使用稳态函数曲线、BIC gain 与 SNR 扫描说明可识别边界。Fig1d-Fig1e 采用 true vs ADSIHT 的实线/虚线对照。Fig1f 使用 igraph 绘制真实与推断网络，促进作用为红色，抑制作用为蓝色；推断网络中 TP/FP/FN 用线型和灰度辅助标记。

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig1a_ode_to_pss_measurement` | 展示 ODE 动力学如何在扰动下产生 PSS 测量值 `(x*, u)`。当前为 additive ODE 与 gLV ODE 两类动力学，包含 baseline、single-node input 和 mixed input。 | `sim_script/01_foundation_recovery/Fig1.R` 中 `Fig1a`、`Fig1a_legend`。 | 初版完成。用于说明 PSS-Net 的观测对象是扰动稳态切片，而非时间序列。 |
| `Fig1b_steady_state_function_identifiability` | 展示 PSS 不能区分等价线性机制，但能区分真正非线性的稳态函数结构；同时用 SNR 扫描展示非线性识别的噪声边界。 | `Fig1.R` 中 `Fig1b_function_shape`、`Fig1b_identifiability`、`Fig1b_noise_snr`、`Fig1b`、`snr_summary`。相关说明：`note/identifiability_glv_vs_additive.md`。 | 初版完成。当前结果支持：linear additive 与 standard gLV 不应被过度区分；真正 nonlinear additive 在足够 SNR 下可被识别。 |
| `Fig1c_adsiht_vs_group_lasso_scaling` | 在 `p = 8, 30, 100` 下比较 node-wise ADSIHT 与 group lasso，并同时覆盖 linear 与 nonlinear truth，展示 MCC、AUPRC、AUROC、`CoefL2` 和 `JacRMSE` 随 `N / (s log p)` 的变化。 | 模拟：`sim_script/01_foundation_recovery/Fig1c_adsiht_group_lasso_scaling.R`。结果：`results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv`。绘图对象：`Fig1.R` 中 `Fig1c_method_scaling`、`Fig1c`。 | 模拟与绘图初版完成。固定 `SNR = 30`，噪声按 `sigma_x = signal_scale / 30` 添加；当前 `R = 30`。 |
| `Fig1d_effect_decomposition` | 在固定 10 节点非线性加性系统中，展示 ADSIHT 推断的 self effect 与 received regulation 沿扰动后轨迹积分后能重构目标节点状态变化。 | `Fig1.R` 中 `Fig1d`；同一代码块中保留禁用的 `Fig1d_dynamics` 用于后续拓展。 | 初版完成。强调函数估计结果可以回到动态效应解释层面，而不仅是静态边集。 |
| `Fig1e_function_shape_recovery` | 展示同一 10 节点系统中，node-wise ADSIHT 对 self feedback 与 cross-node edge functions 的曲线恢复。 | `Fig1.R` 中 `Fig1e`。 | 初版完成。用于支撑 PSS-Net 能恢复稳态函数形状。 |
| `Fig1f_true_vs_inferred_network` | 使用 igraph 并排绘制同一 10 节点系统的真实有向网络与 ADSIHT 推断网络。促进作用红色，抑制作用蓝色；推断图显示 TP、FP、FN 和 MCC。 | `Fig1.R` 中 `Fig1f`；最终总图对象为 `Fig1`。 | 初版完成。与 Fig1d/Fig1e 共用同一 seed、同一 `SNR = 30` 拟合结果，作为函数恢复到网络恢复的视觉总结。 |

## 图 2：样本复杂度与扰动设计

目标：证明 PSS-Net 在高维稀疏设定中有合理样本复杂度，并展示主动扰动设计能提高信息效率。

**贡献边界：**这一部分按 **A（PSS/PSS-Net 的扰动—稳态映射与网络恢复任务）+
B（已有最优实验设计）**定位。D-optimality 与 Fedorov exchange 不作为本文方法创新；
Fig2b/Fig2e 的选点后端直接使用 CRAN `AlgDesign::optFederov()`，项目代码只负责生成 PSS
候选特征并评估网络恢复。算法出处、包版本和实现映射见
`ref/pilot_doptimal_literature.md` 与 `ref/references.bib`。

对应文件夹：`sim_script/02_scaling_design/`

当前视觉叙事依次为：Fig2a 回答“需要多少样本”；Fig2b 用 3D 输入空间说明 random、maximin、oracle D-optimal 和 pilot-estimated D-optimal 如何从共同 pilot 出发继续选点；Fig2c 汇总不同非线性强度下相对 random 的设计增益；Fig2d 把性能差异换算为达到目标 MCC 所需的实验预算；Fig2e 检验未知真实系统时，pilot-estimated D-optimal 能保留多少 oracle 设计价值。`sim_script/02_scaling_design/Fig2.R` 生成 `Fig2a`--`Fig2e` 及总图对象 `Fig2`。

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig2a_highdim_sample_complexity` | 展示 node-wise ADSIHT 的网络恢复是否随重标度样本量 `N / (s log p)` 改善；当前维度为 `p = 8, 50, 100`，其中 `p = 8` 与 Fig1 小网络对应。标注理论参考线 `N = 1, 2, 5 x s log(p)` 以及达到 `MCC > 0.8` 所需的插值预算。 | 模拟：`sim_script/02_scaling_design/Fig2a_highdim_sample_complexity.R`。结果：`results/sim_results/Fig2a_highdim_sample_complexity.csv`。绘图对象：`sim_script/02_scaling_design/Fig2.R` 中 `Fig2a_highdim_sample_complexity`、`Fig2a`。 | 模拟与绘图初版完成。当前 `R = 5`、加性观测噪声 `sigma = 0.03`，只展示 ADSIHT 的 MCC sample-complexity 曲线；主文前需增加重复数，并补充 AUPRC/AUROC 和失败率。 |
| `Fig2b_design_mechanism` | 在同一个 3D 扰动空间中并排展示 random、maximin、oracle D-optimal 和 pilot-estimated D-optimal。四种策略共享 8 个随机 pilot 条件（黑色外圈）；oracle 使用玩具非线性真实响应映射，pilot 版本仅用含噪 pilot 数据经 ridge 回归估计局部响应映射，再由 `AlgDesign::optFederov()` 选择其余 20 个条件。 | 概念模拟：`sim_script/02_scaling_design/Fig2_explain.R`。绘图对象：`Fig2.R` 中 `Fig2b_design_mechanism`、`Fig2b`；独立对象包括 `Fig2_design_concept_oracle` 与 `Fig2_design_concept_pilot`。 | 已完成并纳入主图。当前示意中 oracle 与 pilot 的 20 个后续选点重合 7 个，能够同时显示共同的边界采样倾向与模型估计误差造成的选点差异；该面板只解释机制，不承担性能比较。 |
| `Fig2c_design_gain_map` | 直接回答“经过设计的 `u` 是否比 random 好”。每格基于同一系统 seed、同一预算的配对比较，显示 maximin 或包实现 D-optimal 相对 random 的 mean `Delta MCC`，第二行显示正增益 seed 百分比；圆点标记该预算下平均增益较大的 structured design。 | 数值来源：`Fig2b_design_linear.csv`、`design_nl_comparison.csv`、`design_nl_seq_comparison.csv`；三份模拟均已统一为 `AlgDesign::optFederov(frml = ~ ., criterion="D")`。绘图对象：`Fig2.R` 中 `fig2c_gain`、`Fig2c_design_gain_map`、`Fig2c`。 | 已完成包后端复算。30 个 regime×预算×structured-method 单元中有 26 个平均增益为正；两种非线性 regime 的 20 个单元全部为正。线性 regime 的优势主要位于小预算（maximin `Delta MCC = +0.03` 至 `+0.05`），预算充足后与 random 基本持平。标准 exact D-optimal 模型含 16 个 PSS 特征和截距，因此只在 `N >= 17` 时定义；扫描网格中线性从 `N=17`、非线性从 `N=20` 展示 D-optimal，不用自写 ridge 准则填补欠定区。 |
| `Fig2d_budget_to_target` | 将均值学习曲线换算为达到 `MCC = 0.5` 和 `MCC = 0.6` 所需的插值预算 `N/(s log p)`；箭头从 random 指向 matched structured design，只有左向箭头才表示节省。超过最大扫描预算仍未达标时显式标记 `> max`。 | 数值来源同 Fig2c；绘图对象：`Fig2.R` 中 `fig2d_budget`、`Fig2d_budget_to_target`、`Fig2d`。 | 已完成复算。达到 `MCC=0.5`：普通非线性中 random 约需 `N=60`，maximin/D-optimal 约需 `N=35/38`；强非线性中 random 约需 `N=57`，maximin/D-optimal 约需 `N=29/35`。达到 `MCC=0.6` 时，两种非线性的 random 在扫描上限 `N=60` 内仍未达到，而 maximin/D-optimal 的插值预算分别约为 `48/52` 和 `45/48`。线性系统只显示小预算的轻微 maximin 收益，D-optimal 因可识别性下限并不节省预算。 |
| `Fig2e_oracle_vs_estimated_design` | 区分 oracle 与 pilot-estimated D-optimal，评估用有限 pilot 估计 PSS 稳态映射造成的设计损失。四种策略共享同一个 random pilot，且 pilot 严格计入总预算；候选 PSS 特征生成后，每个总预算分别调用 `AlgDesign::optFederov(..., criterion="D", augment=TRUE)`，把 pilot 作为不可交换的 protected runs。该实现是一次 pilot 后的 two-stage exact augmentation，不是每批重估模型的在线 adaptive Wynn。 | 模拟：`sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R`。结果：`results/sim_results/Fig2e_oracle_vs_estimated_design.csv`。绘图对象：`Fig2.R` 中 `Fig2e_design_flow`、`Fig2e_regret`、`Fig2e_oracle_vs_estimated_design`、`Fig2e`。 | 已用包后端完成重跑：`p = 8`、`R = 20`、`sigma = 0.04`，`N_pilot = 8, 12, 16`，`N_total = 20, 30, 40, 60`，共 960 行且无缺失。pilot D-optimal 在 12 个组合中的 10 个平均 MCC 高于 random；两个例外均在 `N_total=20`（pilot 为 12 或 16）。因此该图支持“多数设定有收益”，不支持“所有设定稳定占优”。 |

## 图 3：鲁棒性、网络结构与外部 benchmark

目标：把 PSS-Net 放入 bioinformatics/network inference 同行语境中比较，说明其优势来自使用扰动输入和稳态方程，同时诚实展示组成型数据和网络结构依赖。

对应文件夹：`sim_script/03_robustness_benchmarks/`

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig3a_external_benchmark_main` | 在匹配模拟 PSS 数据上，比较 PSS-Net 与外部网络推断方法。最低方法集：PSS-Net/ADSIHT、lasso/elastic net、group lasso、STLSQ/SINDy-style、GENIE3/GRNBoost2、correlation/graphical lasso。 | 参考：`ref/external_benchmark_methods.md`。建议模拟：`sim_script/03_robustness_benchmarks/Fig3a_external_benchmark_main.R`。建议绘图：`analysis_script/Fig3a_external_benchmark_main.R`。 | 缺失，优先级最高。指标应包含 MCC、AUPR/AUROC、precision、recall、sign accuracy、runtime、failure rate。 |
| `Fig3b_method_capability_matrix` | 建立方法能力矩阵：输入类型、输出网络、是否有方向、是否估计非线性边函数、是否需要时序、是否处理 compositional data。 | 参考：`ref/external_benchmark_methods.md`。建议绘图：`analysis_script/Fig3b_method_capability_matrix.R`；可先用手工 CSV/tibble。 | 缺失。适合作为 Fig3 的解释性面板。 |
| `Fig3c_compositional_data_limitation` | 展示 absolute abundance、relative abundance、CLR、relative abundance times noisy total 对恢复的影响；当前绘制 MCC、Precision、Recall 的 mean ± SD，并以 absolute abundance 作为 oracle 上界参考。 | 模拟：`sim_script/03_robustness_benchmarks/Fig3c_compositional_data_limitation.R`。结果：`results/sim_results/Fig3c_compositional_data_limitation.csv`。绘图对象：`sim_script/03_robustness_benchmarks/Fig3.R` 中 `Fig3c_compositional_data_limitation`、`Fig3c`。 | 模拟与绘图初版完成。当前 `p = 10`、`N = 200`、`R = 10`、相对测量噪声 `sigma_rel = 0.05`、总量测量 `CV = 0.15`；noisy total 校正未恢复 oracle 表现，后续应补 total-biomass CV sweep。 |
| `Fig3d_scalefree_joint_vs_nodewise` | 在 scale-free / hub-like 网络中比较 joint block-diagonal estimation 与 node-wise estimation。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`。结果：`results/sim_results/scalefree_compare.csv`。建议绘图：`analysis_script/Fig3d_scalefree_joint_vs_nodewise.R`。 | 数值模拟脚本已有，绘图缺失。当前设定为 `p = 50`、平均入度约 2、`R = 5`、`sigma = 0.03`；主文使用前需增加 R 并明确 NA/失败处理。 |
| `Fig3e_scalefree_hub_recovery` | 评估 joint estimation 是否更稳定恢复核心 hub 排名。指标包括 estimated out-degree vs true out-degree Spearman、top-k hit rate。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`。结果：`results/sim_results/scalefree_compare.csv`。建议绘图：`analysis_script/Fig3e_scalefree_hub_recovery.R`。 | 数值模拟已有，绘图缺失。可与 Fig3d 共用同一数据。 |
| `Fig3f_uniform_joint_favorable_setting` | 在均匀入度、同质信号、极小样本和较高噪声这一 joint-friendly 设定中，检验全局模型选择是否比 node-wise estimation 更稳定；该面板不是 hub 依赖性的 negative control。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_joint_smalln.R`。结果：`results/sim_results/joint_smalln.csv`。建议绘图：`analysis_script/Fig3f_uniform_joint_favorable_setting.R`。 | 数值模拟脚本已有，绘图缺失。当前设定为 `p = 50`、固定入度 3、`R = 8`、`sigma = 0.04`。若要证明优势依赖 hub-like 结构，仍需另建真正的结构负对照。 |
| `Fig3g_targeted_design_tradeoff` | 若先验知道 hub，评估集中扰动 hub 是否提高 hub 相关边恢复、是否牺牲整体网络。由图 2 迁入：定位为网络结构依赖性分析（扰动设计 × 网络结构的交互），与 Fig3d/Fig3e 的 hub 恢复主题相邻。 | 现有结论记录：`note/targeted_perturbation_tradeoff.md`。建议模拟：`sim_script/03_robustness_benchmarks/Fig3g_targeted_design_tradeoff.R`。建议绘图：`analysis_script/Fig3g_targeted_design_tradeoff.R`。 | 当前工作区中正式模拟脚本与 `results/sim_results/targeted_design.csv` 均不存在；进入主图前必须恢复可复现模拟，否则放入补充材料或 discussion。 |

## 补充模拟

补充模拟规划放在 `sup/README.md`。这些内容不进入主图，除非它们改变主结论。

优先补充项：

- 隐藏节点或未观测调控因子；
- off-target 或 noisy perturbation input；
- correlated perturbations 与 partial perturbability；
- wrong basis dictionary 或 over-complete dictionary；
- multiple steady states / failed convergence；
- total-biomass CV sensitivity；
- 更大 Monte Carlo 重复数、失败率与 NA 处理规则。

## 真实数据案例

真实数据作为后续 case study 单独规划，不在模拟脚本中体现。选择数据集时应优先满足：

- 有明确扰动目标和扰动强度；
- 有绝对丰度、总量校正或低噪声连续状态；
- 每个条件达到稳态或近似稳态；
- 有可解释的先验网络、文献验证或生物学标签。

真实数据分析代码后续可放入 `analysis_script/` 或单独 case-study 脚本；数据来源说明放入 `data/`。

## 当前最优先缺口

1. `Fig3a_external_benchmark_main`：外部 benchmark 主实验。
2. `Fig3b_method_capability_matrix`：建立同行方法能力矩阵。
3. `Fig3d/Fig3e`：为现有 scale-free 数值模拟补充 joint vs node-wise 与 hub-recovery 正式面板。
4. `Fig2a_highdim_sample_complexity`：增加重复数，并补充 AUPRC/AUROC 与失败率后冻结 Fig2。

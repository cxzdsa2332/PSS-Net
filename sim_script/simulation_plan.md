# PSS-Net 模拟研究规划

本文档用于把正式模拟脚本组织成论文级模拟路线。当前主文建议围绕三张多面板图展开：图 1 证明 PSS-Net 的可识别性与基础恢复能力，并在视觉上明确真实生成机制是 ODE 动力学、PSS 数据只是稳态切片；图 2 证明样本复杂度与扰动设计价值；图 3 证明鲁棒性、结构依赖与外部 benchmark。模型错设与敏感性分析放入 `sup/`；真实数据作为后续 case study，不混入模拟脚本。

## 总体原则

主文模拟需要回答三个审稿人会直接追问的问题：

1. 扰动稳态（perturbed steady state, PSS）数据能识别什么，不能识别什么？
2. PSS-Net 是否能在高维稀疏网络中随样本预算增加而稳定恢复有向耦合网络？
3. 与现有网络推断方法相比，PSS-Net 使用扰动输入 `u` 能带来什么额外信息？

所有正式模拟脚本应只生成数值结果，写入 `results/sim_results/`。主文图形由 `analysis_script/` 读取结果后生成。`sim_script/manual/` 仍保留探索性脚本，可以混合绘图和检查。

## 图 1：基础可识别性与核心估计器

目标：先展示 PSS-Net 的数据来源是 ODE 动力系统，而不是静态相关样本；再说明在扰动后达到稳态时，ODE 可转化为 PSS-Net 使用的稳态约束。图 1 中可用加性 ODE 与 gLV ODE 作为两个示例：真实系统沿轨迹演化，扰动改变最终稳态，实际测量值是每个扰动条件下的 `x*` 与已知 `u`。在此基础上，再展示 PSS-Net 恢复稳态耦合方程、可识别边界和 ADSIHT 的基础恢复能力。

对应文件夹：`sim_script/01_foundation_recovery/`

模拟数据简述：主要使用可控稀疏 ODE 系统生成扰动稳态数据；观测对象为每个扰动条件下的稳态测量值 `x*` 与已知扰动输入 `u`。图 1 的动力学展示包含加性 ODE `dx/dt = F(x) + u` 和乘性 gLV ODE 两个示例。内部基线比较使用同一批稳态样本、同一候选函数库和同一真实网络。可识别性面板展示 PSS 能识别稳态函数形状，但不能单独恢复完整瞬态机制。

视觉规则：Fig1 当前以 `sim_script/01_foundation_recovery/Fig1.R` 为准。脚本使用 `standard0.R` 的 8 节点历史参数，生成加性线性 ODE 与乘性 gLV ODE 的轨迹；baseline 从统一初始状态出发，single-node input 与 mixed input 从各自模型的 baseline steady state 出发。稳态采样用黑色垂直虚线标注，加性 ODE 与 gLV ODE 使用不同浅色背景，节点名直接标在线末端，不使用 ggplot 图例。Fig1b 不加入机制示意，只展示稳态函数形状与定量 BIC 支持。

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig1a_ode_to_pss_measurement` | 展示 ODE 动力学如何在扰动下产生 PSS 测量值 `(x*, u)`。当前为 2 行 x 3 列：上排 `Additive ODE`，下排 `gLV ODE`；列为 `baseline u = 0`、`single-node input`、`mixed input`。single-node 与 mixed input 从 baseline steady state 出发；黑色垂直虚线为 PSS 采样时间 `t = 20`，终点圆点为 `x*`。 | `sim_script/01_foundation_recovery/Fig1.R` 中 `Fig1a`、`Fig1a_legend`。历史探索图仍保留在 `sim_script/manual/plot_pss_concept.R`，但不作为当前 Fig1a 代码依据。 | 初版完成。脚本只生成 R 图对象，暂不保存输出；后续用于 patchwork 拼图。 |
| `Fig1b_steady_state_function_identifiability` | 展示 PSS 不能区分等价线性机制，但能区分真正非线性的稳态函数结构。当前左侧为代表性边的稳态函数形状，右侧为加入二次基函数后的 BIC gain。 | `sim_script/01_foundation_recovery/Fig1.R` 中 `Fig1b_function_shape`、`Fig1b_identifiability`、`Fig1b`。其中 `Fig1b <- Fig1b_function_shape | Fig1b_identifiability`（若 `patchwork` 可用）。相关批量验证：`pss_net_glv_ss.R` 输出 `results/sim_results/glv_ss_verification.csv`。 | 初版完成。当前模拟中 linear additive 与 standard gLV 的 BIC gain 为负，nonlinear additive 为正；脚本只生成 R 图对象，暂不保存输出。 |
| `Fig1c_adsiht_vs_internal_baselines` | 比较 ADSIHT 与 group lasso 等内部结构化稀疏基线；建议旁边放小型组稀疏/设计矩阵示意。 | 模拟：`sim_script/01_foundation_recovery/pss_net_compare.R`。结果：`results/sim_results/mcc_comparison.csv`。现有汇总：`analysis_script/summarize_mcc_comparison.R`。建议绘图：`analysis_script/Fig1c_adsiht_vs_internal_baselines.R`。 | 数值模拟与表格汇总已有，主图绘图缺失。后续补 AUPR/AUROC。 |
| `Fig1d_identifiability_boundary` | 展示 PSS 可识别稳态函数形状，但不能仅凭稳态区分乘性 gLV 与加性线性动力学；可画 linear/quadratic/saturating 边函数与 GOF。 | 模拟：`sim_script/01_foundation_recovery/pss_net_discriminate.R`。结果：`results/sim_results/discriminate_gof.csv`。建议绘图：`analysis_script/Fig1d_identifiability_boundary.R`。 | 数值模拟已有，绘图缺失。需要与 Fig1b 避免信息重复。 |
| `Fig1e_representative_network_recovery` | 固定一个 seed，展示真实网络与推断网络的 TP、FP、FN、方向、符号和核心节点。 | 建议数据脚本：`sim_script/01_foundation_recovery/Fig1e_representative_network_recovery_data.R`。建议绘图：`analysis_script/Fig1e_representative_network_recovery.R`。 | 缺失。需要保存真实邻接、估计邻接、边权、边符号和节点度数。 |

## 图 2：样本复杂度与扰动设计

目标：证明 PSS-Net 在高维稀疏设定中有合理样本复杂度，并展示主动扰动设计能提高信息效率。

对应文件夹：`sim_script/02_scaling_design/`

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig2a_highdim_sample_complexity` | 展示网络恢复是否随重标度样本量 `N / (s log p)` 改善，不同维度 `p` 的曲线是否近似坍缩。 | 模拟：`sim_script/02_scaling_design/pss_net_highdim.R`。结果：`results/sim_results/highdim_recovery.csv`。建议绘图：`analysis_script/Fig2a_highdim_sample_complexity.R`。 | 数值模拟已有，绘图缺失。后续增加重复数、AUPR/AUROC 和失败率。 |
| `Fig2b_design_linear` | 在线性稳态系统中比较 random、maximin、D-optimal 扰动设计的信息效率。 | 模拟：`sim_script/02_scaling_design/pss_net_design.R`。结果：`results/sim_results/design_comparison.csv`。现有绘图：`analysis_script/plot_design_curves.R`。建议正式面板：`analysis_script/Fig2b_design_linear.R`。 | 数值模拟与通用绘图已有；建议拆出正式面板脚本。 |
| `Fig2c_design_nonlinear` | 非线性互作下展示扰动空间均匀采样不等于特征空间信息最优。 | 模拟：`sim_script/02_scaling_design/pss_net_design_nl.R`。结果：`results/sim_results/design_nl_comparison.csv`。现有绘图：`analysis_script/plot_design_curves.R`。建议正式面板：`analysis_script/Fig2c_design_nonlinear.R`。 | 数值模拟与通用绘图已有；建议拆出正式面板脚本。 |
| `Fig2d_design_strong_nonlinear` | 强非线性互作下展示线性化 D-optimal 打分的优势和局限。 | 模拟：`sim_script/02_scaling_design/pss_net_design_nl_seq.R`。结果：`results/sim_results/design_nl_seq_comparison.csv`。现有绘图：`analysis_script/plot_design_curves.R`。建议正式面板：`analysis_script/Fig2d_design_strong_nonlinear.R`。 | 数值模拟与通用绘图已有；建议拆出正式面板脚本。 |
| `Fig2e_oracle_vs_estimated_design` | 区分 oracle D-optimal 与 realistic pilot-estimated D-optimal，评估可行序贯设计损失。 | 建议模拟：`sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R`。建议绘图：`analysis_script/Fig2e_oracle_vs_estimated_design.R`。 | 缺失。需要包含 oracle、pilot-estimated、maximin、random。 |
| `Fig2f_targeted_design_tradeoff` | 若先验知道 hub，评估集中扰动 hub 是否提高 hub 相关边恢复、是否牺牲整体网络。 | 现有结果：`results/sim_results/targeted_design.csv`。结论记录：`note/targeted_perturbation_tradeoff.md`。建议恢复模拟：`sim_script/02_scaling_design/Fig2f_targeted_design_tradeoff.R`。建议绘图：`analysis_script/Fig2f_targeted_design_tradeoff.R`。 | 结果文件存在但正式模拟代码缺失；进入主图前必须恢复脚本，否则放补充或 discussion。 |

## 图 3：鲁棒性、网络结构与外部 benchmark

目标：把 PSS-Net 放入 bioinformatics/network inference 同行语境中比较，说明其优势来自使用扰动输入和稳态方程，同时诚实展示组成型数据和网络结构依赖。

对应文件夹：`sim_script/03_robustness_benchmarks/`

| 子图 | 要做什么 | 对应文件 / 对象 | 当前状态 |
|------|----------|-----------------|----------|
| `Fig3a_external_benchmark_main` | 在匹配模拟 PSS 数据上，比较 PSS-Net 与外部网络推断方法。最低方法集：PSS-Net/ADSIHT、lasso/elastic net、group lasso、STLSQ/SINDy-style、GENIE3/GRNBoost2、correlation/graphical lasso。 | 参考：`ref/external_benchmark_methods.md`。建议模拟：`sim_script/03_robustness_benchmarks/Fig3a_external_benchmark_main.R`。建议绘图：`analysis_script/Fig3a_external_benchmark_main.R`。 | 缺失，优先级最高。指标应包含 MCC、AUPR/AUROC、precision、recall、sign accuracy、runtime、failure rate。 |
| `Fig3b_method_capability_matrix` | 建立方法能力矩阵：输入类型、输出网络、是否有方向、是否估计非线性边函数、是否需要时序、是否处理 compositional data。 | 参考：`ref/external_benchmark_methods.md`。建议绘图：`analysis_script/Fig3b_method_capability_matrix.R`；可先用手工 CSV/tibble。 | 缺失。适合作为 Fig3 的解释性面板。 |
| `Fig3c_compositional_data_limitation` | 展示 absolute abundance、relative abundance、CLR、relative abundance times noisy total 对恢复的影响。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_compositional.R`。结果：`results/sim_results/compositional_recovery.csv`。建议绘图：`analysis_script/Fig3c_compositional_data_limitation.R`。 | 数值模拟已有，绘图缺失。后续补 total-biomass CV sweep。 |
| `Fig3d_scalefree_joint_vs_nodewise` | 在 scale-free / hub-like 网络中比较 joint block-diagonal estimation 与 node-wise estimation。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`。结果：`results/sim_results/scalefree_compare.csv`。建议绘图：`analysis_script/Fig3d_scalefree_joint_vs_nodewise.R`。 | 数值模拟已有，绘图缺失。当前重复数偏少，主文使用前需增加 R 并明确 NA/失败处理。 |
| `Fig3e_scalefree_hub_recovery` | 评估 joint estimation 是否更稳定恢复核心 hub 排名。指标包括 estimated out-degree vs true out-degree Spearman、top-k hit rate。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`。结果：`results/sim_results/scalefree_compare.csv`。建议绘图：`analysis_script/Fig3e_scalefree_hub_recovery.R`。 | 数值模拟已有，绘图缺失。可与 Fig3d 共用同一数据。 |
| `Fig3f_uniform_negative_control` | 证明 joint estimation 的优势依赖 hub-like 结构；在均匀网络下主要是 precision/recall trade-off。 | 模拟：`sim_script/03_robustness_benchmarks/pss_net_joint_smalln.R`。结果：`results/sim_results/joint_smalln.csv`。建议绘图：`analysis_script/Fig3f_uniform_negative_control.R`。 | 数值模拟已有，绘图缺失。 |

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
2. `Fig2e_oracle_vs_estimated_design`：区分 oracle D-optimal 和 realistic pilot-estimated D-optimal。
3. `Fig1e_representative_network_recovery`：保存并绘制代表性网络恢复图。
4. `Fig3b_method_capability_matrix`：建立同行方法能力矩阵。

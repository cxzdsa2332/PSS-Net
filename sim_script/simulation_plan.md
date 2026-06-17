# PSS-Net 模拟研究规划

本文档用于把正式模拟脚本组织成论文级模拟路线。当前主文建议围绕三张多面板图展开：图 1 证明 PSS-Net 的可识别性与基础恢复能力，图 2 证明样本复杂度与扰动设计价值，图 3 证明鲁棒性、结构依赖与外部 benchmark。模型错设与敏感性分析放入 `sup/`；真实数据作为后续 case study，不混入模拟脚本。

## 总体原则

主文模拟需要回答三个审稿人会直接追问的问题：

1. 扰动稳态（perturbed steady state, PSS）数据能识别什么，不能识别什么？
2. PSS-Net 是否能在高维稀疏网络中随样本预算增加而稳定恢复有向耦合网络？
3. 与现有网络推断方法相比，PSS-Net 使用扰动输入 `u` 能带来什么额外信息？

所有正式模拟脚本应只生成数值结果，写入 `results/sim_results/`。主文图形由 `analysis_script/` 读取结果后生成。`sim_script/manual/` 仍保留探索性脚本，可以混合绘图和检查。

## 图 1：基础可识别性与核心估计器

目标：证明 PSS-Net 恢复的是稳态耦合方程，而不是完整瞬态动力学；同时展示 ADSIHT 相对当前内部基线的优势。

对应文件夹：`sim_script/01_foundation_recovery/`

### Fig1a_concept_workflow

- 子图问题：PSS 数据如何从“扰动输入 + 稳态响应”转化为可解释的有向耦合网络？
- 当前代码：`sim_script/manual/plot_pss_concept.R`
- 当前输出：`results/figure/pss_concept_diagram.pdf`、`results/figure/pss_concept_diagram.png`
- 状态：已有探索版。若进入正式主图，建议迁移或重写为 `analysis_script/Fig1a_concept_workflow.R`。

### Fig1b_glv_steady_state_equivalence

- 子图问题：乘性 gLV 在正稳态下是否严格满足 PSS-Net 使用的加性稳态方程？
- 当前代码：`sim_script/01_foundation_recovery/pss_net_glv_ss.R`
- 当前结果：`results/sim_results/glv_ss_verification.csv`
- 建议正式绘图代码：`analysis_script/Fig1b_glv_steady_state_equivalence.R`
- 状态：数值模拟已有，绘图代码缺失。

### Fig1c_adsiht_vs_internal_baselines

- 子图问题：双稀疏 ADSIHT 是否优于 group lasso 等内部结构化稀疏基线？
- 当前代码：`sim_script/01_foundation_recovery/pss_net_compare.R`
- 当前结果：`results/sim_results/mcc_comparison.csv`
- 现有汇总代码：`analysis_script/summarize_mcc_comparison.R`
- 建议正式绘图代码：`analysis_script/Fig1c_adsiht_vs_internal_baselines.R`
- 状态：数值模拟与表格汇总已有，主图绘图代码缺失。
- 后续增强：加入 rank-based AUPR/AUROC；当前结果主要是阈值化 MCC、precision、recall、F1、JacRMSE。

### Fig1d_identifiability_boundary

- 子图问题：PSS 稳态数据能否区分乘性 gLV 与加性线性模型？能否检测非线性边函数？
- 当前代码：`sim_script/01_foundation_recovery/pss_net_discriminate.R`
- 当前结果：`results/sim_results/discriminate_gof.csv`
- 建议正式绘图代码：`analysis_script/Fig1d_identifiability_boundary.R`
- 状态：数值模拟已有，绘图代码缺失。
- 重点解释：PSS 可识别稳态函数形状，但不能仅凭稳态区分乘性 gLV 与加性线性动力学。

### Fig1e_representative_network_recovery

- 子图问题：在一个代表性重复中，真实网络和 PSS-Net 推断网络的 TP、FP、FN、方向、符号如何对应？
- 当前代码：缺失。
- 建议模拟/提取代码前缀：`sim_script/01_foundation_recovery/Fig1e_representative_network_recovery_data.R`
- 建议绘图代码：`analysis_script/Fig1e_representative_network_recovery.R`
- 状态：缺失。需要从 `pss_net_compare.R` 或 `pss_net_glv_ss.R` 中保存单次 seed 的真实邻接、估计邻接、边权和符号。

## 图 2：样本复杂度与扰动设计

目标：证明 PSS-Net 在高维稀疏设定中有合理样本复杂度，并展示主动扰动设计能提高信息效率。

对应文件夹：`sim_script/02_scaling_design/`

### Fig2a_highdim_sample_complexity

- 子图问题：网络恢复是否随重标度样本量 `N / (s log p)` 改善，不同维度 p 的曲线是否有近似坍缩？
- 当前代码：`sim_script/02_scaling_design/pss_net_highdim.R`
- 当前结果：`results/sim_results/highdim_recovery.csv`
- 建议正式绘图代码：`analysis_script/Fig2a_highdim_sample_complexity.R`
- 状态：数值模拟已有，绘图代码缺失。
- 后续增强：增加 R、补 AUPR/AUROC，并记录失败率。

### Fig2b_design_linear

- 子图问题：线性稳态系统中，random、maximin、D-optimal 扰动设计在固定预算 N 下谁更有效？
- 当前代码：`sim_script/02_scaling_design/pss_net_design.R`
- 当前结果：`results/sim_results/design_comparison.csv`
- 现有绘图代码：`analysis_script/plot_design_curves.R`
- 建议正式绘图代码：`analysis_script/Fig2b_design_linear.R`
- 状态：数值模拟与通用绘图已有；建议拆出正式面板脚本或在统一 Fig2 脚本中调用。

### Fig2c_design_nonlinear

- 子图问题：非线性互作下，为什么在扰动空间均匀采样不等于在特征空间信息最优？
- 当前代码：`sim_script/02_scaling_design/pss_net_design_nl.R`
- 当前结果：`results/sim_results/design_nl_comparison.csv`
- 现有绘图代码：`analysis_script/plot_design_curves.R`
- 建议正式绘图代码：`analysis_script/Fig2c_design_nonlinear.R`
- 状态：数值模拟与通用绘图已有；建议拆出正式面板脚本或在统一 Fig2 脚本中调用。

### Fig2d_design_strong_nonlinear

- 子图问题：强非线性互作下，线性化 D-optimal 打分的优势和局限在哪里？
- 当前代码：`sim_script/02_scaling_design/pss_net_design_nl_seq.R`
- 当前结果：`results/sim_results/design_nl_seq_comparison.csv`
- 现有绘图代码：`analysis_script/plot_design_curves.R`
- 建议正式绘图代码：`analysis_script/Fig2d_design_strong_nonlinear.R`
- 状态：数值模拟与通用绘图已有；建议拆出正式面板脚本或在统一 Fig2 脚本中调用。

### Fig2e_oracle_vs_estimated_design

- 子图问题：当前 D-optimal 是否依赖 oracle 真模型？真实可行的 pilot-estimated D-optimal 与 oracle D-optimal 差多少？
- 当前代码：缺失。
- 建议模拟代码前缀：`sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R`
- 建议绘图代码：`analysis_script/Fig2e_oracle_vs_estimated_design.R`
- 状态：缺失。
- 必要设计：
  - oracle D-optimal：使用真实稳态映射或真实局部 Jacobian，作为上界；
  - pilot-estimated D-optimal：先随机施加少量扰动，估计 surrogate/Jacobian，再序贯更新设计；
  - maximin 与 random：非自适应对照。

### Fig2f_targeted_design_tradeoff

- 子图问题：若先验知道 hub，集中扰动 hub 是否提升 hub 相关边恢复？是否牺牲整体网络？
- 当前结果：`results/sim_results/targeted_design.csv`
- 当前代码：缺失，原探索脚本已删除；结论记录在 `note/targeted_perturbation_tradeoff.md`
- 建议恢复代码前缀：`sim_script/02_scaling_design/Fig2f_targeted_design_tradeoff.R`
- 建议绘图代码：`analysis_script/Fig2f_targeted_design_tradeoff.R`
- 状态：结果文件存在，但正式模拟代码缺失。进入主图前必须恢复脚本；否则放补充或 discussion。

## 图 3：鲁棒性、网络结构与外部 benchmark

目标：把 PSS-Net 放入 bioinformatics/network inference 同行语境中比较，说明其优势来自使用扰动输入和稳态方程，同时诚实展示组成型数据和网络结构依赖。

对应文件夹：`sim_script/03_robustness_benchmarks/`

### Fig3a_external_benchmark_main

- 子图问题：在匹配模拟 PSS 数据上，PSS-Net 相比外部网络推断方法表现如何？
- 当前代码：缺失。
- 参考方法清单：`ref/external_benchmark_methods.md`
- 建议模拟代码前缀：`sim_script/03_robustness_benchmarks/Fig3a_external_benchmark_main.R`
- 建议绘图代码：`analysis_script/Fig3a_external_benchmark_main.R`
- 状态：缺失，优先级最高。
- 最低方法集：
  - PSS-Net / ADSIHT；
  - lasso 或 elastic net，同一稳态方程；
  - group lasso，同一稳态方程；
  - STLSQ/SINDy-style library regression；
  - GENIE3 或 GRNBoost2，state-only directed black-box；
  - correlation 或 graphical lasso，association baseline。
- 指标：MCC、AUPR/AUROC、precision、recall、sign accuracy、runtime、failure rate。

### Fig3b_method_capability_matrix

- 子图问题：不同方法使用什么输入、输出什么网络、是否有方向、是否能估计非线性边函数、是否需要时序？
- 当前代码：缺失。
- 参考文档：`ref/external_benchmark_methods.md`
- 建议绘图代码：`analysis_script/Fig3b_method_capability_matrix.R`
- 状态：缺失。可先用手工 CSV 或 tibble 构建。

### Fig3c_compositional_data_limitation

- 子图问题：绝对丰度、相对丰度、CLR、相对丰度乘含噪总量对 PSS-Net 网络恢复有什么影响？
- 当前代码：`sim_script/03_robustness_benchmarks/pss_net_compositional.R`
- 当前结果：`results/sim_results/compositional_recovery.csv`
- 建议正式绘图代码：`analysis_script/Fig3c_compositional_data_limitation.R`
- 状态：数值模拟已有，绘图代码缺失。
- 后续增强：补 total-biomass CV 扫描，建议放 supplement 或作为 Fig3c 的 inset。

### Fig3d_scalefree_joint_vs_nodewise

- 子图问题：scale-free / hub-like 网络中，joint block-diagonal estimation 是否比 node-wise 更好？
- 当前代码：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`
- 当前结果：`results/sim_results/scalefree_compare.csv`
- 建议正式绘图代码：`analysis_script/Fig3d_scalefree_joint_vs_nodewise.R`
- 状态：数值模拟已有，绘图代码缺失。
- 注意：当前 R 偏少，且曾出现 NA；主文使用前需要增加重复数、明确失败处理。

### Fig3e_scalefree_hub_recovery

- 子图问题：joint estimation 是否更稳定地恢复核心 hub 排名？
- 当前代码：`sim_script/03_robustness_benchmarks/pss_net_scalefree.R`
- 当前结果：`results/sim_results/scalefree_compare.csv`
- 建议正式绘图代码：`analysis_script/Fig3e_scalefree_hub_recovery.R`
- 状态：数值模拟已有，绘图代码缺失。
- 指标：estimated out-degree vs true out-degree Spearman、top-k hit rate。

### Fig3f_uniform_negative_control

- 子图问题：joint estimation 的优势是否依赖 hub-like 结构？在均匀网络下是否只是 precision/recall trade-off？
- 当前代码：`sim_script/03_robustness_benchmarks/pss_net_joint_smalln.R`
- 当前结果：`results/sim_results/joint_smalln.csv`
- 建议正式绘图代码：`analysis_script/Fig3f_uniform_negative_control.R`
- 状态：数值模拟已有，绘图代码缺失。

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

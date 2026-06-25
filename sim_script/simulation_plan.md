# PSS-Net 模拟研究规划与当前实现

本文档以 `sim_script/01_foundation_recovery/`、`sim_script/02_scaling_design/` 和 `sim_script/03_robustness_benchmarks/` 中的现行脚本为准，记录 Fig1--Fig3 **目前实际执行的模拟、读取的结果和生成的绘图对象**。规划内容与已完成内容分开标注，避免把“已有脚本”“已有结果文件”和“主图已经可复现”混为一谈。

## 当前代码组织

- Fig1 不是完全分离的“模拟脚本 + 绘图脚本”：`Fig1.R` 内联运行 Fig1a、Fig1b、Fig1d--Fig1f 的模拟，只从 CSV 读取 Fig1c，并在末尾直接写出根目录下的 `Fig1.pdf`。
- Fig2 主要从现有 CSV 读取数值结果；Fig2b 的概念设计由 `Fig2.R` 通过 `sys.source()` 运行 `Fig2_explain.R`。`Fig2.R` 只生成工作区对象，不写图文件。
- Fig3 从三份 CSV 读取数值结果，另在 `Fig3.R` 内生成能力矩阵和拓扑示意。脚本可生成 `Fig3a`--`Fig3e`，但目前没有组装总图对象 `Fig3`，也不写图文件。
- 当前主线 benchmark 使用二阶单项式库 `M_ord = 2`；Fig1g 另以多种 1--3 项函数库做错设敏感性分析。PSS-Net 以 source node 为组进行 node-wise ADSIHT/DSIC 拟合。

## Fig1--Fig3 主叙事主线

整套模拟图的主线应从“为什么 PSS 数据可以用于网络恢复”逐步推进到“什么时候恢复更好、什么时候会失效”，避免把每个子图写成彼此独立的技术展示。

1. **Fig1：建立问题与方法可信度。** PSS-Net 并不直接拟合完整 ODE 轨迹，而是利用扰动后稳态 `(x*, u)` 中保留下来的调控函数信息。Fig1 的核心任务是证明：在合适的稳态可识别条件、噪声水平和函数库下，PSS-Net 能从 PSS 测量中恢复 directed/signed sparse network，并能恢复边函数形状；同时承认可识别边界、噪声和字典错设会限制恢复。
2. **Fig2：说明实验设计为什么重要。** 在同样的 PSS-Net 推断框架下，样本预算大致跟随 sparse scale `N/(s log p)`，但扰动输入 `u` 的选择会改变有效信息量。Fig2 的主叙事应写成：**designed inputs improve recovery especially under nonlinear / low-budget regimes**。D-optimal 只是现成实验设计工具之一，不是本文算法创新；当前 30 次重复结果也显示 maximin 与 pilot/oracle D-optimal 的表现依赖 regime 和预算，因此不应声称 D-optimal 必然最优。
3. **Fig3：把方法放回真实比较与限制中。** Fig3 的核心不是单纯“赢过所有方法”，而是定位 PSS-Net 与现有 MRA/稀疏回归/动力学工具的适用边界：PSS-Net 的优势在于匹配 PSS 输入、directed nonlinear edge functions 和边排序能力；但阈值化 MCC、拓扑结构、joint vs node-wise 假设、组成型测量和绝对尺度缺失都会影响实际表现。

一句话版本：**Fig1 证明 PSS 信息中有可恢复的调控函数；Fig2 证明实验设计能让这些信息更有效，尤其在非线性和小预算下；Fig3 证明该优势需要放在外部方法、拓扑和测量限制中谨慎解释。**

## 图 1：PSS 数据来源、可识别边界与基础恢复

目标：说明真实生成机制是 ODE 动力学，而 PSS-Net 使用的是各扰动条件达到稳态后的 `(x*, u)`；随后展示稳态函数的可识别边界、噪声影响、基础网络恢复、边函数恢复和网络可视化。

主叙事：Fig1 先回答“PSS-Net 到底在利用什么信息”。ODE 轨迹不是直接输入，稳态响应函数才是推断对象；当稳态函数在候选库中可近似、噪声不过高且样本预算足够时，PSS-Net 可以同时恢复 sparse support、edge direction/sign 和局部非线性边函数。Fig1 也需要明确边界：二阶库不是万能函数逼近器，BIC/噪声扫描和 Fig1g 的字典错设结果应作为可信度边界，而不是弱化主结果的附属材料。

主绘图脚本：`sim_script/01_foundation_recovery/Fig1.R`

| 子图 | 当前实际内容 | 数据与参数 | 当前状态 |
|------|--------------|------------|----------|
| `Fig1a_ode_to_pss_measurement` | 对 additive ODE 与 multiplicative gLV 分别绘制 baseline、single-node input、mixed input 的动态轨迹，并在终点标记 PSS 测量。 | 历史 8 节点系统；模拟直接写在 `Fig1.R` 中，不读取结果文件。 | 已实现；对象为 `Fig1a`、`Fig1a_legend`。 |
| `Fig1b_steady_state_function_identifiability` | 比较 additive linear、standard gLV 和 additive nonlinear 三类机制的稳态函数形状与线性/二次基 BIC；另做 SNR 扫描，量化二次基被选中的概率。 | `N_id = 180`；基础示例噪声 `sigma_id = 0.003`；SNR 为 `100, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5, 1`，每档 30 次。模拟直接写在 `Fig1.R`。 | 已实现；对象为 `Fig1b_function_shape`、`Fig1b_identifiability`、`Fig1b_noise_snr`、`Fig1b` 和 `snr_summary`。 |
| `Fig1c_adsiht_vs_group_lasso_scaling` | 比较 node-wise ADSIHT 与 group lasso 在 linear/nonlinear truth 下的 MCC、AUPRC、AUROC、`CoefL2` 和 `JacRMSE`。 | 独立模拟脚本 `Fig1c_adsiht_group_lasso_scaling.R`；`p = 8, 30, 100`，对应 `s_in = 2, 3, 3`；`N/(s log p) = 4, 8, 12, 16`；`M_ord = 2`，`R = 30`，绘图固定 `SNR = 30`。结果 CSV 现有 1,440 行且无 NA。 | 模拟结果和绘图均存在；对象为 `Fig1c_method_scaling`、`Fig1c`。 |
| `Fig1d_effect_decomposition` | 在固定 10 节点非线性加性系统中，将目标节点的 self effect 与 received regulation 分解，并比较 true 与 ADSIHT 恢复。完整动态 overlay 代码仍由 `if (FALSE)` 禁用。 | `p = 10`、`N = 200`、`SNR = 30`、`M_ord = 2`；与 Fig1e/Fig1f 共用一次内联拟合。 | 已实现；对象为 `Fig1d`。 |
| `Fig1e_function_shape_recovery` | 展示 self feedback 与 cross-node edge function 的 true/estimated 曲线。 | 使用 Fig1d 同一 10 节点系统和拟合。当前主图不读取已有的 `Fig1e_function_shape_recovery.csv`。 | 已实现；对象为 `Fig1e`。 |
| `Fig1f_true_vs_inferred_network` | 用 igraph 并排展示真实和推断的有向带符号网络，并区分 TP、FP、FN。 | 使用 Fig1d/Fig1e 同一系统、seed 和拟合结果。 | 已实现；对象为 `Fig1f`。 |
| `Fig1g_basis_robustness` | 交叉比较 quadratic、Monod、sine truth 与 linear/poly2/poly3/Monod/Fourier 拟合库，区分支持恢复稳健性与边函数形状恢复对字典匹配的依赖。 | `Fig1x_basis_misspecification.R`；`p = 20`、`s_in = 2`、`SNR = 30`、`R = 30`，单一充足预算；结果为 `Fig1x_basis_misspecification.csv`，现有 450 行且无 NA。 | 已按 30 次重复重跑；对象为 `Fig1g_basis_robustness`、`Fig1g`。 |

`Fig1.R` 在安装 `cowplot` 时组装 A4 纵向总图 `Fig1`，并调用 `ggsave("Fig1.pdf", ...)`。当前根目录已有重绘后的 `Fig1.pdf`；本轮另导出 `results/figure/Fig1.png` 与 `results/figure/Fig1.pdf` 方便预览和归档。

## 图 2：样本复杂度与扰动设计

目标：展示 PSS-Net 随 `N/(s log p)` 增加的恢复趋势，以及 random、maximin、oracle D-optimal 和 pilot-estimated D-optimal 扰动设计之间的机制与性能差异。

主叙事：Fig2 应聚焦“同样数量的 PSS 实验，怎么选扰动输入更有信息”。Fig2a 给出稀疏恢复的预算尺度；Fig2b 解释 random、space-filling 和模型辅助设计在输入空间中的差异；Fig2c--d 量化 designed inputs 相对 random 的恢复收益和预算节省；Fig2e 说明 pilot-estimated D-optimal 是把现有 D-optimal 设计后端接到 PSS 响应映射上的可行 two-stage 实现。当前主结论建议写成：**designed inputs improve recovery especially under nonlinear / low-budget regimes**。不要把 D-optimal 写成 PSS-Net 的算法创新，也不要声称 D-optimal 总是优于 maximin 或 random；30 次重复下，线性高预算处 random 可与设计法持平或略高。

主绘图脚本：`sim_script/02_scaling_design/Fig2.R`

D-optimal 选点调用 `AlgDesign::optFederov()`；这是现有最优实验设计后端，不作为 PSS-Net 的算法创新。Fig2 各面板当前只保留标题，说明文字放在 caption 中，不再使用副标题。

| 子图 | 当前实际内容 | 数据与参数 | 当前状态 |
|------|--------------|------------|----------|
| `Fig2a_highdim_sample_complexity` | 绘制 node-wise ADSIHT 的 MCC 随 `N/(s log p)` 变化，并标注 `MCC = 0.8` 的插值预算。当前没有在主面板展示 AUPRC/AUROC。 | `Fig2a_highdim_sample_complexity.R`；`p = 8, 50, 100`，`s = 2, 3, 3`；12 档预算；`M_ord = 2`、`R = 30`、`sigma = 0.03`。CSV 现有 1,080 行且无 NA。 | 已按 30 次重复重跑；绘图代码存在。 |
| `Fig2b_design_mechanism` | 在 3D 扰动空间展示 random、maximin、oracle D-optimal 和 pilot-estimated D-optimal。四种策略共享 8 个 pilot，再选 20 个条件。 | `Fig2_explain.R`；候选池 3,000，设计总数 28；oracle 使用玩具非线性稳态映射，pilot 版本用含噪 pilot 经 ridge 估计局部映射。 | 绘图代码存在；当前环境已可调用 `AlgDesign`，对象目标为 `Fig2b_design_mechanism`、`Fig2b`。 |
| `Fig2c_design_gain_map` | 对 linear、nonlinear、strong nonlinear 三类系统，计算 maximin/D-optimal 相对 random 的配对 `Delta MCC` 和胜率，并画成 regime × budget 热图。 | linear：`Fig2b_design_linear.R`，`p = 8`、`s = 2`、`R = 30`、`SNR = 30`、`N = 7, 9, 11, 13, 17, 25, 34`；两类 nonlinear：`pss_net_design_nl.R` 和 `pss_net_design_nl_seq.R`，均为 `p = 8`、`s = 2`、`R = 30`、`sigma = 0.04`、`N = 12, 16, 20, 30, 40, 60`。三份 CSV 均存在且无 NA。 | 已按 30 次重复重跑；图中的 D-optimal 是使用真实响应映射的 oracle 设计。30 次重复下，线性高预算处 random 可与设计法持平或略高，主结论应聚焦 nonlinear/early-budget 的设计收益。 |
| `Fig2d_budget_to_target` | 将 Fig2c 的均值学习曲线插值为达到 `MCC = 0.5/0.6` 所需预算；超出扫描上限时标记 `> max`。 | 使用与 Fig2c 相同的三份 CSV，不运行额外模拟。 | 绘图代码存在；对象为 `Fig2d_budget_to_target`、`Fig2d`。 |
| `Fig2e_oracle_vs_estimated_design` | 在共享 pilot 且 pilot 计入总预算的条件下，比较 random、maximin、oracle D-optimal 和 pilot-estimated D-optimal；可视化 oracle regret 与相对 random 的增益。 | `Fig2e_oracle_vs_estimated_design.R`；`p = 8`、`M_ord = 2`、`R = 30`、`sigma = 0.04`、候选池 2,500；`pilot_n = 8, 12, 16`，`N_total = 20, 30, 40, 60`。CSV 现有 1,440 行且无 NA。 | 已按 30 次重复重跑；当前结果显示 pilot D-opt 相对 random 多数为正，但 oracle regret 有正有负，说明 D-optimal 更应作为现成设计工具/敏感性分析，而非保证最优的创新点。 |

当前 `Fig2.R` 的真实可运行状态：脚本语法正确，当前环境已可完整 source 并生成 `Fig2a`--`Fig2e` 和总图 `Fig2`。本轮已导出 `results/figure/Fig2.png` 与 `results/figure/Fig2.pdf`。

## 图 3：外部 benchmark、拓扑依赖与组成型限制

目标：先用能力矩阵定位 PSS-Net，再比较外部方法的边检测与边排序，随后展示 scale-free/ER 拓扑差异、joint 与 node-wise 恢复的结构依赖，以及组成型测量的局限。

主叙事：Fig3 应把 Fig1--Fig2 的内部模拟结果放进更真实的比较语境。Fig3a 先声明不同方法的输入假设和输出能力不同，避免把不完全同类的方法强行排成单一榜单；Fig3b 用 MCC/AUPRC 展示 PSS-Net 在匹配 PSS 设定下具有强边排序和可竞争的支持恢复，但阈值定义会影响 MCC；Fig3c--d 说明拓扑结构会改变 node-wise 与 joint 策略的相对优势；Fig3e 则提醒组成型或缺少绝对总量的测量会破坏稳态函数恢复。Fig3 的写法应强调“适用边界清晰的优势”，而不是“无条件优于所有基线”。

主绘图脚本：`sim_script/03_robustness_benchmarks/Fig3.R`

### 当前主图结构

| 子图 | 当前实际内容 | 数据与参数 | 当前状态 |
|------|--------------|------------|----------|
| `Fig3a_method_capability_matrix` | 9 种方法 × 7 项能力：匹配的 PSS 输入、directed、signed、nonlinear edge functions、built-in sparse selection、scale-free/hub applicability、compositional aware。MRA 在当前表中记为支持 signed edges。 | 内容在 `Fig3.R` 中手工定义，依据 `ref/external_benchmark_methods.md`。 | 已实现；对象为 `Fig3a_method_capability_matrix`、`Fig3a`。 |
| `Fig3b_benchmark_edge_recovery` | 在代表性预算下用 seed 点和 mean ± SD 比较 MCC 与 AUPRC。MCC 表示默认支持集恢复；AUPRC 表示连续边分数的排序能力。当前主面板不画 AUROC，也不画 Function NRMSE。 | `Fig3b_external_benchmark_main.R`；`p = 30`、`s_in = 2`、`M_ord = 2`、`R = 30`、`SNR = 30`；linear/strong nonlinear；`N = 69, 103, 171, 273`；8 种方法。PSS-Net 与 aiMeRA 均以 row-normalized `abs(link) > 0.05` 判边；PSS-Net 仍保留绝对系数/Jacobian 用于强度误差。现有 CSV 1,888 行；strong nonlinear 的 `N = 69` 有 4 个 seed 因有效稳态不足跳过。 | 已按 30 次重复重跑。全预算平均：linear 下 PSS-Net/aiMeRA 的 MCC 和 AUPRC 均为 1.000；strong nonlinear 下 PSS-Net/aiMeRA 的 MCC 为 0.963/0.948，AUPRC 为 0.988/0.979。 |
| `Fig3c_topology_schematic` | 用两个小型有向图直观对比 scale-free hub topology 与同质 ER topology，节点大小和填色表示 out-degree。 | 直接在 `Fig3.R` 生成；示意图 `p = 18`、`avg_in = 2`，scale-free 示意使用 `pa_power = 3`。 | 已实现；对象为 `Fig3c_topology_schematic`、`Fig3c`。 |
| `Fig3d_structure_results` | 将相同 seed、topology、预算下的 node-wise 与 joint PSS-Net 画成配对散点；横轴 node-wise，纵轴 joint，`y > x` 表示 joint 更好。指标为 edge MCC、out-degree Spearman rho 和 top-k hub hit rate。 | `Fig3c_structure_dependence.R`；`p = 50`、`avg_in = 2`、`M_ord = 2`、`R = 30`、`sigma = 0.03`；`N = 16, 24, 32, 40`；scale-free 模拟实际使用 `pa_power = 2`，并与固定入度 ER 对照。CSV 现有 240 行且无 NA。 | 已按 30 次重复重跑；对象为 `Fig3d_structure_results`、`Fig3d`。总体上 joint 在 scale-free/hub 指标上更稳定，例如 scale-free `N = 40` 的 MCC node-wise/joint 为 0.484/0.520，hub rho 为 0.529/0.641。 |
| `Fig3e_compositional_data_limitation` | 比较 absolute、relative、CLR 和 relative × noisy total 四种输入下的 Precision、Recall、MCC；absolute 是 oracle 上界。 | 当前模拟脚本名为 `Fig3e_compositional_data_limitation.R`，但输出仍沿用历史文件名 `Fig3c_compositional_data_limitation.csv`；`p = 10`、`N = 200`、`M_ord = 2`、`R = 30`、`sigma_rel = 0.05`、`T_cv = 0.15`。CSV 现有 120 行且无 NA。 | 已按 30 次重复重跑。当前平均 MCC：absolute 0.890、relative 0.520、CLR 0.317、relative × noisy total -0.047；当前 noisy-total 校正没有恢复 oracle，且 MCC 可为负，主图 y 轴后续应允许显示负值。 |

### 边强度与函数恢复诊断

`Fig3b_external_benchmark_main.csv` 同时保存 `EdgeJacRMSE/NRMSE`、线性/非线性系数误差、`FuncRMSE/NRMSE`、sign accuracy、runtime 和 failure 状态。`Fig3.R` 目前额外生成补充对象 `Fig3c_benchmark_recovery_curves`：

- 横轴为全预算 `N/(s log p)`；
- 纵向指标为局部边强度 `EdgeJacNRMSE` 和完整边函数 `FuncNRMSE`；
- 仅纳入能够输出绝对尺度边函数的 PSS-Net、Lasso、Elastic net、Group lasso 和 PySINDy STLSQ；
- MRA 因输出是 row-normalized response 且没有完整边函数，不参与绝对尺度恢复误差比较。

该对象是补充诊断，不是当前 Fig3c 主面板。AUROC 也保留在 CSV 中，但不进入 Fig3b 主面板。

### Fig3 当前可运行状态

`Fig3.R` 已在现有环境完整 source 验证，可生成 `Fig3a`、`Fig3b`、`Fig3c`、`Fig3d`、`Fig3e` 和 `Fig3c_benchmark_recovery_curves`。本轮已临时组装总图对象并导出 `results/figure/Fig3.png` 与 `results/figure/Fig3.pdf`；正式组装代码仍可进一步写回 `Fig3.R`。Fig3e 的结果文件仍保留历史 `Fig3c_...csv` 名称，但脚本名和缺文件提示已统一为 `Fig3e_compositional_data_limitation.R`。

## 补充模拟与历史结果

| 内容 | 当前文件 | 当前状态 |
|------|----------|----------|
| scale-free 单拓扑下 node-wise vs joint | `pss_net_scalefree.R` → `scalefree_compare.csv` | 旧模拟脚本和结果均存在；主图已由同时包含 scale-free/ER 的 `Fig3c_structure_dependence.R` 取代。 |
| 同质网络 small-N joint negative control | `pss_net_joint_smalln.R` → `joint_smalln.csv` | 模拟脚本和结果存在，尚无正式补充绘图脚本。 |
| targeted perturbation trade-off | `results/sim_results/targeted_design.csv`、`note/targeted_perturbation_tradeoff.md` | 结果文件存在（30 行），但仓库中没有能生成该 CSV 的 R 脚本；当前属于不可完全复现的历史结果。 |
| benchmark 扩展指标 | `Fig3b_external_benchmark_main.csv` | AUROC、Jacobian/函数误差、sign accuracy、runtime、failure 已存在；GENIE3/GIES 尚未纳入当前结果。 |

其他计划中的敏感性分析仍包括隐藏节点、off-target/noisy perturbation、correlated perturbations、wrong/over-complete basis、multiple steady states、failed convergence 和 total-biomass CV sweep。

## 真实数据案例

真实数据仍作为后续 case study，不混入当前模拟主图。数据集应优先满足：有明确扰动目标和强度；有绝对丰度、总量校正或低噪声连续状态；条件达到稳态或近似稳态；并有可用于解释的先验网络或验证标签。

## 当前最优先缺口

1. 优化 Fig2c/Fig2e 的最终版式密度，尤其是热图单元文字与图例位置；当前 30 次重复结果已可生成主图。
2. 决定 Fig2 是否在补充材料加入 AUPRC/AUROC、失败率和 bootstrap CI，主图目前仅展示 MCC 与设计收益。
3. 为 `Fig3a`--`Fig3e` 增加正式版式组装对象和输出文件；如需彻底重命名 Fig3e 历史 CSV，应一次性同步脚本与结果。
4. 增加组成型 total-biomass CV sweep，确认 noisy-total 校正在什么误差范围内才有恢复价值。
5. 为 `targeted_design.csv` 恢复可复现的生成脚本，或从正式结果中移除该历史文件。

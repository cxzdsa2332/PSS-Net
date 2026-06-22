# Fig1 主要结果解释草稿

本文档根据当前 `sim_script/01_foundation_recovery/Fig1.R` 与 `results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv` 的结果整理，用于后续撰写 Results 小节。当前版本仍是初步模拟，主文前建议增加重复数并固定最终排版。

## 总体信息

图 1 的核心信息是：PSS-Net 并不是从静态相关样本中直接画网络，而是利用扰动稳态数据 `(x*, u)` 约束 ODE 的稳态方程。该图按逻辑分为三层：第一层说明数据如何由 ODE 动力学产生；第二层说明 PSS 数据能识别稳态函数形状，但不能任意区分所有瞬态机制；第三层展示 ADSIHT 在有限噪声和有限样本下可以恢复函数形状与有向带符号网络。

## Fig1a：ODE 动力学产生 PSS 测量

Fig1a 使用历史 8 节点参数，对比 additive ODE 与 multiplicative gLV ODE。在 baseline、single-node input 和 mixed input 三种扰动下，轨迹先演化到稳态，黑色垂直虚线处的终点即实际用于 PSS-Net 的观测 `x*`。该面板的写作重点是：PSS-Net 的输入不是时间序列，而是多个扰动条件下的稳态切片；但这些稳态切片仍然来自明确的 ODE 生成机制。

## Fig1b：稳态函数形状的可识别边界

Fig1b 表明 PSS 数据可以支持稳态函数形状识别，但需要明确边界。在线性 additive ODE 与标准 gLV ODE 中，二者可导出相同或近似等价的线性稳态约束，因此加入二次基函数并不会得到正的 BIC 支持。相反，在真正非线性的 additive ODE 中，二次基函数获得正的 BIC gain，说明 PSS 数据能够检测稳态函数形状的非线性。

SNR 辅助扫描进一步显示，非线性识别依赖测量质量：低 SNR 时二次项选择率明显下降，约 `SNR < 20` 基本不稳定，`20-30` 为过渡区，`SNR >= 50` 较稳定。这个结果适合写成方法边界，而不是过度宣称 PSS 能恢复完整瞬态机制。

## Fig1c：ADSIHT 与 group lasso 的基础恢复比较（线性 vs 非线性）

Fig1c 在 `p = 8, 30, 100` 下比较 node-wise ADSIHT 与 node-wise group lasso，并在**线性**与**非线性**两类真值系统上各跑一遍。两类真值共享同一组 `A`、`r` 与扰动设计，仅区别在二次项 `B` 是否激活：线性真值用闭式稳态，非线性真值（约半数边带 `b != 0` 的弯曲项）通过积分加性 ODE 得到稳态，因此非线性情形同时考察组稀疏（选哪个源）与组内稀疏（选哪个单项式），与 Fig1d/e/f 的模型类一致。样本预算按 `N / (s log p)` 扫描，观测噪声固定为 `SNR = 30`，`sigma_x = signal_scale / 30`。当前结果文件含 1440 行，即 3 个维度 × 4 个样本预算 × 30 个重复 × 2 类真值 × 2 种方法。

绘图上线性用虚线、非线性用实线，方法用颜色区分（ADSIHT 蓝、group lasso 红）；因四条线叠加，已去掉 ±sd 阴影带，仅保留均值线与点。

按所有样本预算和 30 个重复求平均，主要结果如下：

| p | method | truth | MCC | AUPRC | AUROC | CoefL2 |
|---|--------|-------|-----|-------|-------|--------|
| 8 | ADSIHT | linear | 0.886 | 0.916 | 0.975 | 0.357 |
| 8 | ADSIHT | nonlinear | 0.876 | 0.924 | 0.976 | 0.534 |
| 30 | ADSIHT | linear | 0.912 | 0.956 | 0.997 | 0.255 |
| 30 | ADSIHT | nonlinear | 0.900 | 0.954 | 0.994 | 0.464 |
| 100 | ADSIHT | linear | 0.898 | 0.966 | 0.999 | 0.209 |
| 100 | ADSIHT | nonlinear | 0.892 | 0.958 | 0.999 | 0.414 |
| 8 | GroupLasso | linear | 0.438 | 0.719 | 0.847 | 0.534 |
| 8 | GroupLasso | nonlinear | 0.422 | 0.725 | 0.861 | 0.712 |
| 30 | GroupLasso | linear | 0.786 | 0.919 | 0.976 | 0.343 |
| 30 | GroupLasso | nonlinear | 0.741 | 0.886 | 0.971 | 0.768 |
| 100 | GroupLasso | linear | 0.924 | 0.945 | 0.981 | 0.386 |
| 100 | GroupLasso | nonlinear | 0.907 | 0.912 | 0.973 | 0.529 |

结果解读应保持克制，可归纳为三点：

1. **方法排序是主导效应。** ADSIHT 在 p=8、p=30 的 MCC/AUPRC/AUROC 明显优于 group lasso（如 p=8 MCC 约 0.88 vs 0.42–0.44），且优势在线性与非线性两类真值下都成立；这印证 ADSIHT 的组内稀疏能在弯曲边上保留二次项、在线性边上剔除二次项，而 group lasso 在低维过度选择（precision 仅约 0.45）。
2. **线性与非线性曲线接近。** 引入非线性只带来较小的 MCC/precision 代价与更高的 CoefL2（曲率确实更难拟合），方法间相对排序不变。因此 Fig1c 的样本复杂度结论对模型是否线性是**稳健**的——这是展示两类真值的主要理由，建议在正文/图注中以一句话说明。
3. **高维下两法趋同。** p=100 时 group lasso 的 MCC 略高（更高 precision、更低 recall），ADSIHT 反之（更高 recall、更高 AUPRC/AUROC）；高维二分类边集指标应结合阈值与 precision-recall trade-off 解释，而非单看 MCC。

## Fig1d：推断函数可回到动态效应解释

Fig1d 使用固定 10 节点非线性加性 ODE，展示两个目标节点的 integrated self effect、received regulation 以及二者之和。真实曲线与 ADSIHT 推断曲线使用实线/虚线对照。该面板的重点不是证明 PSS 数据直接提供瞬态轨迹，而是说明：一旦从 PSS 方程估计出可解释的加性函数，模型可以把扰动后状态变化拆解为 self feedback 和 received cross-node regulation，从而形成动态层面的解释。

## Fig1e：稳态函数形状恢复

Fig1e 展示同一 10 节点系统中，ADSIHT 对代表性 self feedback 和 cross-node effect functions 的恢复。真实函数与推断函数在多个源-靶对上大体重合，包括线性边和带二次项的非线性边。该面板支持 Fig1b 的正向部分：在有足够扰动覆盖和合理 SNR 时，PSS 数据不仅能恢复边是否存在，还能恢复边函数的稳态形状。

## Fig1f：真实网络与推断网络

Fig1f 将同一 10 节点系统的真实有向网络与 ADSIHT 推断网络并排展示。边方向为 source 到 target；促进作用用红色，抑制作用用蓝色；推断网络中 TP、FP 和 FN 通过线型与灰色辅助显示。该面板把 Fig1d/Fig1e 的函数恢复结果汇总为网络层面的视觉证据，说明 PSS-Net 可以输出有方向、有符号、可追溯到函数形状的调控网络。

## 建议写作口径

图 1 的主结论可以写成：PSS-Net converts perturbed steady-state measurements into interpretable constraints on an underlying ODE system, recovers steady-state coupling functions under identifiable settings, and yields directed signed networks with competitive sparse recovery across dimensions.

需要避免的过度表述：

- 不应说 PSS 数据单独恢复完整瞬态机制；
- 不应把 Fig1c 写成 ADSIHT 在所有指标和所有维度上绝对优于 group lasso；
- Fig1c 现为 `R = 30` 重复、线性与非线性两类真值，可作为稳定趋势引用；其余面板若进入主文仍需对齐重复数；
- Fig1d-Fig1f 是单 seed 展示性结果，应作为代表性机制可视化，而非总体性能统计。

# Fig1 主要结果解释草稿

本文档根据当前 `sim_script/01_foundation_recovery/Fig1.R` 与 `results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv` 的结果整理，用于后续撰写 Results 小节。当前版本仍是初步模拟，主文前建议增加重复数并固定最终排版。

## 总体信息

图 1 的核心信息是：PSS-Net 并不是从静态相关样本中直接画网络，而是利用扰动稳态数据 `(x*, u)` 约束 ODE 的稳态方程。该图按逻辑分为三层：第一层说明数据如何由 ODE 动力学产生；第二层说明 PSS 数据能识别稳态函数形状，但不能任意区分所有瞬态机制；第三层展示 ADSIHT 在有限噪声和有限样本下可以恢复函数形状与有向带符号网络。

## Fig1a：ODE 动力学产生 PSS 测量

Fig1a 使用历史 8 节点参数，对比 additive ODE 与 multiplicative gLV ODE。在 baseline、single-node input 和 mixed input 三种扰动下，轨迹先演化到稳态，黑色垂直虚线处的终点即实际用于 PSS-Net 的观测 `x*`。该面板的写作重点是：PSS-Net 的输入不是时间序列，而是多个扰动条件下的稳态切片；但这些稳态切片仍然来自明确的 ODE 生成机制。

## Fig1b：稳态函数形状的可识别边界

Fig1b 表明 PSS 数据可以支持稳态函数形状识别，但需要明确边界。在线性 additive ODE 与标准 gLV ODE 中，二者可导出相同或近似等价的线性稳态约束，因此加入二次基函数并不会得到正的 BIC 支持。相反，在真正非线性的 additive ODE 中，二次基函数获得正的 BIC gain，说明 PSS 数据能够检测稳态函数形状的非线性。

SNR 辅助扫描进一步显示，非线性识别依赖测量质量：低 SNR 时二次项选择率明显下降，约 `SNR < 20` 基本不稳定，`20-30` 为过渡区，`SNR >= 50` 较稳定。这个结果适合写成方法边界，而不是过度宣称 PSS 能恢复完整瞬态机制。

## Fig1c：ADSIHT 与 group lasso 的基础恢复比较

Fig1c 在 `p = 8, 30, 100` 下比较 node-wise ADSIHT 与 node-wise group lasso。样本预算按 `N / (s log p)` 扫描，观测噪声固定为 `SNR = 30`，其中 `sigma_x = signal_scale / 30`。当前结果文件包含 72 行，即 3 个维度、4 个样本预算、3 个重复、2 种方法。

按所有样本预算和 3 个重复求平均，主要结果如下：

| p | method | MCC | AUPRC | AUROC | Precision | Recall | CoefL2 | JacRMSE |
|---|--------|-----|-------|-------|-----------|--------|--------|---------|
| 8 | ADSIHT | 0.882 | 0.894 | 0.968 | 0.851 | 0.990 | 0.282 | 0.114 |
| 8 | GroupLasso | 0.494 | 0.731 | 0.854 | 0.517 | 1.000 | 0.566 | 0.175 |
| 30 | ADSIHT | 0.915 | 0.973 | 0.998 | 0.855 | 1.000 | 0.247 | 0.060 |
| 30 | GroupLasso | 0.781 | 0.936 | 0.980 | 0.678 | 0.968 | 0.318 | 0.051 |
| 100 | ADSIHT | 0.902 | 0.961 | 0.999 | 0.821 | 0.999 | 0.207 | 0.026 |
| 100 | GroupLasso | 0.922 | 0.944 | 0.981 | 0.893 | 0.961 | 0.382 | 0.035 |

结果解读应保持克制：ADSIHT 在 p=8 和 p=30 的 MCC、AUPRC、AUROC 明显优于 group lasso；在 p=100 下 group lasso 的 MCC 略高，但 ADSIHT 仍有更高的 AUPRC/AUROC 和更低的 CoefL2。当前结论更适合表述为：ADSIHT 在基础稀疏 PSS 方程中表现出稳定的排序与参数恢复优势，但高维二分类边集指标需要结合阈值、precision-recall trade-off 和更多重复数解释。

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
- 当前 `R = 3` 只能作为初版趋势，主文前需要增加 Monte Carlo 重复数；
- Fig1d-Fig1f 是单 seed 展示性结果，应作为代表性机制可视化，而非总体性能统计。

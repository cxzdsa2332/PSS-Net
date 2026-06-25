# Fig1 解释：PSS 数据来源、可识别边界与基础恢复

## 主叙事

Fig1 要回答的是：PSS-Net 到底从什么信息中恢复网络。核心叙事是，PSS-Net 并不直接拟合完整 ODE 时间轨迹，而是利用多个扰动条件达到稳态后的 `(x*, u)`，把这些稳态切片转化为对调控函数的约束。在合适的可识别条件、噪声水平和函数库下，这些稳态信息足以恢复 sparse support、edge direction/sign，并给出可解释的局部边函数形状。

这张图也要明确边界：PSS 数据不能无条件恢复完整瞬态机制；二阶多项式库不是万能函数逼近器；噪声、稳态覆盖范围和字典错设都会影响函数形状恢复。

## Panel a：ODE dynamics generate perturbed steady-state measurements

Fig1a 展示 additive ODE 与 gLV ODE 在 baseline、single-node input 和 mixed input 下如何从动态轨迹收敛到稳态。黑色虚线处的终点才是 PSS-Net 使用的观测。

该 panel 的重点是：

- PSS-Net 的输入不是时间序列，而是扰动后的稳态响应；
- 这些稳态响应仍来自明确的 ODE 生成机制；
- 不同扰动 `u` 改变系统的稳态位置，从而提供恢复调控函数的信息。

## Panel b：steady-state function shape has an identifiable boundary

Fig1b 对比 additive linear、standard gLV 和 additive nonlinear 三类机制的稳态函数形状与 BIC 选择结果。线性 additive ODE 与 standard gLV 在稳态方程层面可能产生近似线性的约束，因此二次基函数不一定得到支持；真正 additive nonlinear ODE 才会让二次项获得稳定的 BIC gain。

SNR 扫描说明非线性识别依赖测量质量。当前结果显示低 SNR 时二次项选择率不稳定，约 `SNR = 20–30` 是过渡区，更高 SNR 下非线性稳态签名更容易被检测到。

写作时应强调：Fig1b 证明的是“稳态函数形状在一定条件下可识别”，而不是“PSS 数据能区分所有可能的动态机制”。

## Panel c：ADSIHT vs group lasso scaling

Fig1c 比较 node-wise ADSIHT 与 group lasso 在 linear/nonlinear truth 下的网络恢复。设置为：

- `p = 8, 30, 100`
- `s_in = 2, 3, 3`
- `N/(s log p) = 4, 8, 12, 16`
- `SNR = 30`
- `R = 30`
- CSV：`results/sim_results/Fig1c_adsiht_group_lasso_scaling.csv`
- 当前结果：1440 行，无 NA

按所有预算平均，ADSIHT 的 MCC 大致为：

| p | linear MCC | nonlinear MCC |
|---|---:|---:|
| 8 | 0.886 | 0.876 |
| 30 | 0.912 | 0.900 |
| 100 | 0.898 | 0.892 |

Group lasso 的平均 MCC 为：

| p | linear MCC | nonlinear MCC |
|---|---:|---:|
| 8 | 0.438 | 0.422 |
| 30 | 0.786 | 0.741 |
| 100 | 0.924 | 0.907 |

解读上有三点：

1. ADSIHT 在低维和中等维度下明显优于 group lasso，尤其 p=8 和 p=30。
2. linear 与 nonlinear truth 的曲线接近，说明 ADSIHT 的 sparse recovery 对是否存在二次边函数相对稳健。
3. p=100 时 group lasso 的 MCC 可略高，但 ADSIHT 的 AUPRC/AUROC 和 recall 仍很强；因此高维处不要只用 MCC 做绝对优劣判断。

## Panel d：effect decomposition links inferred functions to dynamics

Fig1d 在固定 10 节点非线性加性系统中，将目标节点状态变化分解为 self effect 与 received regulation。真实曲线和 ADSIHT 推断曲线相互对照。

该 panel 是代表性机制可视化，不是 Monte Carlo 统计结果。它说明：一旦 PSS-Net 从稳态方程中估计出加性调控函数，就可以把网络边重新解释为动态贡献，而不是只输出一张无方向相关图。

## Panel e：steady-state edge-function recovery

Fig1e 展示代表性 self feedback 和 cross-node edge functions 的 true vs estimated 曲线。它支持 Fig1b 和 Fig1d 的正向部分：在扰动覆盖充分、SNR 合理且函数库匹配时，PSS-Net 不只恢复边是否存在，还能恢复边函数的局部形状。

## Panel f：true vs inferred directed signed network

Fig1f 将同一 10 节点系统的真实网络与推断网络并排展示。促进边和抑制边用颜色区分，TP/FP/FN 用线型区分。

该 panel 的作用是把 Fig1d/e 的函数恢复汇总到网络层面，展示 PSS-Net 输出的是 directed、signed、function-backed network。

## Panel g：basis robustness and misspecification

Fig1g 检查 quadratic、Monod、sine truth 在 linear/poly2/poly3/Monod/Fourier 等拟合库下的表现。设置为：

- `p = 20`
- `s_in = 2`
- `SNR = 30`
- `R = 30`
- 单一充足预算
- CSV：`results/sim_results/Fig1x_basis_misspecification.csv`
- 当前结果：450 行，无 NA

主要结论是：

- support recovery 对 fitted library 的错设相对稳健，MCC 多数约在 `0.79–0.92`；
- edge-function NRMSE 明显依赖字典匹配；
- 例如 poly2 truth 下 poly2 library 的 function NRMSE 最低，Monod truth 下 Monod library 的 function NRMSE 最低。

因此 Fig1g 的主旨应写成：support recovery can be robust to moderate library misspecification, but accurate edge-function recovery benefits from a matched dictionary.

## 建议写作口径

推荐主结论：

> PSS-Net converts perturbed steady-state measurements into constraints on underlying regulatory functions, enabling recovery of directed signed sparse networks and local edge-function shapes when the steady-state response is identifiable and the function library is adequate.

需要避免：

- 不说 PSS 数据单独恢复完整 ODE 瞬态机制；
- 不说二阶库能表示所有非线性；
- 不把 Fig1d/e/f 的单 seed 示例写成总体统计；
- 不把 ADSIHT 描述成在所有维度和所有指标上都绝对优于 group lasso。

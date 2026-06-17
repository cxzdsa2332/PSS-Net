# Note: Fig1c ADSIHT vs group lasso scaling benchmark

**用途**：记录 Fig1c 的维度-样本预算模拟。
**日期**：2026-06-17。相关脚本：
[../sim_script/01_foundation_recovery/Fig1c_adsiht_group_lasso_scaling.R](../sim_script/01_foundation_recovery/Fig1c_adsiht_group_lasso_scaling.R)、
[../sim_script/01_foundation_recovery/Fig1.R](../sim_script/01_foundation_recovery/Fig1.R)。

## 设计原则

Fig1c 先只比较 **node-wise ADSIHT** 与 **node-wise group lasso**，维度取 `p = 8, 30, 100`。
暂时不加入 external benchmark、joint block-diagonal estimation 或 scale-free network，原因是这些因素会把
"基础估计器差异"、"网络拓扑异质性"和"外部方法输入不同"混在一起。当前面板只回答一个主问题：
在普通稀疏网络和相同 PSS 稳态方程下，ADSIHT 是否比 group lasso 更稳定，且这种优势是否随 p 增大仍保留。

模拟使用闭式稳态

$$
(\mathrm{diag}(\gamma)-A)x^* = r + u,
$$

对应加性线性 ODE 和正稳态 gLV 的共享 PSS 方程。每个节点固定少量入边，使用单项式库 `M=2`。
样本预算按 `N / (s log p)` 扫描。观测噪声不再按固定标准差添加，而是按信噪比控制：

$$
\sigma_x = \frac{\mathrm{signal\_scale}}{\mathrm{SNR}},
\qquad
\mathrm{signal\_scale} = \frac{1}{p}\sum_{j=1}^p \mathrm{sd}(x^*_{\cdot j}).
$$

当前固定使用中等噪声 `SNR = 30`，展示维度和样本预算效应。指标包括 MCC、AUPRC、AUROC、precision、recall、sign accuracy、runtime，以及参数恢复误差 CoefL2、JacRMSE、EdgeWeightRMSE。

## Fig1c 初步结果：SNR = 30（3 seeds）

| p | method | MCC | AUPRC | AUROC | CoefL2 | JacRMSE | runtime/sec |
|---|--------|-----|-------|-------|--------|---------|-------------|
| 8 | ADSIHT | 0.882 | 0.894 | 0.968 | 0.282 | 0.114 | 0.003 |
| 8 | GroupLasso | 0.494 | 0.731 | 0.854 | 0.566 | 0.175 | 0.011 |
| 30 | ADSIHT | 0.915 | 0.973 | 0.998 | 0.247 | 0.060 | 0.031 |
| 30 | GroupLasso | 0.781 | 0.936 | 0.980 | 0.318 | 0.051 | 0.097 |
| 100 | ADSIHT | 0.902 | 0.961 | 0.999 | 0.207 | 0.026 | 0.542 |
| 100 | GroupLasso | 0.922 | 0.944 | 0.981 | 0.382 | 0.035 | 0.915 |

## 解释

在 `SNR = 30` 的主展示层，ADSIHT 在 p=8 和 p=30 的 MCC、AUPRC、AUROC 均高于 group lasso；p=100 时 group lasso 的 MCC 略高，但 ADSIHT 的 AUROC、AUPRC 和参数 L2 仍较好。这个结果提示：当前初版不能简单写成“ADSIHT 全面优于 group lasso”，更合理的表述是 ADSIHT 在较低维和中等维度下给出更稳健的网络识别，在高维高稀疏设定中 group lasso 也可能通过强收缩获得较好的二分类边集表现。

参数误差指标仍需要谨慎解读：高维极稀疏矩阵中，大量真零边被强收缩到零会降低整体 L2/RMSE。因此参数误差需要和 AUPRC/MCC、非零边权误差、sign accuracy 一起解释，不能单独作为网络恢复优劣的结论。

这个结果符合方法预期：
ADSIHT 同时利用 source-level group sparsity 和 within-edge basis sparsity；group lasso 只做组选择，
在 `M=2` 库下更容易保留不必要的组内项或在低预算下出现 precision/recall trade-off。

当前结果只应作为 Fig1c 初版：`R=3` 偏少，主文前需要增加重复数，并固定 group lasso 的调参规则
（当前为每个节点沿 grpreg path 用 BIC 选 lambda）。scale-free topology、joint estimation 和 external benchmark
应放到后续 Fig3 或专门 robustness benchmark 中分析。

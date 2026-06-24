# Fig3b：PSS-Net 判边阈值是需要重点调整的方法学环节

记录日期：2026-06-23

> 状态说明（2026-06-24）：主 benchmark 保持 `M_ord=2`；PSS-Net 与 aiMeRA 的主分析
> 现统一使用 row-normalized `abs(link)` 分数及固定 `> 0.05` 判边。minimum-DSIC
> `A_out` 保留为 native-selection 敏感性结果，不再作为 headline MCC。

## 发现

在当前 Fig3b benchmark（`p=30`、每节点入度 `s_in=2`、边密度 6.9%、
`N=69,103,171,273`、`R=10`、`SNR=30`）中，PSS-Net 的连续边强度排序很好，
但原来的二值判边规则把“ADSIHT 选中的任意非零组”都算作边。这会保留一批强度很小的
假阳性，明显压低 MCC。

| PSS-Net 判边规则 | linear MCC | strong-nonlinear MCC |
|---|---:|---:|
| 原始：组范数 `>1e-8` | 0.900 | 0.881 |
| 绝对 Jacobian：`abs(J)>0.05` | 1.000 | 0.983 |
| 按拟合 self-Jacobian 行归一化后：`abs(link)>0.05` | 1.000 | 0.965 |
| aiMeRA：`abs(link)>0.05` | 1.000 | 0.951 |

与 aiMeRA 尺度最接近的是第三行。其 strong-nonlinear MCC 在四个预算下分别为
0.961、0.965、0.967、0.967。PSS-Net 原始 AUPRC 已达到 0.997（linear）和
0.986（strong-nonlinear），进一步支持“问题主要在二值化，而非边排序失败”。

数值来源：

- `results/sim_results/Fig3b_external_benchmark_main.csv`
- `results/sim_results/Fig3b_pss_threshold_sensitivity.csv`
- `sim_script/03_robustness_benchmarks/Fig3b_external_benchmark_main.R`

## 当前解释

1. 这不是新的估计器改进：系数、Jacobian 和 FuncRMSE 都没有变化，只改变了连续 score
   到二值边集的映射。
2. 直接对 `abs(J)` 使用 0.05 与 aiMeRA 不同尺度，不适合做 headline 比较。若要共用
   0.05，应先用拟合 self-Jacobian 做行归一化，使对角尺度对应 -1。
3. 当前结果说明原始 `>1e-8` 规则过松，但不能据此把 0.05 当成普适最优阈值。0.05 是
   从 aiMeRA 分析规则继承的事先给定值；在看到结果后固定它仍有 post-hoc 风险。
4. 当前真边系数与零边分离较明显、噪声较低，可能使 0.05 特别有效。这个现象必须在
   弱边、更高噪声和不同状态尺度下复核。

## 后续重点调整（优先级顺序）

1. **统一 score 定义。** PSS-Net 主分析使用参考稳态处的 `abs(J)`；跨 MRA 比较时另报
   row-normalized link。不要再用“任意非零组”直接作为最终边集。
2. **把 AUPRC 设为阈值无关主指标。** MCC 作为预先定义阈值下的辅助指标，避免算法各自
   使用有利阈值造成不公平。
3. **预注册阈值选择。** 比较固定 0.05、训练集/内层交叉验证、稳定选择和 bootstrap-FDR；
   阈值不得使用测试 truth 或同一批测试 seed 调整。
4. **绘制阈值路径。** 对 `tau=0.005--0.20` 展示 MCC、precision、recall 和预测边数，
   同时画真边/零边 score 分布，判断 0.05 是宽平台还是偶然峰值。
5. **做尺度与难度敏感性分析。** 至少扫描 SNR、弱边比例、`p/s_in`、扰动幅度、状态单位
   重标度及正/负曲率混合；检查 row-normalization 在 self-Jacobian 接近零时是否不稳定。
6. **正文措辞保持克制。** 在上述验证完成前，只能表述为“预设阈值敏感性分析显示
   PSS-Net 的低 MCC 主要来自微小假阳性”，不能宣称 PSS-Net 普遍优于 aiMeRA。

## 建议图形

将该发现放入 Supplement：左面板为 MCC--threshold 曲线，右面板为真边与零边的
row-normalized `abs(link)` 分布；主文 Fig3b 同时报 AUPRC 和预注册阈值下的 MCC。

# Pilot-informed D-optimal design：文献边界与 Fig2e 实现映射

## 结论先行

`pilot D-optimal` 不是一个具有唯一标准定义的算法名称。Fig2e 当前方法更准确的名称是
**pilot-informed locally D-optimal design** 或 **two-stage plug-in D-optimal design**：先用一批随机
pilot PSS 数据估计扰动到稳态的局部响应映射，再把该估计代入 D-optimal 准则选择后续扰动。

它属于以下方法谱系的组合，但**不是对其中任何一篇论文算法的逐字复现**：

1. Chernoff (1953) 的 locally optimal design：非线性模型的最优设计依赖未知参数，只能在一个
   给定/估计参数点处局部最优；
2. Kiefer--Wolfowitz (1960) 的 D/G 等价思想：D-optimal 与控制最大预测方差之间的联系；
3. Fedorov exchange algorithm 与 CRAN `AlgDesign::optFederov()`：在有限候选集中做带受保护
   pilot runs 的 exact D-optimal augmentation；这是当前代码实际调用的设计后端；
4. Wynn (1970) 与 adaptive Wynn 文献：提供序贯/自适应 D-optimal 的理论谱系，但当前代码
   不再自行实现 Wynn 贪心选点；
5. ridge regression：稳定小 pilot 下的局部响应映射估计。

完整书目信息在 [`references.bib`](references.bib) 中，键为
`chernoff1953locally`、`kiefer1960equivalence`、`wynn1970sequential`、
`wheeler2025algdesign`、`fedorov1972theory`、`hoerl1970ridge`、
`ford1992canonical`、`chaloner1995bayesian`、
`freise2021adaptivewynn` 和 `freise2024pstep`。

## 核实后的核心论文

| 文献 | 经核实的主要内容 | 与 Fig2e 的关系 | 不应声称的内容 |
|---|---|---|---|
| Chernoff (1953), DOI: `10.1214/aoms/1177728915` | 参数依赖的 locally optimal design 奠基工作。 | 支持“用 pilot 估计代替未知真参数后做局部最优设计”的基本逻辑。 | 没有提出 PSS 网络或当前 ridge 响应映射。 |
| Kiefer & Wolfowitz (1960), DOI: `10.4153/CJM-1960-030-4` | 经典 D-optimal/G-optimal 等价结果。 | 解释 D-optimal 与控制最大预测方差之间的联系。 | 不等于当前有限候选池 exact augmentation 的直接性能定理。 |
| Fedorov (1972); Wheeler (2025), `AlgDesign` 1.2.1.2, DOI: `10.32614/CRAN.package.AlgDesign` | 有限候选集上的 exchange algorithm；`optFederov(..., augment=TRUE)` 可保护已有设计点并增广 exact D-optimal design。 | 当前 Fig2b/Fig2e 的实际选点后端；PSS 代码只负责生成候选特征矩阵。 | 包不估计 PSS 响应映射，也不提供 PSS-Net 网络恢复理论。 |
| Wynn (1970), DOI: `10.1214/aoms/1177696809` | 逐点生成 D-optimum designs。 | 是 sequential D-optimal 的经典谱系与未来在线版本的依据。 | 当前正式 Fig2e 已改用 Fedorov exchange 包实现，不再声称复现 Wynn 算法。 |
| Ford, Torsney & Wu (1992), DOI: `10.1111/j.2517-6161.1992.tb01897.x` | 非线性问题中 locally optimal design 的构造。 | 支持 oracle/plug-in 局部设计的非线性背景。 | 不提供 Fig2e 的具体网络推断流程。 |
| Chaloner & Verdinelli (1995), DOI: `10.1214/ss/1177009939` | Bayesian experimental design 综述。 | 提醒：若 pilot 参数不确定性很大，应积分参数不确定性，而不只是代入一个点估计。 | Fig2e 当前不是 Bayesian design。 |
| Freise, Gaffke & Schwabe (2021), DOI: `10.1214/20-AOS1974` | adaptive Wynn：利用迄今观测反复估参并选择下一点；在其模型与正则条件下研究一致性和渐近局部 D-optimality。 | 是未来完整闭环版本最直接的理论谱系。 | 其理论假设是特定非线性/GLM 单响应模型，不能直接移植为 PSS-Net 网络恢复保证。 |
| Freise, Gaffke & Schwabe (2024), DOI: `10.1007/s00362-023-01502-4` | 每阶段加入参数维数个点的 p-step-ahead adaptive D-optimal 非线性回归算法。 | 与“分批实验、每批后重估”的未来 Fig2e 扩展最接近。 | 当前 Fig2e 并未实现每批重新估计响应映射。 |
| Hoerl & Kennard (1970), DOI: `10.1080/00401706.1970.10488634` | ridge regression。 | 对应 pilot 样本少时，以岭项稳定 `dx*/du` 的矩阵估计。 | ridge 本身不提供 D-optimal 设计。 |

## Fig2e 当前实现

脚本：[`../sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R`](../sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R)

### 1. 真值系统与观测

模拟 8 节点非线性加性 ODE：

\[
\dot x_j=\mu_j-\gamma_jx_j+\sum_{i\ne j}\{A_{ji}x_i+B_{ji}x_i^2\}+u_j.
\]

每节点有 2 条随机入边，其中至少 1 条含二次项。对每个扰动 `u` 用 `lsoda` 积分到
`t = 2000`，以末状态近似 `x*(u)`，再加入标准差 `sigma = 0.04` 的独立加性高斯测量噪声。

### 2. 公平预算

- 共同随机 pilot：`N_pilot = 8, 12, 16`；第一条为 `u = 0`；
- 总预算：`N_total = 20, 30, 40, 60`；
- 后续可用条件数严格为 `N_adaptive = N_total - N_pilot`；
- 四种策略共享同一真值系统、pilot、2500 个候选扰动和候选点对应的噪声实现。

因此 pilot 不是主动方法获得的“免费额外数据”。

### 3. Pilot-estimated 响应映射

对中心化 pilot 数据拟合

\[
X_c\approx U_cH,
\qquad
\widehat H=(U_c^\top U_c+0.25I)^{-1}U_c^\top X_c.
\]

再预测候选稳态：

\[
\widehat x^*(u)=\bar x+(u-\bar u)\widehat H.
\]

这里 `H` 是经验性的局部 `dx*/du`，不是完整 ODE 参数或网络 Jacobian 的直接估计。

### 4. Oracle 响应映射

oracle 使用真值系统在未扰动稳态的 Jacobian `J_true`，以

\[
x^*_{\mathrm{oracle}}(u)\approx x^{\mathrm{wt}}-J_{\mathrm{true}}^{-1}u
\]

预测候选稳态。它仍是**真 Jacobian 下的局部线性 oracle**，而不是提前计算每个候选扰动的
精确非线性稳态。

### 5. D-optimal 特征与包实现

对预测稳态构造

\[
\phi(x)=(x_1,x_1^2,\ldots,x_8,x_8^2).
\]

候选特征先按列中心化、标准化；若精确稳态关系导致确定性线性依赖，则用 pivoted QR 保留
可估计的满秩列空间。数据矩阵不手工添加常数列，但 `AlgDesign` 公式 `~ .` 会按标准模型加入
一个截距；因此 16 个满秩 PSS 特征对应的 exact-design 模型秩为 17。

随后调用：

```r
AlgDesign::optFederov(
  frml = ~ .,
  data = as.data.frame(rbind(Phi_pilot, Phi_pool)),
  nTrials = N_total,
  criterion = "D",
  augment = TRUE,
  rows = seq_len(N_pilot),
  maxIteration = 100,
  nRepeats = 1
)
```

`augment=TRUE` 使 pilot rows 不参与交换；包使用 Fedorov exchange algorithm 从候选池选择其余
`N_total-N_pilot` 个点。每个总预算单独优化，因此不同预算的设计不强制互为嵌套前缀。

### 6. 四个对照

1. `random`：从共同候选池随机延续；
2. `maximin`：在原始扰动 `u` 空间最大化到已选点的最小距离；
3. `oracle_dopt`：真 Jacobian 局部映射 + D-optimal；
4. `pilot_dopt`：含噪 pilot ridge 局部映射 + D-optimal。

### 7. 网络恢复与评价

真实稳态模拟完成后，以 PSS 方程把 `-u_j` 回归到所有 `(x_i,x_i^2)`，用 ADSIHT 选择源节点组，
再用 Precision、Recall 和 MCC 比较推断边集与真网络。`R = 20` 个随机系统共生成 960 行结果。

## 当前实现与完整 adaptive Wynn 的关键差异

当前 Fig2e 是**一次 pilot 后固定响应映射的包增广设计**：

1. 用 pilot 得到一次 `H_hat`；
2. 按 `H_hat` 预测所有候选特征；
3. 将预测特征交给 `AlgDesign::optFederov()`，对每个总预算计算 exact augmentation；
4. 新选点的真实响应不会被用于重新估计 `H_hat`。

而 Freise et al. (2021, 2024) 所讨论的 adaptive 方法会在获得新响应后重新估参，再决定后续点。
因此论文中应写：

> We use a pilot-informed, locally D-optimal perturbation design based on a
> ridge-estimated local steady-state response map.

不应直接写：

> We implement the adaptive Wynn algorithm of Freise et al.

除非后续代码真正加入“观测一批 -> 重估响应映射 -> 再选下一批”的循环。

## 后续可扩展方向

1. **完整批量 adaptive Wynn**：每获得一批真实 PSS 后重估 `H` 或 PSS-Net，再重新打分；
2. **Bayesian D-optimal**：对响应映射/Jacobian 后验积分，避免 plug-in 设计忽略 pilot 不确定性；
3. **稳健/maximin local design**：在多个可能 Jacobian 上优化最坏情况，而非仅使用一个点估计；
4. **稀疏目标设计准则**：从全参数 `log det` 改为更直接服务边支持恢复的 coherence、restricted
   eigenvalue 或目标节点加权准则。

# PSS-Net: Sparse Additive Nonparametric ODE Inference from Perturbed Steady States

**方法说明文档** · PSS-Net v3 · 对应实现：`script/sindy_ss_v3.R`

---

## 1. 问题设定

### 1.1 复杂系统动力学模型

考虑由 $p$ 个相互耦合的状态变量构成的复杂系统。状态变量可表示组分浓度、活性水平、负载、库存、细胞状态比例或其他非负系统量；其动力学由如下耦合常微分方程（ODE）系统刻画：

$$\frac{dx_j(t)}{dt} = \mu_j + f_{jj}\!\left(x_j(t)\right) + \sum_{i \neq j} f_{ji}\!\left(x_i(t)\right) + u_j(t), \qquad j = 1, \ldots, p \tag{1}$$

其中：
- $x_j(t) \geq 0$ 为节点/变量 $j$ 在时刻 $t$ 的状态；
- $\mu_j \in \mathbb{R}$ 为节点 $j$ 的基线漂移或内禀趋势（baseline rate）；
- $f_{ji} : \mathbb{R}_{\geq 0} \to \mathbb{R}$ 为**未知**的单变量效应函数；$i\neq j$ 时描述变量 $i$ 对变量 $j$ 的非线性耦合效应，$i=j$ 时描述变量 $j$ 对自身的反馈/调节效应；
- $u_j(t)$ 为施加于节点 $j$ 的外部扰动或控制输入。

**示例：广义 Lotka–Volterra（gLV）系统。** gLV 可作为非线性复杂系统中的一个典型 benchmark。其动力学是**乘性**的：节点 $j$ 的变化率正比于自身状态 $x_j$ 乘以一个状态依赖的相对增长/变化率。为清楚区分"其他变量对 $j$ 的影响"与"$j$ 对自身的影响"，可将二者分开写：

$$\frac{dx_j}{dt} = x_j\Big(r_j + \sum_{i \neq j} A_{ji} x_i - \gamma_j x_j + u_j\Big) \tag{2}$$

其中 $r_j$ 为基线变化率；$A_{ji}$（$i\neq j$）为**跨节点耦合**强度，即变量 $i$ 对 $j$ 的影响；$\gamma_j > 0$ 为**自反馈/自阻尼**强度，对应自影响项 $-\gamma_j x_j$；$u_j$ 为作用于括号内相对变化率的外部扰动。

**与紧凑矩阵写法的关系**：若定义

$$B_{ji}=A_{ji}\ (i\neq j), \qquad B_{jj}=-\gamma_j \tag{2b}$$

则 (2) 可压缩为标准矩阵形式 $\dot x_j=x_j\big(r_j+\sum_{i=1}^p B_{ji}x_i+u_j\big)$，即 $\dot x_j=x_j\big(r_j+(Bx)_j+u_j\big)$。下文默认使用拆开的 $A_{ji}$ 与 $\gamma_j$ 记号，避免把自反馈和跨节点耦合混在同一个求和号里。扰动 $u_j$ 作用于括号内的相对变化率；只有这样除以 $x_j^*$ 才得线性稳态 (2$'$)。若改为括号外的绝对流入 $+u_j$，稳态将含 $u_j/x_j^*$，加性回归不再精确。

展开 (2) 含跨节点耦合的**双线性项** $A_{ji}x_jx_i$（$i\neq j$）与自反馈项 $-\gamma_j x_j^2$，故 gLV 在**动力学层面不属于**式 (1) 的加性类。

**稳态约化：gLV $\xrightarrow{\ \div\, x_j\ }$ 加性线性关系。** PSS 推断只用稳态。在可行（正）稳态 $x_j^* > 0$ 处，(2) 的括号项必为零：

$$r_j + \sum_{i \neq j} A_{ji} x_i^* - \gamma_j x_j^* + u_j = 0 \tag{2$'$}$$

这恰是通用加性模型 (1) 在

$$f_{ji}(x_i)=A_{ji}x_i\ (i\neq j), \qquad f_{jj}(x_j)=-\gamma_j x_j, \qquad \mu_j=r_j \tag{2a}$$

下的稳态形式：跨节点耦合与自反馈分别进入加性代理，且**双线性项 $x_jx_i$ 因除以 $x_j$ 被消去**。由此：

- 乘性 gLV (2) 与加性 ODE (1) **共享同一组稳态方程 (2$'$)**，稳态值 $x^*$ 完全相同；
- 当真实系统是 gLV 时，PSS 的稳态回归（§3–§5）仍然**严格有效**，无需显式建模 $x_jx_i$。此即稳态线性化思想，本文将其用于更一般的非参数耦合函数 $f_{ji}$。

**适用范围说明**：式 (1) 是本文的主要建模对象，适用于稳态方程可写成加性单变量效应之和的复杂系统；gLV 只是其中一个可被稳态约化到该形式的示例。在**动力学瞬态**上，加性代理与真 gLV 不同（仅稳态重合），故 §8 的轨迹对比验证的是加性代理而非 gLV 瞬态。若需在动力学上表达真双线性/高阶耦合，须引入成对交互函数 $f_{ij}(x_i,x_j)$ 或更高阶项（后续扩展）。模拟真值提供两种实现：加性线性（`sindy_ss_*.R`）与**真乘性 gLV**（`sim_script/01_foundation_recovery/pss_net_glv_ss.R`，后者仅作为示例 benchmark，用于验证稳态回归在乘性系统上的网络恢复与稳态等价性）。

### 1.2 耦合网络定义

定义有向图 $G = (V, E)$，节点集 $V = \{1, \ldots, p\}$，有向边集

$$E = \{(j \leftarrow i) : f_{ji} \not\equiv 0,\ i\neq j\}$$

即变量 $i$ 对变量 $j$ 存在非零跨节点效应。自反馈 $f_{jj}$ 单独估计，不计入跨节点网络边集。**网络重构**的目标是：从观测数据中估计 $\hat{E}$ 使其尽可能接近真实 $E$。

### 1.3 可识别性约束

为消除截距与效应函数之间的不可识别性（$\mu_j$ 与 $f_{ji}$ 的常数项无法独立估计），施加如下约束：

$$f_{ji}(0) = 0, \qquad \forall j, i \tag{3}$$

**建模含义**：变量 $i$ 处于零基准状态时，对目标变量 $j$ 不产生额外效应；常数偏移统一由 $\mu_j$ 表示。在约束 (3) 下，$\mu_j$ 唯一可识别。在 gLV 示例中，$\mu_j$ 对应基线变化率 $r_j$。

**实现**：选取满足 $\psi_m(0) = 0$（$\forall m$）的基函数族，即单项式基 $\psi_m(x) = x^m$（$m = 1, \ldots, M$）天然满足此约束。

---

## 2. 扰动稳态实验（PSS）

### 2.1 稳态方程

在稳态 $\dot{x}_j = 0$ 时，式 (1) 退化为代数约束：

$$\mu_j + f_{jj}\!\left(x_j^*\right) + \sum_{i \neq j} f_{ji}\!\left(x_i^*\right) + u_j = 0 \tag{4}$$

对乘性 gLV 示例，此式即 §1.1 的 (2$'$)（已除以 $x_j^*>0$）；故下文回归对加性 ODE 与可稳态约化到加性形式的乘性系统通用。

**核心优势**：无需对时序数据进行数值微分。观测到的稳态状态 $x^*$ 直接进入线性回归，避免了导数估计引入的额外噪声。

### 2.2 扰动设计

对系统施加 $N$ 个独立扰动条件 $\{u^{(k)}\}_{k=1}^N$，每次等待系统达到新稳态 $x^{(k)*}$，记录观测对：

$$\mathcal{D} = \left\{\left(\tilde{x}^{(k)},\ u^{(k)}\right)\right\}_{k=1}^{N}$$

其中 $\tilde{x}^{(k)} = x^{(k)*} + \varepsilon^{(k)}$ 为含加性测量噪声的观测值（$\varepsilon^{(k)} \sim \mathcal{N}(0, \sigma^2 I)$）。

**实现细节**（`sindy_ss_v3.R` Step 1）：

- $N_{\text{cond}} = 300$ 个扰动条件；扰动分量 $u_j^{(k)} \sim \text{Uniform}(-0.4,\ 0.8)$，独立采样（在 gLV benchmark 中作用于括号内相对变化率；幅度需保证稳态仍位于有效状态空间）；
- 第一个条件设为无扰动参考状态：$u^{(1)} = \mathbf{0}$；
- 每个条件以随机初值 $x_0 \sim \text{Uniform}(0.5, 2.0)^p$ 积分至 $t = 5000$（`rtol=1e-9`、`atol=1e-11`），取末值作为稳态估计；
- 过滤掉未收敛（含 `NaN`/`Inf` 或负值）的条件，保留有效条件 $N \leq N_{\text{cond}}$；
- 加入测量噪声 $\sigma = 0.03$，以 $x_{\min} = 10^{-6}$ 截断防止零值。

### 2.3 无扰动参考稳态的无噪声估计

Jacobian 评估和动力学重建需要精确的无扰动参考稳态 $x^{\mathrm{ref}}$。直接使用观测值 $\tilde{x}^{(1)}$ 会将 $\sigma = 0.03$ 的测量噪声引入后续所有非线性运算，产生系统性偏差。

**解决方案**：以 $u = \mathbf{0}$ 独立运行高精度 ODE 积分（$t_{\max} = 10^4$，`rtol=1e-12`，`atol=1e-14`），获得无噪声稳态：

$$x^{\mathrm{ref}} = \lim_{t \to \infty} x(t;\, u = \mathbf{0}) \tag{5}$$

$x^{\mathrm{ref}}$ 仅用于 Jacobian 计算（Step 4）和动力学对比（Step 6–7），不参与回归设计矩阵的构建。

---

## 3. 非参数基展开

### 3.1 函数近似

采用 $M$ 阶截断单项式基近似未知效应函数：

$$f_{ji}(a) \approx \sum_{m=1}^{M} \theta_{jim} \cdot \psi_m(a) = \psi(a)^\top \theta_{ji}, \qquad \psi_m(a) = a^m \tag{6}$$

其中 $\theta_{ji} = (\theta_{ji1}, \ldots, \theta_{jiM})^\top \in \mathbb{R}^M$ 为待估系数向量。

**为何选用单项式基**：

1. 自然满足约束 $\psi_m(0) = 0$（式 3），无需额外处理；
2. 对任意线性稳态真值 $f_{ji}(x)=C_{ji}x$，仅第一阶系数非零（$\theta_{ji1}=C_{ji}$，其余高阶系数为 0），真实稀疏结构与基展开完全兼容；gLV 示例对应 $C_{ji}=A_{ji}$（$i\neq j$）与 $C_{jj}=-\gamma_j$；
3. 在 $M = 2$、状态变量位于适中正值范围时，列标准化后设计矩阵条件数通常可控，不需要正交多项式（后者会将线性函数的信号分散至所有 $M$ 列，破坏组稀疏结构）。

### 3.2 逐节点设计矩阵

将式 (6) 代入稳态方程 (4)，对目标节点 $j$ 在 $N$ 个扰动条件下堆叠：

$$\Psi \theta_j + \mu_j \mathbf{1}_N = -u_{\cdot j} + \varepsilon_j \tag{7}$$

其中设计矩阵 $\Psi \in \mathbb{R}^{N \times pM}$ 定义为：

$$\Psi_{k,\, (i-1)M+m} = \psi_m\!\left(\tilde{x}_i^{(k)}\right) = \left(\tilde{x}_i^{(k)}\right)^m \tag{8}$$

全局参数向量 $\theta_j = [\theta_{j1}^\top, \ldots, \theta_{jp}^\top]^\top \in \mathbb{R}^{pM}$，响应向量 $u_{\cdot j} = (u_j^{(1)}, \ldots, u_j^{(N)})^\top$。

**实现原则**：对每个目标节点 $j$ 独立构建 $\Psi$，逐节点循环回归；不构造 $p$ 个节点合并的块对角大矩阵（理论上等价，计算上高效）。

---

## 4. 中心化与列标准化

### 4.1 列中心化（消除截距）

直接用式 (7) 估计 $\mu_j$ 与 $\theta_j$ 时，由于 $\psi_m(\cdot)$ 的列均值非零，会导致两者联合估计数值不稳定。标准处理为列均值中心化：

$$\bar{\psi}_{im} = \frac{1}{N} \sum_{k=1}^{N} \psi_m\!\left(\tilde{x}_i^{(k)}\right), \qquad \Psi_c = \Psi - \mathbf{1}_N \bar{\Psi}^\top \tag{9}$$

$$\bar{u}_j = \frac{1}{N} \sum_{k=1}^{N} u_j^{(k)}, \qquad u_{c,j} = u_{\cdot j} - \bar{u}_j \mathbf{1}_N \tag{10}$$

中心化后，式 (7) 化为无截距线性模型：

$$\Psi_c \theta_j = -u_{c,j} + \varepsilon_j \tag{11}$$

### 4.2 截距恢复

对式 (7) 两端取均值，得：

$$\hat{\mu}_j = -\bar{u}_j - \bar{\Psi} \hat{\theta}_j \tag{12}$$

其中 $\bar{\Psi} = [\bar{\psi}_{11}, \ldots, \bar{\psi}_{pM}] \in \mathbb{R}^{1 \times pM}$。在线性稳态真值下，$\hat{\mu}_j$ 收敛到目标节点的基线漂移；在 gLV 示例中对应 $r_j$。

### 4.3 列标准化（改善条件数）

单项式基 $x$ 与 $x^2$ 在有限状态范围内可能高度共线，导致中心化设计矩阵 $\Psi_c$ 列间尺度差异显著。进一步对每列除以样本标准差：

$$\hat{\sigma}_{im} = \operatorname{std}\!\left(\Psi_{c,\, :,\, (i-1)M+m}\right), \qquad \Psi_{cs} = \Psi_c \operatorname{diag}(\hat{\sigma})^{-1} \tag{13}$$

稀疏回归在标准化矩阵 $\Psi_{cs}$ 上进行。估计完毕后，系数反标准化恢复至原始尺度：

$$\hat{\theta}_j = \hat{\theta}_j^{(s)} \oslash \hat{\sigma} \tag{14}$$

其中 $\oslash$ 表示逐元素除法。**实现保护**：$\hat{\sigma}_{im} \leftarrow \max(\hat{\sigma}_{im},\, 10^{-10})$，防止常数列导致除零。

---

## 5. 双稀疏结构与回归

### 5.1 双稀疏假设

参数集 $\{\theta_{jim}\}$ 具有两层自然稀疏性：

**组间稀疏**（S1）：真实耦合网络稀疏，即对每个目标节点 $j$，有效输入源的数目远小于 $p$：

$$\left|\{i\neq j : \|\theta_{ji,\cdot}\|_2 > 0\}\right| \leq s \ll p \tag{S1}$$

**组内稀疏**（S2）：每条跨节点边的函数形式相对简单（如线性、低阶非线性或初级饱和），仅需少数基函数：

$$\|\theta_{ji,\cdot}\|_0 \leq s_0 \ll M, \qquad \forall (j, i) \in E \tag{S2}$$

联合约束 (S1)+(S2) 称为**双稀疏**（double sparsity），由 ADSIHT 算法严格支持。

### 5.2 ADSIHT（首选方法）

**算法**：Adaptive Double Sparse Iterative Hard Thresholding（Zhang et al. 2023）。对如下组合约束优化问题：

$$\min_{\theta_j}\ \frac{1}{2N}\|\Psi_{cs}\,\theta_j + u_{c,j}\|_2^2 \quad \text{s.t. (S1) and (S2)} \tag{15}$$

ADSIHT 通过迭代硬阈值步骤在组间和组内同时施加稀疏性，理论上达到 minimax 最优估计率。

**调用方式**：

```r
library(ADSIHT)
group   <- rep(seq_len(p), each = M)
fit     <- ADSIHT(Psi_cs, matrix(rhs_c), group, ic.type = "dsic")
best    <- which.min(fit$ic)
theta_j <- fit$beta[, best] / Psi_sd   # 反标准化
```

**关键参数**：
- `group = rep(1:p, each=M)`：将每个源变量的 $M$ 个基列归为同一组；
- `ic.type = "dsic"`：双稀疏信息准则，自动选择最优稀疏度组合 $(s, s_0)$，无需手动调参；
- `kappa = 0.9`（默认）：每个激活组内保留 $\lfloor \kappa M \rfloor$ 个系数（$M = 2$ 时保留 1 个）；
- ADSIHT 产生**精确的组零**，可用浮点安全阈值 $\tau = 10^{-10}$ 直接识别。

### 5.3 Group Lasso（次选方法）

作为对比基线，采用 `grpreg` 包的 Group Lasso（grLasso）：

$$\min_{\theta_j}\left\{\frac{1}{2N}\|\Psi_{cs}\,\theta_j + u_{c,j}\|_2^2 + \lambda \sum_{i=1}^p \|\theta_{ji,\cdot}\|_2\right\} \tag{16}$$

**调用方式**：

```r
library(grpreg)
cv      <- cv.grpreg(Psi_cs, rhs_c, group = group,
                     penalty = "grLasso", nfolds = 5)
theta_j <- coef(cv)[-1] / Psi_sd
```

Group Lasso 的惩罚为连续收缩，不产生精确的组零；边判定需使用相对阈值（见 §6.3）。

---

## 6. Jacobian 提取与网络重构

### 6.1 Jacobian 矩阵

由拟合的加性模型，在无扰动参考稳态 $x^{\mathrm{ref}}$（式 5）处计算 Jacobian：

$$J_{ji} = \left.\frac{\partial}{\partial x_i}\left[\mu_j + f_{jj}(x_j) + \sum_{i' \neq j} f_{ji'}(x_{i'})\right]\right|_{x = x^{\mathrm{ref}}} = \psi'(x_i^{\mathrm{ref}})^\top \hat{\theta}_{ji} \tag{17}$$

对单项式基，导数解析可得 $\psi'_m(x) = m x^{m-1}$，因此：

$$J_{ji} = \sum_{m=1}^{M} m \cdot (x_i^{\mathrm{ref}})^{m-1} \cdot \hat{\theta}_{jim} \tag{18}$$

### 6.2 线性效应与 gLV 参数提取

对一般复杂系统，$J_{ji}$ 即参考稳态附近变量 $i$ 对变量 $j$ 的局部线性效应。若使用 gLV 作为 benchmark，则可进一步解释为

$$\hat{A}_{ji} = J_{ji}\ (i \neq j), \qquad \hat{\gamma}_j = -J_{jj} \tag{19}$$

$J_{ji}>0$ 表示变量 $i$ 在参考稳态附近对变量 $j$ 有正向效应；$J_{ji}<0$ 表示负向效应。gLV 示例中，非自环的 $J_{ji}$ 对应 $\hat{A}_{ji}$。

### 6.3 边判定

采用群 L2 范数阈值：

$$\hat{E} = \left\{(j \leftarrow i) : \|\hat{\theta}_{ji,\cdot}\|_2 \geq \tau,\ i \neq j\right\} \tag{20}$$

| 方法 | 阈值 | 依据 |
|------|------|------|
| ADSIHT | $\tau = 10^{-10}$ | IHT 产生精确组零；浮点安全阈值直接识别非零组 |
| grLasso | $\tau_j = 0.01 \cdot \max_{i\neq j} \|\hat{\theta}_{ji,\cdot}\|_2$ | 连续收缩不产生精确零；按非自环候选边的行最大范数 1% 自适应截断 |

---

## 7. 评估指标

### 7.1 二元边分类

$$\text{Precision} = \frac{\text{TP}}{\text{TP} + \text{FP}}, \quad \text{Recall} = \frac{\text{TP}}{\text{TP} + \text{FN}}, \quad F_1 = \frac{2 \cdot \text{Pr} \cdot \text{Re}}{\text{Pr} + \text{Re}}$$

### 7.2 Matthews 相关系数（MCC）

MCC 是不平衡分类问题的推荐指标，取值范围 $[-1, 1]$（$+1$ 为完美分类）：

$$\text{MCC} = \frac{\text{TP} \cdot \text{TN} - \text{FP} \cdot \text{FN}}{\sqrt{(\text{TP}+\text{FP})(\text{TP}+\text{FN})(\text{TN}+\text{FP})(\text{TN}+\text{FN})}} \tag{21}$$

### 7.3 系数 L2 误差

$$\text{CoefL2} = \frac{1}{p} \sum_{j=1}^p \|\hat{\theta}_j - \theta_j^{\mathrm{true}}\|_2 \tag{22}$$

其中线性稳态真值为 $\theta_{ji}^{\mathrm{true}}=(C_{ji},0,\ldots,0)^\top$。在 gLV benchmark 中，

$$\theta_{ji}^{\mathrm{true}} =
\begin{cases}
(A_{ji}, 0, \ldots, 0)^\top, & i\neq j,\\
(-\gamma_j, 0, \ldots, 0)^\top, & i=j,
\end{cases} \tag{22a}$$

即跨节点耦合和自反馈都只在第一阶非零，但含义不同。

### 7.4 Jacobian RMSE

对一般线性 benchmark，非自环 Jacobian 误差为

$$\text{JacRMSE} = \sqrt{\frac{1}{p(p-1)} \sum_{j \neq i} (\hat{J}_{ji} - J_{ji}^{\mathrm{true}})^2} \tag{23}$$

在 gLV benchmark 中，$J_{ji}^{\mathrm{true}}=A_{ji}$（$i\neq j$），因此式 (23) 等价于比较 $\hat{A}_{ji}$ 与 $A_{ji}$。

---

## 8. 动力学验证与效应分解

### 8.1 重建 ODE

基于推断系数 $\hat{\theta}_j$ 和截距 $\hat{\mu}_j$，构建重建 ODE 系统：

$$\frac{d\hat{x}_j}{dt} = \hat{\mu}_j + \psi(x_j)^\top \hat{\theta}_{jj} + \sum_{i\neq j} \psi(x_i)^\top \hat{\theta}_{ji} + u_j \tag{24}$$

以 $x^{\mathrm{ref}}$ 为初值，在测试扰动条件下积分至 $t = 200$，与真实 ODE 轨迹对比。

### 8.2 效应分解

沿重建轨迹，分解单个源变量贡献随时间的变化。估计效应为

$$\hat{f}_{ji}(x_i(t)) = \sum_{m=1}^{M} \hat{\theta}_{jim} \cdot x_i(t)^m \tag{25}$$

若使用 gLV benchmark，则真值效应为

$$f_{ji}^{\mathrm{true}}(x_i(t)) =
\begin{cases}
A_{ji} \cdot x_i(t), & i\neq j,\\
-\gamma_j \cdot x_j(t), & i=j,
\end{cases} \tag{25a}$$

---

## 9. 网络可视化配色规范

| 类型 | 颜色 | R 色名 | 含义 |
|------|------|------|------|
| 正向边（$J_{ji} > 0,\ i\neq j$） | 红色 | `tomato3` | 源变量升高会增加目标变量的局部变化率 |
| 负向边（$J_{ji} < 0,\ i\neq j$） | 蓝色 | `steelblue3` | 源变量升高会降低目标变量的局部变化率 |
| FP 误报边 | 橙色 | `orange2`（虚线） | 推断多出的假阳性边 |
| FN 遗漏边 | 灰色 | `grey60`（点虚线叠加） | 推断漏掉的假阴性边 |

---

## 10. 完整推断流程

```
输入: {x̃^(k), u^(k)}, k=1..N
         ↓
[Step 1] 模拟扰动稳态
         X_obs (N×p, 含噪), U_obs (N×p)
         独立积分 u=0 → x_ref (无噪声参考稳态, eq.5)
         ↓
[Step 2] 构建单项式设计矩阵
         Ψ_{k,(i-1)M+m} = x̃_i^(k)^m  (eq.8)
         列中心化: Ψ_c = Ψ − Ψ̄       (eq.9)
         列标准化: Ψ_cs = Ψ_c / σ̂    (eq.13)
         ↓
[Step 3] for j = 1..p:
           rhs_c = −(u_j − ū_j)
           ADSIHT(Ψ_cs, rhs_c, group)  → θ̂_j^(s)  (DSIC 自动选模, eq.15)
           θ̂_j = θ̂_j^(s) / σ̂                       (反标准化, eq.14)
           μ̂_j = −ū_j − Ψ̄ θ̂_j                     (截距恢复, eq.12)
         ↓
[Step 4] Jacobian at x_ref
           J_{ji} = Σ_m m·(x_i^ref)^{m-1}·θ̂_{jim}  (eq.18)
           gLV benchmark: Â_{ji} = J_{ji} (i≠j),  γ̂_j = −J_{jj}  (eq.19)
         ↓
[Step 5] 边判定 (eq.20)
           ADSIHT: ‖θ̂_{ji}‖₂ ≥ 10^{-10}
           grLasso: ‖θ̂_{ji}‖₂ ≥ 0.01·max_{i≠j}‖θ̂_{ji}‖₂
         ↓
[Step 6] 动力学验证: x_ref → 积分重建 ODE (eq.24) vs 真实轨迹
         ↓
[Step 7] 效应分解: f̂_{ji}(x_i(t)); 若为 gLV benchmark,
                 i≠j: 对比 A_{ji}·x_i(t), i=j: 对比 −γ_j·x_j(t) (eq.25–25a)
         ↓
[Step 8] igraph 网络对比图 (TP/FP/FN 分色)
```

---

## 11. 可识别性与样本复杂度（理论保证）

> **说明**：本节速率与条件**照抄已发表结论**（Cai–Zhang–Zhou 2022；Zhang–Li–Liu–Yin
> 2024 即 ADSIHT；Lasso/group-Lasso 经典文献），仅做记号映射到 PSS-Net 设定，
> **未自行推导新定理**。标注 ⚠️ 处为需对照原文核实的细节。

### 11.1 双稀疏参数空间（照抄 ADSIHT 原文 Definition 1）

Zhang et al. (2024, ADSIHT) 对线性模型 $y = X\beta^* + \varepsilon$（$\varepsilon$ 为尺度
$\sigma^2$ 的 sub-Gaussian 噪声）定义**双稀疏**：系数 $\beta^*\in\mathbb{R}^p$ 称为
$(s, s_0)$-**sparse**，若

$$\|\beta^*\|_{0,2} := \sum_{j=1}^{m} \mathbb{I}(\beta^*_{G_j}\neq 0)\le s,
  \qquad
  \|\beta^*\|_{0} := \sum_{i=1}^{p}\mathbb{I}(\beta^*_i\neq 0)\le s\,s_0, \tag{D1}$$

其中 $\{G_j\}_{j=1}^m$ 为 $m$ 个不重叠分组。即**至多 $s$ 个非零组、总非零元至多 $s\,s_0$**。

**映射到 PSS-Net（逐节点回归 $\Psi_c\theta_j=-u_{c,j}+\varepsilon_j$，式 11）**：

| 通用记号 | PSS-Net 含义 |
|----------|--------------|
| 组数 $m$ | 节点/变量数 $p$（每个源变量 $i$ 的 $M$ 个基列为一组） |
| 组大小 | 基维数 $M$ |
| 样本数 $n$ | 扰动条件数 $N$ |
| 非零组数 $s$ | 目标节点 $j$ 的**有效入边数**（对应 §5.1 的 S1） |
| 组内稀疏 $s_0$ | 每条边的**有效基函数数**（对应 S2，低阶耦合 $s_0\ll M$） |

故 PSS-Net 的双稀疏结构 (S1)+(S2) **恰好**是 (D1)，可直接套用其理论。

### 11.2 双稀疏估计的 minimax 速率（照抄 Cai–Zhang–Zhou 2022, Thm 5）

设 $s_g$ 为非零组数、$s$ 为总非零元数、$d$ 为组数、$b$ 为最大组大小。该文给出估计误差的
**minimax 下界**（其 Theorem 5）

$$\inf_{\hat\beta}\sup_{\beta^*}\mathbb{E}\|\hat\beta-\beta^*\|_2^2
  \;\asymp\; \frac{\sigma^2}{n}\Big( s_g\log\frac{d}{s_g} + s\log\frac{e\,s_g\,b}{s}\Big). \tag{M1}$$

代入 PSS-Net 记号（$s_g\!\to\! s$ 入边数，$d\!\to\! p$，总非零 $s\!\to\! s\,s_0$，$b\!\to\! M$，$n\!\to\! N$）：

$$\boxed{\;\inf_{\hat\theta_j}\sup\ \mathbb{E}\|\hat\theta_j-\theta_j^*\|_2^2
  \;\asymp\; \frac{\sigma^2}{N}\Big( \underbrace{s\log\frac{p}{s}}_{\text{选对哪些边}}
  + \underbrace{s\,s_0\log\frac{eM}{s_0}}_{\text{选对哪些基并估值}}\Big).\;} \tag{M1'}$$

两项分别为"在 $p$ 个候选源中定位 $s$ 条边"与"在 $M$ 维基中定位 $s_0$ 个并估计"的代价。

**ADSIHT 的最优性**：Zhang et al. (2024) 证明其自适应 IHT 过程的上界与上述下界**相匹配**
（在适当条件下达到 minimax 最优），且**无需预知 $s, s_0, \sigma$**——这正是本项目选 ADSIHT
为首选求解器的理论依据。⚠️（上界的精确常数与对数因子以原文 Theorem 为准。）

### 11.3 设计条件（决定可识别性，连接最优扰动设计）

上述速率成立需对设计矩阵施加条件：

- **Cai–Zhang–Zhou (2022)** 假设 sub-Gaussian 设计：$X$ 各行 i.i.d. 中心化 sub-Gaussian、
  协方差 $\Sigma$ 特征值有界（含 Gaussian 设计），采用 approximate dual certificate 技术
  （**非** RIP/RE）。
- **IHT 类方法**（ADSIHT 及 best-subset 文献）通常要求设计满足**稀疏特征值 / RIP 型条件**
  （限制等距），即对所有 $(s,s_0)$-稀疏向量 $\theta$，
  $(1-\delta)\|\theta\|_2^2\le \tfrac1N\|\Psi_c\theta\|_2^2\le(1+\delta)\|\theta\|_2^2$。⚠️（具体常数见原文。）

**这是 PSS-Net 的关键衔接点**：设计矩阵 $\Psi_c$ 由**扰动稳态** $x^{*(k)}=x^*(u^{(k)})$ 经基展开
生成（式 8–9），故其 RIP/稀疏特征值**由扰动设计 $\{u^{(k)}\}$ 决定**。

> **可识别性命题（条件式陈述，非新证）**：若所选扰动 $\{u^{(k)}\}_{k=1}^N$ 使中心化设计
> $\Psi_c$ 满足 §11.3 的稀疏特征值/RIP 条件，则由 Zhang et al. (2024)，逐节点 ADSIHT 估计
> $\hat\theta_j$ 以高概率达到速率 (M1')；进而当**信号强度满足 beta-min 条件**
> $\min_{(j\leftarrow i)\in E}\|\theta^*_{ji,\cdot}\|_2 \gtrsim \sigma\sqrt{(\log p)/N}$
> 时，边集 $\hat E$ 与真值 $E$ 一致（精确支撑恢复）。

由此得**样本复杂度的量级**（照抄 Cai et al. 2022 Thm 1 的精确恢复条件，映射后）：

$$N \;\gtrsim\; s\log\frac{p}{s} \;+\; s\,s_0\log(e\,s\,M). \tag{SC}$$

即每个目标节点所需扰动条件数随**有效入边数 $s$**、**基复杂度 $s_0$** 线性增长，随节点数 $p$
仅**对数**增长——稀疏性使高维（$p$ 大）可控。**创新点 A（最优扰动设计）的作用**正是：在固定
预算 $N$ 下选 $\{u^{(k)}\}$ 以改善 $\Psi_c$ 的稀疏特征值/RIP，从而以更小的 $N$ 满足 (SC)。

### 11.4 基线方法速率（对照）

| 方法 | $\ell_2$ 估计误差（平方）量级 | 条件 / 出处 |
|------|------------------------------|-------------|
| **Lasso** | $\dfrac{\sigma^2\, k\log p}{N}$（$k=$ 总非零数）| 限制特征值 (RE)；Bickel–Ritov–Tsybakov 2009 |
| **Group Lasso** | $\dfrac{\sigma^2\,(sM + s\log p)}{N}$（$s$ 活跃组、组大小 $M$）| 组-RE；Lounici–Pontil–van de Geer–Tsybakov 2011 |
| **ADSIHT（双稀疏）** | $\dfrac{\sigma^2}{N}\big(s\log\frac{p}{s}+s\,s_0\log\frac{eM}{s_0}\big)$ | (M1')；Zhang et al. 2024 |

对比可见：当组内真稀疏（$s_0\ll M$，低阶耦合）时，**双稀疏速率的第二项
$s\,s_0\log(eM/s_0)$ 远小于 group Lasso 的 $sM$**——这定量解释了 §5 与模拟中 ADSIHT
（双稀疏）优于 group Lasso 的现象。

---

## 参考文献

1. Henderson, J. & Michailidis, G. (2014). Network reconstruction using nonparametric additive ODE models. *PLOS ONE*.
2. Wu, S., et al. (2014). Sparse additive ODEs for dynamic gene regulatory network recovery. *Journal of the American Statistical Association*.
3. Zhang, X., et al. (2023). Minimax optimal estimation in linear regression via adaptive double sparse iterative hard thresholding.
4. Barzel, B. & Barabási, A.-L. (2013). Universality in network dynamics. *Nature Physics*, 9, 673–681.
5. Meister, A., et al. (2013). Learning a nonlinear dynamical system model of gene regulation: A perturbed steady-state approach.
6. Zhang, Y., Li, Z., Liu, S. & Yin, J. (2024). A minimax optimal approach to high-dimensional double sparse linear regression. *Journal of Machine Learning Research*, 25. arXiv:2305.04182.（即 ADSIHT；§11 双稀疏定义 D1 与最优性）
7. Cai, T. T., Zhang, A. R. & Zhou, Y. (2022). Sparse group Lasso: Optimal sample complexity, convergence rate, and statistical inference. *IEEE Transactions on Information Theory*, 68(9). arXiv:1909.09851.（§11 双稀疏 minimax 速率 M1 / 样本复杂度 SC）
8. Bickel, P. J., Ritov, Y. & Tsybakov, A. B. (2009). Simultaneous analysis of Lasso and Dantzig selector. *The Annals of Statistics*, 37(4), 1705–1732.（§11.4 Lasso 速率与 RE 条件）
9. Lounici, K., Pontil, M., van de Geer, S. & Tsybakov, A. B. (2011). Oracle inequalities and optimal inference under group sparsity. *The Annals of Statistics*, 39(4), 2164–2204.（§11.4 group Lasso 速率）

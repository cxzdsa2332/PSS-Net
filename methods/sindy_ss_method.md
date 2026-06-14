# PSS-Net: Sparse Additive Nonparametric ODE Inference from Perturbed Steady States

**方法说明文档** · PSS-Net v3 · 对应实现：`script/sindy_ss_v3.R`

---

## 1. 问题设定

### 1.1 生态动力学模型

考虑由 $p$ 个物种构成的微生物群落，其群落动力学由如下耦合常微分方程（ODE）系统刻画：

$$\frac{dx_j(t)}{dt} = \mu_j + \sum_{i=1}^{p} f_{ji}\!\left(x_i(t)\right) + u_j(t), \qquad j = 1, \ldots, p \tag{1}$$

其中：
- $x_j(t) \geq 0$ 为物种 $j$ 在时刻 $t$ 的丰度；
- $\mu_j \in \mathbb{R}$ 为物种 $j$ 的内禀增长/死亡率（intrinsic rate）；
- $f_{ji} : \mathbb{R}_{\geq 0} \to \mathbb{R}$ 为**未知**的单变量互作函数，描述物种 $i$ 对物种 $j$ 的非线性效应；
- $u_j(t)$ 为施加于物种 $j$ 的外部扰动（如抗生素、底物浓度或人工添加/去除）。

本文采用**广义 Lotka-Volterra（GLV）模型**作为模拟真值：

$$\frac{dx_j}{dt} = r_j + \sum_{i \neq j} A_{ji} x_i - \gamma_j x_j + u_j \tag{2}$$

其中 $r_j$ 为内禀增长率，$A_{ji}$ 为物种 $i$ 对物种 $j$ 的线性互作强度，$\gamma_j > 0$ 为自我调节（死亡/竞争）率。式 (2) 是式 (1) 在 $f_{ji}(x_i) = A_{ji} x_i$（$i \neq j$）、$f_{jj}(x_j) = -\gamma_j x_j$ 下的特例。

### 1.2 互作网络定义

定义有向图 $G = (V, E)$，节点集 $V = \{1, \ldots, p\}$，有向边集

$$E = \{(j \leftarrow i) : f_{ji} \not\equiv 0\}$$

即物种 $i$ 对物种 $j$ 存在非零因果效应。**网络重构**的目标是：从观测数据中估计 $\hat{E}$ 使其尽可能接近真实 $E$。

### 1.3 可识别性约束

为消除截距与互作函数之间的不可识别性（$\mu_j$ 与 $f_{ji}$ 的常数项无法独立估计），施加如下约束：

$$f_{ji}(0) = 0, \qquad \forall j, i \tag{3}$$

**生态学含义**：物种 $i$ 丰度为零时对物种 $j$ 无互作贡献，符合生态直觉。在约束 (3) 下，$\mu_j$ 唯一可识别，且等于系统内禀增长率 $r_j$。

**实现**：选取满足 $\psi_m(0) = 0$（$\forall m$）的基函数族，即单项式基 $\psi_m(x) = x^m$（$m = 1, \ldots, M$）天然满足此约束。

---

## 2. 扰动稳态实验（PSS）

### 2.1 稳态方程

在稳态 $\dot{x}_j = 0$ 时，式 (1) 退化为代数约束：

$$\mu_j + \sum_{i=1}^{p} f_{ji}\!\left(x_i^*\right) + u_j = 0 \tag{4}$$

**核心优势**：无需对时序数据进行数值微分。观测到的稳态丰度 $x^*$ 直接进入线性回归，避免了导数估计引入的额外噪声。

### 2.2 扰动设计

对系统施加 $N$ 个独立扰动条件 $\{u^{(k)}\}_{k=1}^N$，每次等待系统达到新稳态 $x^{(k)*}$，记录观测对：

$$\mathcal{D} = \left\{\left(\tilde{x}^{(k)},\ u^{(k)}\right)\right\}_{k=1}^{N}$$

其中 $\tilde{x}^{(k)} = x^{(k)*} + \varepsilon^{(k)}$ 为含加性测量噪声的观测值（$\varepsilon^{(k)} \sim \mathcal{N}(0, \sigma^2 I)$）。

**实现细节**（`sindy_ss_v3.R` Step 1）：

- $N_{\text{cond}} = 300$ 个扰动条件；扰动分量 $u_j^{(k)} \sim \text{Uniform}(-0.4,\ 0.8)$，独立采样；
- 第一个条件设为野生型：$u^{(1)} = \mathbf{0}$；
- 每个条件以随机初值 $x_0 \sim \text{Uniform}(0.5, 2.0)^p$ 积分至 $t = 5000$（`rtol=1e-9`、`atol=1e-11`），取末值作为稳态估计；
- 过滤掉未收敛（含 `NaN`/`Inf` 或负值）的条件，保留有效条件 $N \leq N_{\text{cond}}$；
- 加入测量噪声 $\sigma = 0.03$，以 $x_{\min} = 10^{-6}$ 截断防止零值。

### 2.3 野生型稳态的无噪声估计

Jacobian 评估和动力学重建需要精确的野生型稳态 $x^{\mathrm{wt}}$。直接使用观测值 $\tilde{x}^{(1)}$ 会将 $\sigma = 0.03$ 的测量噪声引入后续所有非线性运算，产生系统性偏差。

**解决方案**：以 $u = \mathbf{0}$ 独立运行高精度 ODE 积分（$t_{\max} = 10^4$，`rtol=1e-12`，`atol=1e-14`），获得无噪声稳态：

$$x^{\mathrm{wt}} = \lim_{t \to \infty} x(t;\, u = \mathbf{0}) \tag{5}$$

$x^{\mathrm{wt}}$ 仅用于 Jacobian 计算（Step 4）和动力学对比（Step 6–7），不参与回归设计矩阵的构建。

---

## 3. 非参数基展开

### 3.1 函数近似

采用 $M$ 阶截断单项式基近似未知互作函数：

$$f_{ji}(a) \approx \sum_{m=1}^{M} \theta_{jim} \cdot \psi_m(a) = \psi(a)^\top \theta_{ji}, \qquad \psi_m(a) = a^m \tag{6}$$

其中 $\theta_{ji} = (\theta_{ji1}, \ldots, \theta_{jiM})^\top \in \mathbb{R}^M$ 为待估系数向量。

**为何选用单项式基**：

1. 自然满足约束 $\psi_m(0) = 0$（式 3），无需额外处理；
2. 对线性 GLV 真值（$f_{ji}(x) = A_{ji} x$），仅第一阶系数非零（$\theta_{ji1} = A_{ji}$，$\theta_{ji2} = \cdots = \theta_{jiM} = 0$），真实稀疏结构与基展开完全兼容；
3. 在 $M = 2$、$x \in [0.5, 3]$ 的典型生物丰度范围内，列标准化后设计矩阵条件数可控（$\leq 10^3$），不需要正交多项式（后者会将线性函数的信号分散至所有 $M$ 列，破坏组稀疏结构）。

### 3.2 逐节点设计矩阵

将式 (6) 代入稳态方程 (4)，对目标物种 $j$ 在 $N$ 个扰动条件下堆叠：

$$\Psi \theta_j + \mu_j \mathbf{1}_N = -u_{\cdot j} + \varepsilon_j \tag{7}$$

其中设计矩阵 $\Psi \in \mathbb{R}^{N \times pM}$ 定义为：

$$\Psi_{k,\, (i-1)M+m} = \psi_m\!\left(\tilde{x}_i^{(k)}\right) = \left(\tilde{x}_i^{(k)}\right)^m \tag{8}$$

全局参数向量 $\theta_j = [\theta_{j1}^\top, \ldots, \theta_{jp}^\top]^\top \in \mathbb{R}^{pM}$，响应向量 $u_{\cdot j} = (u_j^{(1)}, \ldots, u_j^{(N)})^\top$。

**实现原则**：对每个目标物种 $j$ 独立构建 $\Psi$，逐节点循环回归；不构造 $p$ 个节点合并的块对角大矩阵（理论上等价，计算上高效）。

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

其中 $\bar{\Psi} = [\bar{\psi}_{11}, \ldots, \bar{\psi}_{pM}] \in \mathbb{R}^{1 \times pM}$。在线性 GLV 真值下，$\hat{\mu}_j \to r_j$（内禀增长率），可用于验证推断质量。

### 4.3 列标准化（改善条件数）

单项式基 $x$ 与 $x^2$ 在典型丰度范围内高度共线，导致中心化设计矩阵 $\Psi_c$ 列间尺度差异显著。进一步对每列除以样本标准差：

$$\hat{\sigma}_{im} = \operatorname{std}\!\left(\Psi_{c,\, :,\, (i-1)M+m}\right), \qquad \Psi_{cs} = \Psi_c \operatorname{diag}(\hat{\sigma})^{-1} \tag{13}$$

稀疏回归在标准化矩阵 $\Psi_{cs}$ 上进行。估计完毕后，系数反标准化恢复至原始尺度：

$$\hat{\theta}_j = \hat{\theta}_j^{(s)} \oslash \hat{\sigma} \tag{14}$$

其中 $\oslash$ 表示逐元素除法。**实现保护**：$\hat{\sigma}_{im} \leftarrow \max(\hat{\sigma}_{im},\, 10^{-10})$，防止常数列导致除零。

---

## 5. 双稀疏结构与回归

### 5.1 双稀疏假设

参数集 $\{\theta_{jim}\}$ 具有两层自然稀疏性：

**组间稀疏**（S1）：真实互作网络稀疏，即对每个目标物种 $j$，有效调控源的数目远小于 $p$：

$$\left|\{i : \|\theta_{ji,\cdot}\|_2 > 0\}\right| \leq s \ll p \tag{S1}$$

**组内稀疏**（S2）：每条互作边的函数形式简单（线性竞争、初级饱和），仅需少数基函数：

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
- `group = rep(1:p, each=M)`：将每个物种的 $M$ 个基列归为同一组；
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

由拟合的加性模型，在野生型稳态 $x^{\mathrm{wt}}$（式 5）处计算 Jacobian：

$$J_{ji} = \left.\frac{\partial}{\partial x_i}\left[\mu_j + \sum_{i'} f_{ji'}(x_{i'})\right]\right|_{x = x^{\mathrm{wt}}} = \psi'(x_i^{\mathrm{wt}})^\top \hat{\theta}_{ji} \tag{17}$$

对单项式基，导数解析可得 $\psi'_m(x) = m x^{m-1}$，因此：

$$J_{ji} = \sum_{m=1}^{M} m \cdot (x_i^{\mathrm{wt}})^{m-1} \cdot \hat{\theta}_{jim} \tag{18}$$

### 6.2 参数提取

$$\hat{A}_{ji} = J_{ji}\ (i \neq j), \qquad \hat{\gamma}_j = -J_{jj} \tag{19}$$

$\hat{A}_{ji} > 0$ 表示物种 $i$ 对物种 $j$ 有促进（promotion）效应；$\hat{A}_{ji} < 0$ 表示抑制（inhibition）效应。

### 6.3 边判定

采用群 L2 范数阈值：

$$\hat{E} = \left\{(j \leftarrow i) : \|\hat{\theta}_{ji,\cdot}\|_2 \geq \tau,\ i \neq j\right\} \tag{20}$$

| 方法 | 阈值 | 依据 |
|------|------|------|
| ADSIHT | $\tau = 10^{-10}$ | IHT 产生精确组零；浮点安全阈值直接识别非零组 |
| grLasso | $\tau_j = 0.01 \cdot \max_i \|\hat{\theta}_{ji,\cdot}\|_2$ | 连续收缩不产生精确零；按行最大范数的 1% 自适应截断 |

---

## 7. 评估指标

### 7.1 二元边分类

$$\text{Precision} = \frac{\text{TP}}{\text{TP} + \text{FP}}, \quad \text{Recall} = \frac{\text{TP}}{\text{TP} + \text{FN}}, \quad F_1 = \frac{2 \cdot \text{Pr} \cdot \text{Re}}{\text{Pr} + \text{Re}}$$

### 7.2 Matthews 相关系数（MCC）

MCC 是不平衡分类问题的推荐指标，取值范围 $[-1, 1]$（$+1$ 为完美分类）：

$$\text{MCC} = \frac{\text{TP} \cdot \text{TN} - \text{FP} \cdot \text{FN}}{\sqrt{(\text{TP}+\text{FP})(\text{TP}+\text{FN})(\text{TN}+\text{FP})(\text{TN}+\text{FN})}} \tag{21}$$

### 7.3 系数 L2 误差

$$\text{CoefL2} = \frac{1}{p} \sum_{j=1}^p \|\hat{\theta}_j - \theta_j^{\mathrm{true}}\|_2 \tag{22}$$

其中 $\theta_{ji}^{\mathrm{true}} = (A_{ji}, 0, \ldots, 0)^\top \in \mathbb{R}^M$（线性 GLV 下仅第一阶非零）。

### 7.4 Jacobian RMSE

$$\text{JacRMSE} = \sqrt{\frac{1}{p(p-1)} \sum_{j \neq i} (\hat{A}_{ji} - A_{ji})^2} \tag{23}$$

---

## 8. 动力学验证与效应分解

### 8.1 重建 ODE

基于推断系数 $\hat{\theta}_j$ 和截距 $\hat{\mu}_j$，构建重建 ODE 系统：

$$\frac{d\hat{x}_j}{dt} = \hat{\mu}_j + \sum_{i=1}^{p} \psi(x_i)^\top \hat{\theta}_{ji} + u_j \tag{24}$$

以 $x^{\mathrm{wt}}$ 为初值，在测试扰动条件下积分至 $t = 200$，与真实 ODE 轨迹对比。

### 8.2 效应分解

沿重建轨迹，分解单物种贡献随时间的变化：

$$f_{ji}^{\mathrm{true}}(x_i(t)) = A_{ji} \cdot x_i(t), \qquad \hat{f}_{ji}(x_i(t)) = \sum_{m=1}^{M} \hat{\theta}_{jim} \cdot x_i(t)^m \tag{25}$$

---

## 9. 网络可视化配色规范

| 类型 | 颜色 | R 色名 | 含义 |
|------|------|------|------|
| 促进边（$A_{ji} > 0$） | 红色 | `tomato3` | 正互作，互利/共生 |
| 抑制边（$A_{ji} < 0$） | 蓝色 | `steelblue3` | 负互作，竞争/寄生 |
| FP 误报边 | 橙色 | `orange2`（虚线） | 推断多出的假阳性边 |
| FN 遗漏边 | 灰色 | `grey60`（点虚线叠加） | 推断漏掉的假阴性边 |

---

## 10. 完整推断流程

```
输入: {x̃^(k), u^(k)}, k=1..N
         ↓
[Step 1] 模拟扰动稳态
         X_obs (N×p, 含噪), U_obs (N×p)
         独立积分 u=0 → x_wt (无噪声野生型稳态, eq.5)
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
[Step 4] Jacobian at x_wt
           J_{ji} = Σ_m m·(x_i^wt)^{m-1}·θ̂_{jim}  (eq.18)
           Â_{ji} = J_{ji} (i≠j),  γ̂_j = −J_{jj}   (eq.19)
         ↓
[Step 5] 边判定 (eq.20)
           ADSIHT: ‖θ̂_{ji}‖₂ ≥ 10^{-10}
           grLasso: ‖θ̂_{ji}‖₂ ≥ 0.01·max_i‖θ̂_{ji}‖₂
         ↓
[Step 6] 动力学验证: x_wt → 积分重建 ODE (eq.24) vs 真实轨迹
         ↓
[Step 7] 效应分解: f̂_{ji}(x_i(t)) vs A_{ji}·x_i(t) (eq.25)
         ↓
[Step 8] igraph 网络对比图 (TP/FP/FN 分色)
```

---

## 参考文献

1. Henderson, J. & Michailidis, G. (2014). Network reconstruction using nonparametric additive ODE models. *PLOS ONE*.
2. Wu, S., et al. (2014). Sparse additive ODEs for dynamic gene regulatory network recovery. *Journal of the American Statistical Association*.
3. Zhang, X., et al. (2023). Minimax optimal estimation in linear regression via adaptive double sparse iterative hard thresholding.
4. Barzel, B. & Barabási, A.-L. (2013). Universality in network dynamics. *Nature Physics*, 9, 673–681.
5. Meister, A., et al. (2013). Learning a nonlinear dynamical system model of gene regulation: A perturbed steady-state approach.

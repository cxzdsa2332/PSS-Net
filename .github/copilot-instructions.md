# PSS-Net

## 项目目标

对稳态与时序 ODE 数据进行稀疏非线性**微生物互作网络**重构，核心框架为**稀疏加性非参 ODE 模型**（Sparse Additive Nonparametric ODE），结合 B 样条基展开与双稀疏推断（ADSIHT）。

## Agent 指令同步规则

- `CLAUDE.md` 与 `.github/copilot-instructions.md` 是同一套项目级 agent 指令。
- 修改任一文件时，必须同步更新另一文件，保持内容一致。
- 完成修改后运行 `cmp -s CLAUDE.md .github/copilot-instructions.md` 检查两者是否完全一致。

## 目录结构

```
PSS-Net/
├── data/             # 真实数据（实验观测、稳态丰度矩阵等）
├── sim_script/       # 模拟与推断脚本（pss_net_*.R, sindy_ss_*.R）
├── analysis_script/  # 分析脚本（ODE推断、网络重构、可视化）
├── methods/          # 方法说明文档
├── results/          # 输出结果（mcc_comparison.csv 等）
├── ref/              # 参考文献（PDF、MD）
├── CLAUDE.md
└── .github/copilot-instructions.md
```

## 技术栈

- **主语言**：R（分析、建模）和 Python（辅助计算、优化）
- **绘图**：一律使用 `ggplot2`（R）；Python 端如需出图，输出数据后交 R 处理
- **网络可视化配色**（igraph / ggplot2 网络图统一规范）：
  - 促进边（$A_{ji} > 0$，正效应）：**红色** `tomato3`
  - 抑制边（$A_{ji} < 0$，负效应）：**蓝色** `steelblue3`
  - FP 误报边：**橙色** `orange2`（虚线）
  - FN 遗漏边：**灰色** `grey60`（点虚线叠加）
- **ODE 求解**：R 使用 `deSolve`；Python 使用 `scipy.integrate`
- **基展开**：R 使用 `splines::bs()`（B 样条，首选）或 Legendre 多项式；约束 $f_{ji}(0)=0$ → 基函数无截距（`intercept=FALSE`）
- **建模范式**：稀疏加性非参 ODE + 逐节点回归 + 双稀疏推断
  - 对每个目标节点 $j$ **独立**构建设计矩阵 $\Psi_j$（规避块对角大矩阵），逐节点循环推断
  - 对每个源变量 $x_i$ 构建 $M$ 个基函数 $\psi(x_i) \in \mathbb{R}^M$ 作为一个 group，天然支持 ADSIHT 的组稀疏结构（组大小 $M \geq 2$ 保证 `floor(kappa×M)≥1`）
  - **稳态路径**（PSS 数据）：$\dot{x}_j=0$ 给出 $\Psi_j \theta_j = -u_j$，无需数值微分
  - **积分路径**（时序数据）：对 ODE 两端积分消除导数噪声，转化为线性方程
  - **弱形式**（高噪声）：乘以测试函数后积分，进一步压制噪声
- **稀疏回归**（优先级依次降低）：
  1. **ADSIHT**（CRAN 包，`library(ADSIHT)`）— 首选；`ADSIHT(x, y, group, ic.type="dsic")`，`group=rep(1:p, each=M)`，DSIC 自动选模型；双稀疏：组间（哪些边存在）+ 组内（哪些基函数显著）
  2. **Sparse Group Lasso** — 次选，R 使用 `sparsegl` / `gglasso`
  3. **Lasso / Group Lasso** — 基线对比，R 使用 `glmnet`

## 编码规范

- 模拟与推断脚本放入 `sim_script/`，分析脚本放入 `analysis_script/`
- 探索性、一次性或历史草稿脚本放入 `sim_script/manual/`；该目录下脚本不必遵循“分析代码与绘图代码分离”规则，但不得作为正式结果生成入口
- 真实数据放入 `data/`，参考文献放入 `ref/`，输出结果放入 `results/`
- R 脚本使用 `snake_case` 命名
- 每个脚本顶部注明用途、输入输出
- 图形输出保存为 PDF 或 PNG，命名与脚本对应

### 分析代码与绘图代码分离（强制）

为保证透明与可复现，**计算与绘图必须分两步、两脚本**：

1. **正式分析/模拟脚本**（`sim_script/` 根目录、`analysis_script/`）：只负责计算，把**所有数值结果写入
   `results/sim_results/`**（CSV 等中间数据文件，命名与脚本对应，如 `design_nl.csv`）。
   不在分析脚本里直接出图。
2. **绘图脚本**（`analysis_script/` 下独立的 `plot_*.R`）：只负责读取
   `results/sim_results/` 中的数据并出图，输出到 `results/figure/`。**绘图脚本不重新跑模拟**。
3. **探索性脚本例外**：`sim_script/manual/` 中的脚本可混合计算、绘图和交互检查；这些脚本仅用于人工探索，不作为论文正式结果的复现入口。

理由：模拟可能耗时且随机，固化中间数据后，改图样式无需重算；他人只需 `results/sim_results/`
即可复现全部图表。`results/table/` 留作论文用的格式化表格。

> 若某个探索性脚本要升级为正式结果生成脚本，必须移出 `sim_script/manual/`，并拆分为计算脚本与 `analysis_script/plot_*.R` 或表格汇总脚本。

## 核心方法

### 1. ODE 模型

考虑 $p$ 个物种组成的微生物群落，其动力学由耦合 ODE 系统描述：

$$\frac{dx_j(t)}{dt} = \mu_j + \sum_{i=1}^{p} f_{ji}\!\left(x_i(t)\right), \quad j = 1,\ldots,p \tag{1}$$

其中 $x_j(t) \geq 0$ 为物种 $j$ 的丰度，$\mu_j$ 为内禀增长/死亡率，$f_{ji}(\cdot)$ 为物种 $i$ 对物种 $j$ 的**未知单变量互作函数**。

**网络定义**：有向图 $G = (V,E)$，$V = \{1,\ldots,p\}$；边 $(j \leftarrow i) \in E$ 当且仅当 $f_{ji} \not\equiv 0$（物种 $i$ 对 $j$ 有非零互作效应）。

**稀疏性假设**：真实互作网络稀疏，即对每个目标物种 $j$，有效调控源 $\{i : f_{ji} \not\equiv 0\}$ 数目远小于 $p$。

**可识别性约束**（Meister et al. 2013）：令

$$f_{ji}(0) = 0 \tag{2}$$

物种 $i$ 丰度为零时对 $j$ 无互作效应，生态上合理。该约束消除了 $\mu_j$ 与各 $f_{ji}$ 截距之间的不可识别性：$\mu_j = \mu_{j0} + \sum_i \mu_{ji}$ 中各 $\mu_{ji}$ 无法独立估计，约束 (2) 使得 $\mu_j$ 唯一且可估计。

---

### 2. 扰动稳态（PSS）实验设计

#### 2.1 稳态方程

在稳态 $\dot{x}_j = 0$ 时，(1) 式给出代数约束：

$$\mu_j + \sum_{i=1}^{p} f_{ji}\!\left(x_i^*\right) = 0 \tag{3}$$

#### 2.2 扰动设计

对系统施加外部扰动 $u_j^{(k)}$（如物种 $j$ 的添加/去除、抑制剂浓度），得到扰动动力学：

$$\frac{dx_j}{dt} = \mu_j + \sum_{i=1}^{p} f_{ji}(x_i) + u_j \tag{4}$$

在第 $k$ 个扰动条件下等待系统达到新稳态 $x^{(k)*}$，则：

$$\mu_j + \sum_{i=1}^{p} f_{ji}\!\left(x_i^{(k)*}\right) + u_j^{(k)} = 0 \tag{5}$$

收集 $N$ 个扰动条件（含野生型 $u^{(1)}=0$），得到观测数据矩阵：

$$X = \{x^{(k)}\}_{k=1}^N \in \mathbb{R}^{N \times p}, \quad U = \{u^{(k)}\}_{k=1}^N \in \mathbb{R}^{N \times p}$$

由 (5) 式，对每个目标物种 $j$，回归响应为 $-u_j^{(k)}$，无需数值微分。

---

### 3. 非参数基展开

#### 3.1 函数近似

用 $M$ 维截断基近似未知互作函数（Henderson & Wu 2014, Wu et al. 2014）：

$$f_{ji}(a) = \psi(a)^\top \theta_{ji} + \delta_{ji}(a), \quad \theta_{ji} \in \mathbb{R}^M \tag{6}$$

其中 $\psi(a) = (\psi_1(a), \ldots, \psi_M(a))^\top$ 为 $M$ 维基函数向量，$\delta_{ji}$ 为截断残差。

约束 (2) 要求 $\psi(0)^\top \theta_{ji} = 0$；最简单的满足方式是令 $\psi_m(0) = 0\ \forall m$，即**无截距基**：

$$\text{B 样条（首选）：} \psi_m = \text{第 } m \text{ 个 B 样条基},\ \mathtt{bs(x, df=M, intercept=FALSE)}$$
$$\text{Legendre 多项式（备选）：} \psi_m = P_m(x),\ m=1,\ldots,M,\ P_m(0)=0\ \text{（奇次项）}$$

#### 3.2 设计矩阵构造

将 (6) 代入 (5)，对第 $k$ 个扰动条件、目标物种 $j$：

$$\mu_j + \sum_{i=1}^{p} \psi\!\left(x_i^{(k)}\right)^\top \theta_{ji} = -u_j^{(k)} + \varepsilon_j^{(k)} \tag{7}$$

对 $N$ 个条件堆叠，引入**逐节点设计矩阵**：

$$\Psi = \bigl[\psi(x_1^{(1)}),\ldots,\psi(x_p^{(1)});\ldots;\psi(x_1^{(N)}),\ldots,\psi(x_p^{(N)})\bigr]$$

具体地，令

$$\Psi \in \mathbb{R}^{N \times pM}, \quad \Psi_{k,\,(i-1)M+m} = \psi_m\!\left(x_i^{(k)}\right)$$

参数向量 $\theta_j = [\theta_{j1}^\top, \ldots, \theta_{jp}^\top]^\top \in \mathbb{R}^{pM}$，则 (7) 化为线性模型：

$$\boxed{\Psi\, \theta_j + \mu_j \mathbf{1} = -u_{\cdot j} + \varepsilon_j} \tag{8}$$

**实现原则**：对每个目标物种 $j$ **独立**构建 $\Psi$，**不构造** $p$ 个节点合并的块对角大矩阵 $\mathbf{X} = \mathrm{diag}(\Psi,\ldots,\Psi) \in \mathbb{R}^{Np \times p^2M}$（该形式仅用于理论分析，计算中逐节点循环等价且高效）。

---

### 4. 中心化消除截距

$\Psi$ 的各列在约束 (2) 下不含常数列，但不同条件间的 $\psi_m(x_i^{(k)})$ 均值非零，导致设计矩阵满秩时 $\mu_j$ 与 $\theta_j$ 联合估计数值不稳定。标准处理：

**列中心化**：令 $\bar{\psi}_{im} = \frac{1}{N}\sum_k \psi_m(x_i^{(k)})$，定义

$$\Psi_c = \Psi - \bar{\Psi}, \qquad u_{c,j} = u_{\cdot j} - \bar{u}_j \mathbf{1}$$

(8) 式中心化后变为无截距形式：

$$\Psi_c\, \theta_j = -u_{c,j} + \varepsilon_j \tag{9}$$

**截距恢复**：由 (8) 式取均值：

$$\mu_j = -\bar{u}_j - \bar{\Psi}\, \theta_j \tag{10}$$

其中 $\bar{\Psi} = [\bar{\psi}_{1,1},\ldots,\bar{\psi}_{p,M}] \in \mathbb{R}^{1 \times pM}$。

---

### 5. 双稀疏回归

#### 5.1 双稀疏结构

参数集 $\{\theta_{jim}\}$ 具有两层稀疏性（对应 ADSIHT 的理论框架）：

**组间稀疏**（哪些物种对存在互作）：

$$\sum_{i=1}^{p} \mathbb{I}\!\left(\|\theta_{ji,\cdot}\|_2 \neq 0\right) \leq s \tag{S1}$$

边 $(j \leftarrow i)$ 存在 $\Leftrightarrow$ $\|\theta_{ji,\cdot}\|_2 \neq 0$（物种 $i$ 的全部基系数非零）

**组内稀疏**（每条边的互作函数形式简单）：

$$\sum_{m=1}^{M} \mathbb{I}(\theta_{jim} \neq 0) \leq s_0, \quad \forall (j,i) \tag{S2}$$

低阶生态互作（线性竞争、初级饱和）只需少数几个基函数。(S1)+(S2) 合称 **double sparsity**。

| 稀疏层次 | 约束 | 生态含义 |
|----------|------|----------|
| 组间（S1） | $\|\cdot\|_{2,0} \leq s$ | 每物种有效互作伙伴稀少 |
| 组内（S2） | $\|\cdot\|_0 \leq s_0$ per group | 互作函数多为低阶（线性竞争/初级饱和） |

#### 5.2 目标函数

对 (9) 求解：

$$\min_{\theta_j}\ \frac{1}{2N}\|\Psi_c\,\theta_j + u_{c,j}\|_2^2 \quad \text{s.t.}\ \text{(S1) and (S2)} \tag{10}$$

**Sparse Group Lasso 松弛**（次选，凸）：

$$\min_{\theta_j}\left\{\frac{1}{2N}\|\Psi_c\,\theta_j + u_{c,j}\|_2^2 + \lambda_1 \sum_{i=1}^p \|\theta_{ji,\cdot}\|_2 + \lambda_2 \|\theta_j\|_1\right\}$$

**ADSIHT**（首选，非凸但理论最优）：对 (10) 的组合约束问题使用迭代硬阈值算法。

#### 5.3 ADSIHT 调用方式

```r
library(ADSIHT)
# Psi_c: N × pM 中心化设计矩阵；rhs = -u_{c,j}
# group: rep(1:p, each=M)，物种 i 的 M 个基列归为 group i
fit  <- ADSIHT(Psi_c, matrix(rhs), group, ic.type = "dsic")
best <- which.min(fit$ic)
theta_j_hat <- fit$beta[, best]   # pM 维稀疏系数向量
```

关键参数：`ic.type="dsic"` 使用双稀疏信息准则自动选模型；`kappa=0.9`（默认）控制组内保留比例，组大小 $M=4$ 时每激活组保留 $\lfloor 0.9\times 4\rfloor = 3$ 个基。

---

### 6. 预估计（去噪）⚠️ 慎用

含噪观测 $\tilde{x}_j^{(k)} = x_j^{(k)*} + \varepsilon_j^{(k)}$ 直接代入基展开会将噪声放大。

**⚠️ 注意（实验验证）**：在多物种互作系统中，以下 B 样条平滑方案经 10 次重复实验证实**显著降低推断质量**（ADSIHT MCC 从 0.91 降至 0.25），**不建议在 PSS 推断中使用**。原因：真实稳态 $x_j^*$ 由整个 ODE 系统共同决定，并非 $u_j$ 的简单光滑函数；逐物种回归 $\hat{x}_j \sim \text{bs}(u_j)$ 会压缩跨物种协变结构，破坏设计矩阵的辨识性。

**（不推荐）B 样条平滑**（原理上适用于单物种或独立系统）：

$$\hat{x}_j = \text{fitted}\!\left(\text{lm}(\tilde{x}_j \sim \mathtt{bs}(u_j,\, \mathrm{df}=6))\right), \quad j=1,\ldots,p$$

**推荐替代**：直接使用含噪观测 $\tilde{X}$ 构建设计矩阵；在 $N \geq 100$、$\sigma \leq 0.05$ 的典型 PSS 条件下，ADSIHT 的 DSIC 准则对该噪声水平具有足够鲁棒性。

---

### 7. Jacobian 提取与网络重构

#### 7.1 Jacobian 矩阵

由拟合的加性模型，在野生型稳态 $x^\mathrm{wt}$（$u=0$ 条件的稳态值）处计算 Jacobian：

$$J_{ji} = \frac{\partial}{\partial x_i}\left[\mu_j + \sum_{i'} f_{ji'}(x_{i'})\right]_{x=x^\mathrm{wt}} = \frac{d f_{ji}}{d a}\bigg|_{a=x_i^\mathrm{wt}} = \psi'(x_i^\mathrm{wt})^\top \hat{\theta}_{ji} \tag{11}$$

其中 $\psi'(x) = \frac{d\psi}{dx}$ 用中心差分逼近：

$$\psi'(x) \approx \frac{\psi(x+\varepsilon) - \psi(x-\varepsilon)}{2\varepsilon}, \quad \varepsilon = 10^{-5}$$

#### 7.2 参数提取

$$\hat{A}_{ji} = J_{ji}\ (i \neq j) \qquad \text{（互作强度矩阵）}$$
$$\hat{\gamma}_j = -J_{jj} \qquad \text{（物种 } j \text{ 死亡/降解率）}$$

#### 7.3 网络判定

设检测阈值 $\tau > 0$：

$$\hat{E} = \left\{(j \leftarrow i) : \|\hat{\theta}_{ji,\cdot}\|_2 \geq \tau\right\}$$

正值 $\hat{A}_{ji} > 0$ 表示促进（互利/共生），负值 $\hat{A}_{ji} < 0$ 表示抑制（竞争/寄生）。

---

### 8. 完整流程

```
数据: {x̃^(k), u^(k)}, k=1..N
  ↓
[Step 1] B 样条去噪: x̂_j = fitted(lm(x̃_j ~ bs(u_j, df=6)))
  ↓
[Step 2] for j = 1..p:
    构造 Ψ ∈ R^{N×pM}: Ψ_{k,(i-1)M+m} = ψ_m(x̂_i^(k))
    中心化: Ψ_c = Ψ - Ψ̄,  rhs = -(u_j - ū_j)
    ADSIHT: θ̂_j = ADSIHT(Ψ_c, rhs, group=rep(1:p, each=M))$beta[best]
    截距:   μ̂_j = -ū_j - Ψ̄ θ̂_j
  ↓
[Step 3] Jacobian: J_ji = ψ'(x_i^wt)ᵀ θ̂_ji  (数值中心差分)
  ↓
[Step 4] Â = J (off-diagonal),  γ̂_j = -J_jj
  ↓
[Step 5] 网络边: Ê = {(j←i) : ‖θ̂_ji‖₂ ≥ τ}
```

## 参考文献

- Henderson, J. & Michailidis, G. (2014). *Network reconstruction using nonparametric additive ODE models.*
- Wu, S. et al. (2014). *Sparse additive ODEs for dynamic gene regulatory network recovery.*
- Zhang, X. et al. *Minimax optimal estimation via ADSIHT* (zhang2023minimax).
- Barzel, B. & Barabási, A.-L. (2013). *Universality in network dynamics.*
- `ref/Learning a nonlinear dynamical system model of gene regulation- A perturbed steady-state approach.pdf`
- `ref/ODE_solve.md`

---

## 9. 评估指标与基准

网络重构按**边的二分类**评估（边存在/不存在），统一报告以下指标，结果写入 `results/`：

- **混淆矩阵**：TP / FP / FN / TN（基于 $\|\hat\theta_{ji}\|_2 \ge \tau$ 的边判定）
- **Precision / Recall / F1**
- **MCC**（Matthews 相关系数，主指标；类别不平衡下比 F1 更稳健）
- **系数误差**：CoefL2（$\|\hat\theta - \theta^*\|_2$）
- **Jacobian 误差**：JacRMSE（$\hat J$ 与真值 $J^*$ 的 RMSE）

每个实验设置须**多 seed 重复**（当前 10 次）后报告均值±标准差。基准对照至少包含 ADSIHT（首选）vs Group Lasso（`grpreg::grLasso`，次选），见 `sim_script/pss_net_compare.R`。

## 10. 实验结论（持续更新）

- **预平滑有害**：B 样条预平滑在多物种系统中显著降低 MCC（0.91→0.25），PSS 推断直接使用含噪观测（详见 §6）。
- **ADSIHT 优于 Group Lasso**：在线性 GLV（v3）无平滑方案下，ADSIHT 平均 MCC ≈ 0.91，grLasso ≈ 0.55；FP 数显著更低（双稀疏抑制误报）。
- **非线性互作**（v1，含 $x^2$ 项）仍是难点：当前单项式 M=2 基下 MCC 偏低，待改进基函数与扰动设计。

---

## 论文撰写规范（manuscript/）

撰写或修改 `manuscript/` 下文稿时，须遵守学术写作的常见排版/术语约定。
**该规范已配置为 skill `manuscript-style`**——动笔前先调用它核对，主要包括：

- **基因/物种学名斜体**：基因名、物种拉丁学名（如 *Bacteroides fragilis*）用斜体；
  蛋白名一般正体。
- **缩写首次出现给全称**：每个缩写在正文首次出现时写"全称（缩写）"，此后用缩写
  （如"generalized Lotka–Volterra (gLV)"）。
- **数学/记号一致**：与 `methods/` 文档及 `preamble.tex` 宏保持同一套记号。
- **数字与单位**、**图表标题自洽**、**引用 `\citep`/`\citet` 用法** 等，详见 skill。

> 具体清单与检查方式见 `.claude/skills/manuscript-style/SKILL.md`。

---

## TODO List

### A. 细化模拟（sim_script/）
- [ ] **基函数对比**：单项式 vs B 样条（`bs`, `intercept=FALSE`）vs Legendre，在相同非线性互作下比较 MCC/JacRMSE。
- [ ] **非线性互作恢复**：改进 v1 场景（饱和型 Monod、Hill 函数等生态常见形式），评估各基函数的逼近能力。
- [ ] **维度/样本量扫描**：系统改变 $(p, N, s, \sigma)$，绘制 MCC 关于 $N/p$、噪声 $\sigma$ 的相图（ggplot2）。
- [ ] **扰动设计敏感性**：单物种 vs 多物种联合扰动、连续 vs 离散扰动幅度对可识别性的影响。
- [ ] **ADSIHT 超参**：`kappa`、组大小 $M$、`ic.type` 的稳健性测试。
- [ ] **稳健性**：稳态未完全收敛、组成型数据（CLR 归一化）、缺失条件下的表现。

### B. 实际数据分析（analysis_script/）
- [ ] 接入首选数据集 **Clark 2021**（Zenodo，N=1850，26 物种）：下载、CLR/参考菌归一化、构建 $(X, U)$。
- [ ] 接入 **Venturelli 2018**（OD600 校正，近绝对丰度，12 物种）作为可解释性强的对照。
- [ ] 编写统一的真实数据 PSS 推断 pipeline（数据 → 设计矩阵 → ADSIHT → Jacobian → 网络）。
- [ ] 网络可视化脚本（igraph + ggplot2，遵循配色规范），输出推断互作网络图。
- [ ] 与已发表 gLV/互作结果横向比较（Stein 2013、MDSINE2 作为基准）。

### C. 论文撰写（manuscript/，待建）
- [ ] 方法部分：整理 §1–§7 的公式推导与可识别性论证。
- [ ] 模拟实验部分：汇总 A 的相图与基准表。
- [ ] 实际数据部分：汇总 B 的网络重构结果与生态学解读。
- [ ] 相关工作与贡献定位（相对 SINDy、gLV 回归、Henderson & Michailidis 的改进）。
- [ ] 整理 `ref/references.bib`，统一引用格式。

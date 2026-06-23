# PSS-Net 扰动设计模块：将现有 D-optimal design 接入 PSS 推断

**方法说明** · 当前正式 pilot/oracle 实现：
`sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R`

> **贡献边界**：D-optimality、locally optimal design、Fedorov exchange 和 adaptive design
> 均为已有方法。这里最多应表述为 **A（PSS 稳态响应模型）+ B（已有最优实验设计）** 的组合应用。
> 项目特有部分是把 PSS 扰动到稳态的预测特征交给标准设计后端，并评估它对网络恢复的价值；
> 不应把 D-optimal 算法本身作为新方法贡献。文献与实现边界见
> `ref/pilot_doptimal_literature.md`。

---

## 1. 动机

PSS 数据中外部扰动 $u$ 是可设计变量，因此可以把已有 optimal experimental design 工具用于
选择后续扰动。在固定实验预算 $N$ 下，本模块以 PSS 模型预测候选稳态特征，再调用标准
D-optimal design 后端选择信息量较高的条件。本文需要证明的是这种**已有设计方法在 PSS 网络恢复
场景中的适用性和实验预算收益**，而不是重新提出 D-optimality。

---

## 2. 问题表述

### 2.1 稳态映射

加性 ODE 在稳态 $\dot x=0$：$\mu + f(x^*) + u = 0$。对线性 GLV 真值
$f_{ji}(x)=A_{ji}x,\ f_{jj}(x)=-\gamma_j x$，记 $B=A_\text{off}-\operatorname{diag}(\gamma)$（Hurwitz），
稳态对扰动有**闭式线性映射**：

$$x^*(u) = -B^{-1}(\mu+u) = x^{\mathrm{wt}} - B^{-1}u,\qquad x^{\mathrm{wt}}=-B^{-1}\mu. \tag{1}$$

**关键**：尽管 $x^*$ 对 $u$ 线性，设计矩阵的特征 $\psi_m(x^*)=（x^*)^m$ 对 $u$ **非线性**，
故扰动设计非平凡。

### 2.2 设计矩阵与信息

逐节点回归共享同一设计矩阵 $\Psi\in\mathbb{R}^{N\times pM}$，
$\Psi_{k,(i-1)M+m}=\psi_m(x_i^{*(k)})$；仅响应 $-u_{\cdot j}$ 随节点变化。
**故一次设计服务全部 $p$ 个节点。** 当前实现先对特征列中心化、标准化，并用 pivoted QR
保留可估计列，记所得矩阵为 $\Psi_c$。传给包的数据不手工重复常数列，但显式公式 `~ .`
会按标准线性模型加入一个截距，记 $\tilde\Psi_c=[\mathbf 1,\Psi_c]$。用于精确 D-optimal
augmentation 的信息矩阵为

$$M_N=\tilde\Psi_c^\top\tilde\Psi_c. \tag{2}$$

这里不向 D-optimal 准则加入项目自定义的岭项；岭回归只用于从含噪 pilot PSS 稳定估计
$dx^*/du$，随后把预测特征交给包实现。

### 2.3 设计准则

| 准则 | 目标 | 含义 |
|------|------|------|
| **D-最优** | $\max\ \log\det M_N$ | 最小化 $\hat\theta$ 置信椭球体积 |
| A-最优 | $\min\ \operatorname{tr}(M_N^{-1})$ | 最小化平均方差 |
| E-最优 | $\max\ \lambda_{\min}(M_N)$ | 最坏方向最优 |
| 互相干 | $\min\ \mu(\Psi)$ | 直接控制稀疏恢复（RIP/coherence） |

当前 Fig2e 采用 **pilot-informed exact D-optimal augmentation**：pilot 数据用于估计局部稳态响应
映射；候选特征构造完成后，exact augmentation 交由 CRAN `AlgDesign::optFederov()` 完成。

---

## 3. 当前包实现

对 pilot 和候选扰动预测稳态特征后，先按列中心化、标准化，并用 pivoted QR 去除确定性线性
依赖，只在可估计特征空间内优化。随后调用：

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

- `criterion="D"`：优化 D 准则；
- `augment=TRUE`：把 pilot 条件作为不可交换的 protected runs；
- `nTrials=N_total`：pilot 严格计入总实验预算；
- 其余点由包内 Fedorov exchange algorithm 从有限候选池选择。

每个总预算单独调用包优化，因此设计不强制形成嵌套前缀。项目不维护自写的 D-optimal 优化
核心；PSS 特有代码只负责 `u -> x* -> Phi` 映射。完整算法出处和包引用见
`ref/pilot_doptimal_literature.md` 与 `ref/references.bib`。

---

## 4. 对照基线

1. **random**（现状）：$u\sim\text{Uniform}(-0.4,0.8)^p$，取前 $N$。
2. **maximin（空间填充）**：在 $u$ 空间贪心最大化最小间距（Latin-hypercube 风格）。
3. **D-opt（已有方法）**：`AlgDesign::optFederov()` exact D-optimal augmentation。

---

## 5. 模拟实验方案

- **系统**：$p=8$ 非线性加性 ODE，每节点 2 条入边，其中至少 1 条含二次项。
- **基**：单项式 $M=2$。**回归**：ADSIHT（`ic.type="dsic"`），逐节点 + 中心化 + 标准化。
- **pilot**：$N_{\mathrm{pilot}}\in\{8,12,16\}$，含噪 pilot 用 ridge 估计局部 $dx^*/du$。
- **总预算**：$N_{\mathrm{total}}\in\{20,30,40,60\}$，其中 pilot 计入总数。
- **重复**：每个组合 $R=20$ 个系统种子。
- **指标**：MCC（主）、Precision、Recall、失败节点数。
- **产出**：oracle regret 与 pilot D-optimal 相对 random 的 MCC 增益。

---

## 6. 在论文中的角色

- 作为**应用与系统整合结果**：说明标准 D-optimal augmentation 接入 PSS 稳态预测后，是否能
  减少网络恢复所需实验条件；
- 作为**可行性检查**：量化真 Jacobian oracle 与含噪 pilot plug-in 之间的损失；
- 不将 D-optimality、Fedorov exchange、locally optimal design 或 ridge estimation 声称为本文创新；
- 除非后续提出并证明新的 PSS 专用设计准则，否则该模块最多是 A（PSS-Net）+ B（已有设计方法）。

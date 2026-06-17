# PSS-Net 创新点 A：最优扰动设计（Optimal Perturbation Design）

**方法说明** · 对应实现：`sim_script/02_scaling_design/pss_net_design.R`

---

## 1. 动机

PSS 框架的独占杠杆是**外部扰动 $u$**。现有工作（Xiao 2017、idopNetwork、SINDy、Meister 2013）
要么被动观测稳态，要么随机施扰，**无人形式化"为可识别性该施加哪些扰动"**。
本节把网络重构从"被动推断"升级为"**主动设计 + 推断**"闭环：在固定实验预算 $N$ 下，
选择信息量最大的扰动序列 $\{u^{(k)}\}$，使稀疏非参 ODE 的边恢复质量最大化。

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
**故一次设计服务全部 $p$ 个节点。** 含截距的增广设计 $\tilde\Psi=[\mathbf 1,\Psi]$，
（岭正则）信息矩阵

$$M_N=\tilde\Psi^\top\tilde\Psi+\lambda I\in\mathbb{R}^{(pM+1)\times(pM+1)}. \tag{2}$$

### 2.3 设计准则

| 准则 | 目标 | 含义 |
|------|------|------|
| **D-最优** | $\max\ \log\det M_N$ | 最小化 $\hat\theta$ 置信椭球体积 |
| A-最优 | $\min\ \operatorname{tr}(M_N^{-1})$ | 最小化平均方差 |
| E-最优 | $\max\ \lambda_{\min}(M_N)$ | 最坏方向最优 |
| 互相干 | $\min\ \mu(\Psi)$ | 直接控制稀疏恢复（RIP/coherence） |

本实现采用 **序贯 D-最优（贪心主动学习）**，与稀疏恢复理论（创新点 B）天然衔接。

---

## 3. 序贯 D-最优算法

D-最优的贪心增量有解析形式：向已选集加入特征行 $\tilde\psi_c$ 后

$$\log\det\!\big(M+\tilde\psi_c\tilde\psi_c^\top\big)-\log\det M
  =\log\!\big(1+\tilde\psi_c^\top M^{-1}\tilde\psi_c\big), \tag{3}$$

即**最大化候选行的预测方差** $\tilde\psi_c^\top M^{-1}\tilde\psi_c$（经典序贯 D-最优）。

```
输入：候选扰动池 U_pool（大），预算 N，岭参数 λ
1. 种子：u=0（野生型）+ 少量随机扰动；由 (1) 得 x*，建 Ψ̃，M = Ψ̃ᵀΨ̃ + λI
2. while |selected| < N:
     对每个候选 u_c ∈ U_pool：
        预测 x*(u_c)（闭式 (1)，或对真实非线性系统数值积分）
        构造增广行 ψ̃_c
        score(u_c) = ψ̃_cᵀ M⁻¹ ψ̃_c        # 预测方差，eq.(3)
     选 u* = argmax score；加入 selected；
     秩一更新 M ← M + ψ̃_* ψ̃_*ᵀ，M⁻¹ 用 Sherman–Morrison 更新
3. 输出 {u^(k)} 及对应 {x*^(k)}
```

**复杂度**：每步 $O(|U_\text{pool}|\cdot (pM)^2)$；Sherman–Morrison 免重复求逆。
**主动学习版**：真实系统无闭式 (1) 时，候选打分用线性化 $x^*\approx x^{\mathrm{wt}}-\hat J^{-1}u$，
仅对入选者做高精度积分——对应可在湿实验中序贯执行的自适应扰动设计。

---

## 4. 对照基线

1. **random**（现状）：$u\sim\text{Uniform}(-0.4,0.8)^p$，取前 $N$。
2. **maximin（空间填充）**：在 $u$ 空间贪心最大化最小间距（Latin-hypercube 风格）。
3. **D-opt（本方法）**：序贯 D-最优。

---

## 5. 模拟实验方案

- **系统**：$p=8$ 物种线性 GLV，每节点入边 $\sim2$，强自调节保证 Hurwitz 与 $x^*>0$。
- **基**：单项式 $M=2$。**回归**：ADSIHT（`ic.type="dsic"`），逐节点 + 中心化 + 标准化。
- **关键设置**：在**小预算** $N\in\{12,16,20,30,40,60\}$ 扫描——此区间随机设计欠定/勉强可识别，
  最能体现设计价值（大 $N$ 时各法均饱和）。
- **重复**：每 $(N,\text{strategy})$ 跑 $R=20$ 个种子。
- **指标**：MCC（主）、Precision/Recall、CoefL2、JacRMSE。
- **产出**：学习曲线 MCC–$N$（ggplot2），预期 **D-opt 在小 $N$ 显著占优**
  （达到目标 MCC 所需条件数更少）。输出 `results/design_comparison.csv` 与
  `results/design_mcc_curve.pdf`。

---

## 6. 预期贡献

- 把 PSS-Net 从"又一推断器"升级为 Xiao/idopNetwork/SINDy 都缺失的**设计+推断**维度；
- 序贯 D-最优 + Sherman–Morrison 高效、可解析，且主动学习版可直接落地湿实验；
- 与可识别性/样本复杂度理论（创新点 B）形成闭环。

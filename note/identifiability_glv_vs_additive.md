# Note: 乘性 gLV vs 加性模型的可识别性边界

**用途**：记录"PSS 稳态数据能/不能区分哪种生成机制"的思考与判别策略。
**日期**：2026-06-16。相关：[../methods/sindy_ss_method.md](../methods/sindy_ss_method.md) §1.1、
[../sim_script/01_foundation_recovery/pss_net_glv_ss.R](../sim_script/01_foundation_recovery/pss_net_glv_ss.R)、
[../sim_script/01_foundation_recovery/pss_net_discriminate.R](../sim_script/01_foundation_recovery/pss_net_discriminate.R)。

## 核心结论

PSS 只用扰动稳态 $\{(x^{*(k)}, u^{(k)})\}$，由此：

1. **乘性 gLV vs 加性线性模型：不可识别。**
   二者在内部正稳态满足同一代数方程 $r_j+\sum_{i\neq j}A_{ji}x_i^*-\gamma_j x_j^*+u_j=0$
   （gLV 除以 $x_j$）。`glv_ss_verification.csv` 实测 $\max|x_{\text{gLV}}-x_{\text{lin}}|\approx4\times10^{-16}$。
   ⇒ 任何稳态统计量对两者给出相同拟合/残差/系数，无法区分**动力学形式**。
   要区分需**时序/瞬态数据**（gLV 弛豫率 $\propto x_j$）或**机制先验**。

2. **线性 vs 非线性互作：可识别。**
   通过恢复的 $f_{ji}$ 形状：若仅需一阶项（组内 $s_0=1$、M=1 vs M≥2 残差/DSIC 无显著改善）
   ⇒ 与线性 gLV 相容；若高阶基系数显著（$s_0\ge2$、M≥2 残差显著下降）
   ⇒ 互作非线性，**排除标准 gLV**。

3. **边界灭绝行为：乘性结构的旁证。**
   gLV 中 $x_j=0$ 为吸收边界；加性线性模型在强负扰动下给出负稳态。
   观察到干净灭绝（$x_j\to0$ 锁定）而非负值/线性外推 ⇒ 乘性结构信号。

## 给论文的 caveat 一句话

> PSS 稳态可识别"互作是否非线性"，但**不能**识别动力学是乘性还是加性；
> 二者共享稳态方程，需时序或先验方能区分机制。

## 判别策略（操作）

- 嵌套比较 M=1（线性/gLV 假设）vs M≥2（加二次/样条）：用 DSIC、BIC 或残差平方和下降。
- 看 ADSIHT 选出的组内稀疏 $s_0$。
- 画 $\hat f_{ji}(x_i)$ 看曲率。
- 见判别模拟脚本 [../sim_script/01_foundation_recovery/pss_net_discriminate.R](../sim_script/01_foundation_recovery/pss_net_discriminate.R)
  与 `results/sim_results/discriminate_gof.csv`。

## 判别模拟初步结果（2026-06-16，p=6，N=200，10 seeds）

判别量：逐节点 OLS 的 `relRSS`=(RSS_{M1}−RSS_{M2})/RSS_{M1}、`bicM2`=BIC 选 M=2 的节点比例。

| regime | truth | relRSS | bicM2 | MCC |
|--------|-------|--------|-------|-----|
| narrow ($u\in[-0.3,0.5]$) | glv_lin   | 0.032 | 0.0 | 0.93 |
|                            | add_lin   | 0.031 | 0.0 | 0.89 |
|                            | add_quad  | 0.032 | 0.0 | 0.64 |
|                            | add_monod | 0.032 | 0.0 | 0.72 |
| wide ($u\in[-1.0,2.0]$, 2× 非线性) | glv_lin   | 0.030 | 0.0 | 0.95 |
|                            | add_lin   | 0.030 | 0.0 | 0.91 |
|                            | **add_quad**  | **0.238** | **0.8** | 0.90 |
|                            | add_monod | 0.041 | 0.0 | 0.94 |

稳态差异 glv_lin vs add_lin：narrow $3\times10^{-16}$、wide $3\times10^{-9}$。

**初步结论**：
1. **机制不可识别**：两个 regime 下 glv_lin 与 add_lin 的 relRSS/bicM2/MCC 几乎相同、
   稳态差异≈机器精度 ⇒ 证实乘性/加性不可从稳态区分。
2. **非线性可识别但依赖扰动范围**：narrow 下连真二次互作都查不出（relRSS≈0.03 仅为加
   参数的噪声拟合，BIC 不选 M=2）；wide 下 add_quad 被清晰标出（relRSS 0.24、80% 节点
   BIC→M2），而 glv_lin/add_lin 仍在 0.03。⇒ **检测非线性需足够宽的扰动激发 $x^*$ 曲率**，
   直接呼应最优扰动设计（创新点 A）的"为判别目标设计扰动"。
3. **温和饱和（Monod）更难**：即便 wide 也只 relRSS 0.04、BIC 不选 M=2 ⇒ 缓曲率在单项式
   M=2 下信号弱，检测能力依赖非线性类型与强度（待用样条/更宽扰动改进）。

## 噪声敏感性：BIC 判别非线性在真实数据下会变弱（2026-06-17）

判别"互作是否非线性"靠的是**嵌套模型比较**：在已含全部线性源项之后，加二次基能否把 RSS
降到足以抵过 BIC 罚项 $(\text{ncol}_2-\text{ncol}_1)\log N$。它问的是 $u\!\to\!x^*$ 关系**有无曲率**，
而非绝对拟合优度——所以"调控变量多、线性模型已拟合得很好"并不妨碍判别：线性/gLV 情形二次项
真系数≈0，真非线性情形线性模型会留下**系统性弯曲残差**供二次项压低。

但这套干净分离强烈依赖 **信噪比**。Fig1b（[../sim_script/01_foundation_recovery/Fig1.R](../sim_script/01_foundation_recovery/Fig1.R)）
用的是 $\sigma=0.003$、$N=180\gg16$ 参数的近理想条件，所以线性/gLV 的二次项几乎拿不到任何
虚假信号、add_quad 被干净标出。**真实（微生物组/生态）数据噪声会大得多**：

- 大噪声下线性模型的残差被噪声而非曲率主导，真非线性的 relRSS/BIC gain 被淹没 →
  **假阴性**（漏判非线性）；
- 同时二次基偶然拟合噪声，可能在某些节点给出虚假 BIC gain → **假阳性**；
- 净效应：判别非线性的功效随 $\sigma$ 上升、随 $N$ 下降而衰减，且 Monod 这类缓曲率最先失效
  （见上表已是 narrow 即查不出）。

⇒ **讨论部分应明确**：PSS 对"线性 vs 非线性互作"的识别力是噪声/样本量/扰动幅度的联合函数，
低信噪比真实数据需更大 $N$、更宽扰动（呼应创新点 A 的判别导向扰动设计）、或更稳健的基/惩罚
（样条 + 组稀疏而非单项式 OLS）才能维持。这是从模拟到真实数据的一个核心 caveat。

### Fig1b SNR 扫描补充（2026-06-17）

在当前 `Fig1.R` 的 8 节点设定中，又按 SNR 做了一次辅助分析：先固定无噪声 PSS
矩阵 `X_id_linear`、`X_id_glv`、`X_id_nonlinear`，然后对观测稳态值加噪声

$$
X_{\mathrm{obs}} = X^* + E,\qquad E_{ki}\sim N(0,\sigma_x^2),
$$

其中

$$
\mathrm{SNR}=\frac{\mathrm{mean}_i\{\mathrm{sd}(X^*_{\cdot i})\}}{\sigma_x}.
$$

判别仍然用逐节点嵌套 BIC：线性库 `X` 与二次库 `[X, X^2]` 比较，记录
`bic_gain_per_sample = (BIC_M1 - BIC_M2) / N` 和 `bic_selects_quadratic`。数学上，
二次项能被选中需要它带来的 RSS 下降超过复杂度罚项：

$$
N\log(RSS_1/RSS_2) > (d_2-d_1)\log N.
$$

噪声增大时，`RSS_1` 和 `RSS_2` 同时被测量误差主导，曲率造成的系统性 RSS 下降被稀释，
所以真非线性会首先表现为 **BIC gain 变负**，随后 `bic_selects_quadratic` 接近 0。

当前模拟（`N=180`、二次项强度 `0.42 * sign(A)`、30 次噪声重复）给出：

| SNR | additive linear 选择二次项 | gLV 选择二次项 | nonlinear additive 选择二次项 | nonlinear BIC gain/N |
|-----|----------------------------|----------------|-------------------------------|----------------------|
| 10  | 0.00 | 0.00 | 0.00 | -0.147 |
| 15  | 0.00 | 0.00 | 0.04 | -0.113 |
| 20  | 0.00 | 0.00 | 0.28 | -0.051 |
| 30  | 0.00 | 0.00 | 0.70 |  0.103 |
| 50  | 0.00 | 0.00 | 1.00 |  0.448 |
| 100 | 0.00 | 0.00 | 1.00 |  1.298 |

粗略阈值：

- **SNR < 20**：基本无法稳定识别非线性，BIC gain 为负或接近 0；
- **SNR ≈ 20-30**：过渡区，节点层面的非线性选择率从约 0.28 升到约 0.70；
- **SNR ≥ 50**：当前设定下识别较稳定，二次项几乎总被选择。

线性加性 ODE 与 standard gLV 在所有 SNR 下都保持 0 的二次项选择率，说明这里的失败主要是
**真非线性被噪声淹没的假阴性**，而不是线性机制产生明显假阳性。这个结论依赖当前扰动范围、
非线性强度、`N=180` 和 BIC 罚项；更弱曲率、更窄扰动或更小样本量会把阈值推到更高 SNR。

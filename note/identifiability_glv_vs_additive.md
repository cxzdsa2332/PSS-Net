# Note: 乘性 gLV vs 加性模型的可识别性边界

**用途**：记录"PSS 稳态数据能/不能区分哪种生成机制"的思考与判别策略。
**日期**：2026-06-16。相关：[../methods/sindy_ss_method.md](../methods/sindy_ss_method.md) §1.1、
[../sim_script/pss_net_glv_ss.R](../sim_script/pss_net_glv_ss.R)、
[../sim_script/pss_net_discriminate.R](../sim_script/pss_net_discriminate.R)。

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
- 见判别模拟脚本 [../sim_script/pss_net_discriminate.R](../sim_script/pss_net_discriminate.R)
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

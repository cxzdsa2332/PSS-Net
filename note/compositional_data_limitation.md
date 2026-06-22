# Note: 组成型（相对丰度 / CLR）数据下 PSS 网络恢复的局限

**用途**：记录"PSS 加性稳态回归在组成型数据上为何失效、何种预处理可行"的实验发现。
**日期**：2026-06-16。相关：[../sim_script/03_robustness_benchmarks/Fig3c_compositional_data_limitation.R](../sim_script/03_robustness_benchmarks/Fig3c_compositional_data_limitation.R)、
`results/sim_results/Fig3c_compositional_data_limitation.csv`、[../data/datasets.md](../data/datasets.md)、
[../methods/sindy_ss_method.md](../methods/sindy_ss_method.md) §2。

## 背景

PSS 稳态关系 (2′) $r_j+\sum_i A_{ji}x_i^*-\gamma_j x_j^*+u_j=0$ 写在**绝对**丰度上。
真实 16S 测序只给**相对丰度**（组成型，总量未知）。本实验量化组成型化的破坏并比较预处理。

## 模拟设置

p=10、s=3 入边、N=200、10 seeds；闭式稳态（对角占优、正值）。
含噪绝对丰度 $X_{\rm abs}=x^*\cdot e^{\mathcal N(0,0.05^2)}$；相对 $Z=X_{\rm abs}/\sum X_{\rm abs}$；
CLR$(Z)=\log Z-\overline{\log Z}$；实测总量 $\hat T=T\cdot e^{\mathcal N(0,0.15^2)}$，$X_{\rm hat}=Z\hat T$。
逐节点 ADSIHT、单项式 M=2。

## 结果（mean MCC，10 seeds）

| 输入 | Precision | Recall | MCC |
|------|-----------|--------|-----|
| abs（绝对，oracle 上界） | 0.94 | 0.93 | **0.896** |
| rel（相对，朴素闭合） | 0.73 | 0.61 | **0.522** |
| clr（中心对数比） | 0.52 | 0.59 | **0.300** |
| rel_x_T（相对 × 含噪总量） | 0.33 | 0.60 | **−0.02** |

## 结论

1. **绝对丰度是刚需**：朴素相对丰度掉约 40% MCC（0.90→0.52）。
2. **CLR 反而更糟（0.30）**：CLR 把关系变成 log-ratio，与本方法"绝对 $x$ 上的加性"模型
   失配。CLR 适合 SPIEC-EASI 等协方差/log-ratio 方法，**不适合**加性-$x$ 稳态回归。
3. **靠含噪总量重建绝对会崩（MCC≈0）**：$Z\hat T=x^* e^{\eta_k}$ 引入**逐条件整体尺度误差**，
   且该误差与 $u$ 相关（$\propto(c_k-1)(-u_j-r_j)$），系统性污染回归；15% 总量 CV 即足以摧毁。
   ⚠️ 此项对总量噪声 CV 敏感，数值依设定；定性结论稳健。

## 启示（写入论文 limitations / 未来工作）

- 方法本质需要**低噪声的绝对定量**（OD600 / qPCR 校正），如 Venturelli 2018、Stein 2013
  （绝对丰度），优先于纯组成型 16S。呼应 [../data/datasets.md](../data/datasets.md) 的数据排序。
- 朴素相对丰度或 CLR 都不足以直接用于本方法。
- 未来工作：发展**组成型稳健的 PSS 变体**（如把总量作为隐变量联合估计、或在 log-ratio
  几何下重写稳态约束），是真问题而非工程细节。

# Note: 联合(块对角)估计 vs 逐节点——优劣取决于网络结构

**用途**：记录"何时联合求解优于逐节点"的探索结论。
**日期**：2026-06-17。相关：[../sim_script/pss_net_scalefree.R](../sim_script/pss_net_scalefree.R)
（scale-free 正面）、[../sim_script/pss_net_joint_smalln.R](../sim_script/pss_net_joint_smalln.R)
（同质负对照）、`results/sim_results/scalefree_compare.csv`、CLAUDE.md Method Conventions。

## 设定

- 两种估计器解同一组扰动稳态回归（共享中心化+标准化设计 $\Psi_{cs}$）：
  - **逐节点**：每个目标节点单独跑 ADSIHT（p 次，各 N 样本）；
  - **联合**：块对角 $X=I_p\otimes\Psi_{cs}$，**p·p 组、每组重复 M**（v0.1.txt 规则），
    一次 ADSIHT/DSIC，在 $p\cdot N$ 残差上联合估计稀疏度/噪声。
- 单项式 M=2，p=50，逐边 MCC + hub 识别（估计出度 vs 真出度 Spearman）。

## 关键结论：结构依赖

**scale-free 网络（偏好连接，幂律出度，真出度峰值≈10，5 seeds）—— 联合稳定更优：**

| N | MCC 逐节点 | MCC 联合 | 出度Spearman 逐节点 | 出度Spearman 联合 |
|---|-----------|---------|---------------------|-------------------|
| 16 | 0.070 | 0.073 | 0.01 | **0.21** |
| 24 | 0.196 | **0.231** | 0.25 | **0.36** |
| 40 | 0.379 | **0.415** | 0.43 | **0.46** |

**均匀入度 + 同质信号（负对照）—— 两者打平**：联合 MCC ≈ 或略低于逐节点，
仅表现为 precision↑ / recall↓ 的保守取舍（见 `joint_smalln`）。

## 机理

联合的 `p·p` 分组**不在任务间共享支撑**，只在 DSIC 里池化一个全局稀疏度/噪声水平。
- 在 **scale-free** 下，少数 hub 源出现在**很多**目标的真实支撑中——这是跨目标**反复出现
  的共享弱信号**；全局池化把它**累加**成强信号，于是更稳地选中 hub 边、更早排出 influence
  排名。逐节点各自为战，小 N 下抓不住反复出现的弱信号。
- 在 **均匀/同质** 网络下没有"反复出现的共享源"，无可借之力 → 无增益。

## 实务结论

- **默认逐节点**（快几百倍；联合是稠密 $I_p\otimes\Psi$，p=100 内存不可行）。
- **当网络异质/scale-free（有 hub）且 N 较小时，联合值得用**：逐边 MCC 与 hub 识别都更好。
  真实生物网络多为 scale-free，故联合在真实场景有实际价值。
- 想要更强的联合增益，仍需真正的跨任务支撑共享（多样本/多个体共享网络，轴 A）。

## 备注

- 早期用"源跨目标"分组（组大小 $p\cdot M$）会让 ADSIHT 组内稀疏涌出假阳性、MCC 崩到 ~0；
  必须用 `p·p` 分组。已写入 CLAUDE.md Method Conventions。
- 已删除两个被取代的探索脚本：`pss_net_jointsolve.R`（随机网络）、
  `pss_net_multitask.R`（ad-hoc hub）。

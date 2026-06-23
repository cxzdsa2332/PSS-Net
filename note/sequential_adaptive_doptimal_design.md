# Note: 序贯/自适应 D-optimal 扰动设计——分批重估的取舍

**用途**：记录"先粗略设计几个扰动、再根据结果进一步设计"这一序贯（adaptive）思路相对
当前 two-stage pilot-estimated D-optimal 的机制、收益、代价与方法边界，供 Fig2e 后续扩展或
discussion 引用。暂不实现模拟。
**日期**：2026-06-23。

相关：[../methods/optimal_perturbation_design.md](../methods/optimal_perturbation_design.md)、
[../ref/pilot_doptimal_literature.md](../ref/pilot_doptimal_literature.md)、
模拟脚本 [../sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R](../sim_script/02_scaling_design/Fig2e_oracle_vs_estimated_design.R)、
规划 [../sim_script/simulation_plan.md](../sim_script/simulation_plan.md)、
另一条设计 note [targeted_perturbation_tradeoff.md](targeted_perturbation_tradeoff.md)。

## 当前 Fig2e 是 two-stage（K = 1），不是多轮

现状（见 `Fig2e_oracle_vs_estimated_design.R` 主循环）：

1. 随机 pilot `pilot_n`（8/12/16）个扰动 `u`，首点为 WT（`u = 0`），测稳态（含噪 σ=0.04）。
2. 从 noisy pilot 用 ridge 估计局部响应映射 `Ĥ ≈ dx*/du`（`predict_pool_from_pilot`，λ=0.25）。
3. 用 `Ĥ` 把候选池每个 `u` 预测成 PSS 特征 `[x, x²]`，**一次性**调用
   `AlgDesign::optFederov(criterion="D", augment=TRUE, rows=pilot)` 把剩余 `N_total − pilot_n`
   个 `u` 全部选完，pilot 作为不可交换的 protected runs。

关键：第 3 步只用最初那批 pilot 的估计就把全部后续 `u` 定死，**中途不回看结果**。
`N_adaptive` 列只是 `n_add`。这是 two-stage exact augmentation，不是在线 adaptive Wynn。

## pilot-estimated D-optimal 如何影响实际 `u`

`Ĥ` 决定"信息在 `u` 空间哪个方向最缺"。optFederov 最大化 `|ΨᵀΨ|`，于是把 `u` 放在让 PSS 特征
`[x, x²]` 张得最开的位置（通常靠近 box 边界/角，且沿 `Ĥ` 映射后曲率大的方向）。因此：

- `Ĥ` 准 → 选出的 `u` 接近 oracle 的选择；
- `Ĥ` 差（pilot 太小）→ 选点跑偏，这就是 Fig2e 量化的**可行设计损失**，在 `N_total=20`、
  pilot 较大时最严重（把预算过早压在烂估计上，反而不如 model-free 的 maximin）。

## 序贯/自适应版本 = batched Fedorov–Wynn

把第 2–3 步做成多轮循环即得：

```
随机 pilot → 测 → 估计 Ĥ → 设计下一小批 u → 测 → 重估 Ĥ → 再设计 …（K 轮）
```

当前 two-stage 是它 K = 1 的特例。这在文献里是 batched/sequential（adaptive）optimal design，
也更贴近真实湿实验（本就分批做）。

### 预期收益
- 每轮用更多数据重估 `Ĥ` 再决定下一批，能吃掉一部分小预算区的可行损失——adaptive 会自我纠偏，
  two-stage 把一次烂估计锁死。
- 收益集中在**小预算/早期**：Fig2e 已显示预算一大，pilot ≈ oracle，自适应边际递减。

### 代价 / 权衡
- **轮次成本**：每轮都要等系统到稳态再测，湿实验周转/批次开销高；轮数 K 与每批大小是核心旋钮。
- **探索-利用**：早期既要"探清 `Ĥ`"又要"对当前估计最优"，纯贪心 myopic 批选择可能偏向利用、
  欠探索，需要权衡（甚至前几轮保留一定 space-filling 成分）。
- **边际递减**：N 大时不划算；K 过大只是徒增实验轮次。
- **稳健性**：早期单批烂估计仍可能误导一轮，但多轮能自我修正，这点优于 two-stage。

## 方法边界（为什么默认仍是 two-stage）

- two-stage 的选点后端正好落在 CRAN `AlgDesign::optFederov()`，方法创新边界干净
  （A：PSS 特征构造 + 网络恢复评估；B：已有最优设计），见
  [../ref/pilot_doptimal_literature.md](../ref/pilot_doptimal_literature.md)。
- 在线 adaptive Wynn 需要自写多轮重估循环，会把"已有工具"边界往外推，需在文中明确这是
  feasible adaptive baseline 而非新算法。
- 历史：本条线最早的 Fig2e 曾是"每 6 点重估一次"的批式 adaptive，后重构为现在的干净 two-stage。

## 决策

- **默认保持 two-stage 作为 Fig2e 主面板**（干净、可复现、方法边界清楚）。
- 序贯/自适应作为**可选扩展或补充**：若要落地，建议对照
  `two-stage pilot D-opt vs K-轮 adaptive D-opt vs oracle vs random/maximin`，
  主轴放在小预算区"分批自适应挽回多少可行损失"；计算较重（多轮 ODE）。
- 进入主图前需先决定它是主文 Fig2 子面板还是 `sup/` 补充。

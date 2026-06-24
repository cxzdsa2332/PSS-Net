# Scale-free external benchmark: exploratory finding

Purpose: record what a scale-free / hub variant of the external benchmark
(`Fig3b_external_benchmark_main.R`) showed, and why we did **not** adopt it as a
figure panel. The simulation code was a one-off and has been removed; this note
preserves the conclusion so we do not re-run it blindly.
Date: 2026-06-24.

## What was tried

A drop-in clone of `Fig3b_external_benchmark_main.R` with the **only** change
being the network generator: `make_system` (fixed in-degree `s_in = 2`, uniform
random sources) replaced by a preferential-attachment generator
(`make_system_scalefree`) that gives a few high-out-degree hub sources. Expected
in-degree was held at `s_in` so edge density and the `N / (s log p)` budget axis
stayed comparable. Everything else — methods, metrics, regimes
(`linear`, `strong_nonlinear`), budget grid `N ≈ {69, 103, 171, 273}`,
`p = 30`, `R = 10` seeds — was identical to the homogeneous benchmark.

## Main result

Scale-free numbers are **nearly identical to homogeneous** at every budget.
Representative budget `N = 171` (`N / [s log p] ≈ 25`):

| method | regime | MCC | AUPRC | FuncNRMSE |
|---|---|---|---|---|
| PSS-Net | linear | 0.90 | 1.00 | **0.031** |
| Lasso | linear | 0.82 | 1.00 | 0.089 |
| aiMeRA | linear | **1.00** | 1.00 | N/A |
| PSS-Net | strong-nl | 0.86 | 0.97 | **0.092** |
| Lasso | strong-nl | 0.69 | 0.99 | 0.173 |
| aiMeRA | strong-nl | 0.93 | 0.99 | N/A |

MCC vs N (all methods, both regimes) showed the same ordering as homogeneous:
aiMeRA tops MCC (linear 1.00 across N; strong-nl 0.90→0.95), PSS-Net stable at
0.85–0.91, Lasso the strongest same-library competitor (0.65–0.88),
ElasticNet / GroupLasso / association methods / PySINDy clearly lower. PSS-Net
keeps its decisive **FuncNRMSE** lead (2–3× better than Lasso); aiMeRA / Cor /
PartialCor have no edge function (N/A).

## Why it is not a panel (and not a replacement for Fig3b)

1. **Zero information gain as a headline.** Ranking, magnitude and winners are
   unchanged from the homogeneous benchmark, so swapping topology re-draws the
   same story while discarding the value of the neutral (no-topology-assumption)
   baseline. Not cherry-picked for PSS-Net either — aiMeRA leads MCC on both.
2. **The hub "selling point" was never exercised here.** This run had no
   joint-vs-node-wise PSS-Net variant, no hub-recovery metric, and the budget
   grid stayed at `N > p = 30`, so MRA's response inversion never became
   ill-posed. The regime where hubs matter — CLAUDE.md's note that the joint
   block-diagonal solve beats node-wise on scale-free at small `N` (e.g.
   `p = 50, N = 24`), and where MRA (needs `N > p`) degrades — was not tested.

## What it would take to earn a scale-free panel

Before scale-free deserves figure space it must carry something homogeneous does
not. Minimum additions: (a) PSS-Net **joint vs node-wise** variants; (b) budget
extended into **`N < p`** small-sample regime; (c) **hub-recovery** metric
(out-degree Spearman, top-k hub hit rate); (d) competitor field trimmed to the
mechanism-isolating set (same-library Lasso / GroupLasso + MRA where `N > p`),
dropping the undirected association methods. See
`sim_script/03_robustness_benchmarks/pss_net_scalefree.R` (joint-vs-node-wise on
hubs) and the `Fig3d_scalefree_structure_dependence` row of
`sim_script/simulation_plan.md`, which is the intended home for this analysis.

# External Benchmark Methods for PSS-Net

Purpose: collect candidate external methods for publication-level PSS-Net simulations. The benchmark must separate scientific targets: directed perturbation-response coupling, directed predictive regulator ranking, sparse equation discovery, and undirected association networks.

## Method Classes

### 1. Same Steady-State Equation, Weaker Sparsity Structure

These are the fairest statistical baselines because they can use the same response `-u_j` and the same library `Psi(X)`.

- Lasso: node-wise sparse regression with element-wise sparsity, no source-level grouping.
- Elastic net: node-wise sparse regression with correlated-feature stabilization.
- Group lasso: source-level group sparsity, but no double sparsity within each edge function.
- STLSQ / SINDy-style library regression: sparse thresholded regression on the same candidate library.

Recommended role: main benchmark for demonstrating the value of double sparsity.

Citation keys: `friedman2010glmnet` for lasso / elastic net via `glmnet`;
`yang2015groupbmd` and `lounici2011group` for group-lasso implementations and
theory; `brunton2016sindy` and `kaptanoglu2022pysindy` for SINDy / PySINDy.

### 2. Black-Box Directed Network Inference from State Data

These methods typically use expression/state data only and do not use perturbation inputs. They provide an important bioinformatics comparator but answer a different question.

- GENIE3: tree-ensemble regulator ranking; directed by treating each target as a supervised prediction problem.
- GRNBoost2 / Arboreto: scalable gradient-boosting-style GRN inference.
- dynGENIE3: time-series variant; include only if time-resolved simulated trajectories are generated.

Recommended role: show how much perturbation input improves directed coupling recovery compared with state-only prediction.

Citation keys: `huynhthu2010genie3` for GENIE3 and `moerman2019grnboost2` for
GRNBoost2 / Arboreto.

### 3. Dynamical-System Discovery Methods

These methods infer equations from time-series derivatives or derivative surrogates. PSS-Net avoids numerical differentiation by using steady states, so comparisons require careful data matching.

- SINDy / STLSQ: sparse identification with polynomial or custom libraries.
- SINDYc: SINDy with external inputs/control.
- Implicit-SINDy or SINDy-PI: relevant for rational or implicit dynamics, but may be overkill for main simulations.

Recommended role: compare only in a clearly defined setting, such as simulated time series versus PSS data under equal or documented measurement budgets.

### 4. Association and Compositional Network Methods

These methods usually infer undirected associations rather than directed perturbation-response coupling.

- Pearson/Spearman correlation: minimal association baseline.
- Partial correlation or graphical lasso: conditional association baseline.
- SPIEC-EASI: sparse inverse covariance framework for compositional microbiome data.
- SparCC: correlation estimation for compositional microbial data.

Recommended role: compositional or microbiome-style stress tests; report that the target is association, not directed causal coupling.

Citation keys: `kim2015ppcor` for partial correlation implementation,
`kurtz2015sparse` for SPIEC-EASI, and `friedman2012sparcc` for SparCC.

### 5. Perturbation-Using Directed / Causal Inference

These methods, like PSS-Net, exploit that the perturbation input `u` is known and
controlled. They are the most direct competitors for the central PSS-Net claim
("using `u` recovers directed coupling") and were under-represented in the earlier
set. Self-ablations of PSS-Net are intentionally excluded here; these are external
methods only.

- Modular Response Analysis (MRA): the canonical method for inferring a directed,
  row-normalized local-response matrix from steady-state perturbation responses.
  The benchmark now calls the published `aiMeRA` R package. Because the simulation
  uses continuous multivariate perturbations rather than one intervention per
  module, it first estimates the square global response `dX/dU` from all conditions
  and then calls `aiMeRA::mra(..., Rp=TRUE)`. The former self-contained
  `LinearPSS` regression has been removed from the formal benchmark.
- Interventional causal discovery (GIES): greedy interventional equivalence search
  that consumes intervention/knockout data and returns a directed graph. Faithful
  comparator for the causal-discovery community. Caveat: GIES (and DAG learners such
  as NOTEARS-MLP / DAG-GNN) assume acyclicity, whereas PSS networks allow feedback
  loops; include partly to show PSS-Net represents cyclic feedback that DAG methods
  cannot.
- SINDy with control input (SINDYc): sparse identification with the control `u`. On
  PSS data the steady-state algebraic equation `0 = Theta(x) xi + u` makes SINDYc
  with `u` as the response mathematically coincide with steady-state STLSQ on the
  same library (Class 1); the time-series-derivative SINDYc is a different-data
  variant that needs simulated trajectories and belongs in a data-matched supplement.

Recommended role: head-to-head directed-coupling comparators that also use `u`;
report alongside the linear / nonlinear truth split so the linear-only methods (MRA)
are judged where they are valid and where they break.

Citation keys: `kholodenko2002untangling` for the classical MRA response-analysis
framework, `jimenezdominguez2021aimera` for the package-backed aiMeRA implementation,
and `hauser2015gies` for GIES.

## Default-Parameter Policy

All benchmark methods are run with package defaults or a single documented default
threshold; no per-method hyperparameter search is performed. Where a method has a
built-in selection path (cross-validation, information criterion, sequential
thresholding) the default rule is used; dense-score methods (correlation, partial
correlation, MRA) use a fixed significance/FDR rule documented in the script. The
threshold-free ranking metrics (AUROC, AUPRC) are the primary fair comparison; MCC
and friends are reported under each method's documented selection. In the main
PSS-Net versus aiMeRA comparison, both methods use the absolute row-normalized
local-response score and the fixed `> 0.05` cutoff. PSS-Net's minimum-DSIC `A_out`
support is retained as an audited native-selection sensitivity result rather than
the headline graph. The shared equation-regression library is currently
`[x, x^2]`, matching the maximum order
in the quadratic simulation truth. Higher-order library robustness is treated as
a separate sensitivity analysis rather than mixed into the main benchmark.

## Environment Availability (current workspace)

Cleanly available R packages: `ADSIHT`, `glmnet`, `gglasso`, `deSolve`, `igraph`,
`aiMeRA` 0.99.0, `reticulate`, and `ppcor` 1.1. The script now requires the
original-author Python `pysindy` package at commit
`c4421fcec275c8f4cc5c1e93bebb961b212067ae`. A standalone Python 3.13 runtime and
project venv are installed under `.python-fig3b/` and `.venv-fig3b-standalone/`;
the latter is the script default and can be overridden with `PSSNET_PYTHON`. Not installed:
`grpreg` (use `gglasso` for group lasso), `GENIE3` /
`randomForest` (state-only tree ensemble -- wired but skipped until installed), and
`pcalg` (GIES -- wired but skipped until installed). Partial correlation now calls
`ppcor::pcor()`; Pearson correlation uses base R `stats` functions.

PySINDy setup is intentionally explicit rather than allowing `reticulate` to
silently create an environment. The benchmark script automatically discovers
`.venv-fig3b-standalone`, `.venv-fig3b`, or `.python-fig3b` under the project
root, using `Scripts/python.exe` on Windows and `bin/python` on macOS/Linux.
`PSSNET_PYTHON` remains the highest-priority override.

macOS/Linux setup from the project root:

```bash
python3 -m venv .venv-fig3b
./.venv-fig3b/bin/python -m pip install --upgrade pip
./.venv-fig3b/bin/python -m pip install -r requirements/fig3b-pysindy.txt
FIG3B_R=1 FIG3B_K=10 Rscript sim_script/03_robustness_benchmarks/Fig3b_external_benchmark_main.R
```

Windows PowerShell setup from the project root:

```pwsh
py -3 -m venv .venv-fig3b
.\.venv-fig3b\Scripts\python.exe -m pip install --upgrade pip
.\.venv-fig3b\Scripts\python.exe -m pip install -r requirements\fig3b-pysindy.txt
$env:FIG3B_R = "1"
$env:FIG3B_K = "10"
Rscript sim_script/03_robustness_benchmarks/Fig3b_external_benchmark_main.R
```

The environment variable is only needed for a nonstandard location:

- macOS/Linux: `PSSNET_PYTHON=/path/to/python Rscript ...`
- Windows PowerShell: `$env:PSSNET_PYTHON = "C:\path\to\python.exe"`

For exact reproduction of the currently audited package source:

```r
remotes::install_github("bioinfo-ircm/aiMeRA@86cabc21e8ed124ce372c2fe8e62b47503c2a22b")
```

## Reporting Requirements

For each benchmark method, report:

- inputs used: `X`, `u`, time series, total biomass, count table;
- output type: directed edge ranking, undirected association, signed/unsigned edges, function estimates;
- tuning protocol: cross-validation, information criterion, default package settings, or oracle threshold;
- metrics: MCC for selected networks, AUROC/AUPR for rankings, sign accuracy where
  available, total/linear/nonlinear local-Jacobian RMSE and edge-function
  RMSE/NRMSE where coefficients are available, and hub recovery where relevant;
- runtime and failure rate for larger p.

## Minimum Main-Text Benchmark Set

Run over a linear and a nonlinear truth regime (shared `A`, `r` and perturbation
design, differing only by the quadratic term), so linear-only methods are judged
both where they are valid and where they break. A compact main-text set should
include:

1. PSS-Net / ADSIHT (node-wise double sparsity).
2. Lasso and elastic net on the same steady-state library.
3. Group lasso on the same steady-state library (`gglasso` here, not `grpreg`).
4. `PySINDy_STLSQ`: official PySINDy optimizer on the same steady-state library.
   This is PSS algebraic STLSQ, not time-series SINDYc.
5. `aiMeRA`: package-backed classical normalized Modular Response Analysis.
6. Correlation and package-backed `ppcor` partial correlation as association baselines.
7. (When installed) GENIE3 / GRNBoost2 as state-only directed black-box inference,
   and GIES as interventional causal discovery.

SPIEC-EASI and SparCC should be included when presenting microbiome/compositional
simulations or real-data case studies. The main benchmark script is
`sim_script/03_robustness_benchmarks/Fig3b_external_benchmark_main.R`; it reuses the
linear/nonlinear simulation and the node-wise PSS-Net inference from
`sim_script/01_foundation_recovery/Fig1c_adsiht_group_lasso_scaling.R`.

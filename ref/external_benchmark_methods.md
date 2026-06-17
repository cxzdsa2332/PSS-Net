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

### 2. Black-Box Directed Network Inference from State Data

These methods typically use expression/state data only and do not use perturbation inputs. They provide an important bioinformatics comparator but answer a different question.

- GENIE3: tree-ensemble regulator ranking; directed by treating each target as a supervised prediction problem.
- GRNBoost2 / Arboreto: scalable gradient-boosting-style GRN inference.
- dynGENIE3: time-series variant; include only if time-resolved simulated trajectories are generated.

Recommended role: show how much perturbation input improves directed coupling recovery compared with state-only prediction.

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

## Reporting Requirements

For each benchmark method, report:

- inputs used: `X`, `u`, time series, total biomass, count table;
- output type: directed edge ranking, undirected association, signed/unsigned edges, function estimates;
- tuning protocol: cross-validation, information criterion, default package settings, or oracle threshold;
- metrics: MCC for selected networks, AUROC/AUPR for rankings, sign accuracy where available, hub recovery where relevant;
- runtime and failure rate for larger p.

## Minimum Main-Text Benchmark Set

A compact main-text set should include:

1. PSS-Net / ADSIHT.
2. Lasso or elastic net on the same steady-state library.
3. Group lasso on the same steady-state library.
4. STLSQ/SINDy-style sparse library regression on the same steady-state library.
5. GENIE3 or GRNBoost2 as state-only directed black-box inference.
6. Correlation or graphical lasso as a minimal association baseline.

SPIEC-EASI and SparCC should be included when presenting microbiome/compositional simulations or real-data case studies.

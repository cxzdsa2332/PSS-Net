# PSS-Net

**Perturbed steady-state inference for sparse nonlinear dynamical systems**

PSS-Net is a statistical workflow for learning sparse directed coupling networks
from perturbed steady-state observations of complex nonlinear systems.  The core
idea is to use equilibrium responses to interventions, rather than noisy
time-derivative estimates, to fit sparse additive nonparametric ODE models.

Ecological and microbiome dynamics are an important motivating application, and
generalized Lotka--Volterra (gLV) systems are used as a key simulation benchmark.
They are examples, not the sole scope of the project.

## Method Overview

For $p$ coupled state variables, PSS-Net uses the working model

$$
\frac{dx_j(t)}{dt}
=
\mu_j + f_{jj}(x_j(t)) + \sum_{i\neq j} f_{ji}(x_i(t)) + u_j(t),
\qquad j=1,\ldots,p.
$$

At a perturbed steady state, $\dot x_j=0$, giving an algebraic regression
constraint:

$$
\mu_j + f_{jj}(x_j^*) + \sum_{i\neq j} f_{ji}(x_i^*) + u_j = 0.
$$

The main estimation steps are:

1. Build a no-intercept basis expansion for each univariate effect $f_{ji}$.
2. Stack perturbed steady-state observations into node-wise regressions.
3. Center and scale the design matrix.
4. Estimate sparse grouped coefficients with ADSIHT.
5. Recover directed edges from group norms and local effects from the Jacobian.

Detailed mathematical notes are in `methods/sindy_ss_method.md`; the manuscript
method section is in `manuscript/sections/02_method.tex`.

## Directory Structure

| Path | Purpose |
|------|---------|
| `sim_script/` | Formal simulation and inference scripts. Scripts here generate numeric CSV outputs only. |
| `sim_script/manual/` | Exploratory, one-off, or historical scripts. These may mix calculation and plotting and are not formal reproduction entry points. |
| `analysis_script/` | Plotting and summary scripts that read `results/sim_results/` and write figures/tables. |
| `methods/` | Method notes and design documents. |
| `manuscript/` | LaTeX manuscript project. |
| `data/` | Dataset notes and future real-data inputs. |
| `ref/` | References, BibTeX, and literature notes. |
| `results/` | Local generated outputs. This directory is ignored by Git for now. |
| `CLAUDE.md` | Project-level agent instructions. |
| `.github/copilot-instructions.md` | Copilot mirror of `CLAUDE.md`; keep identical. |

## Reproducible Workflow

Formal simulations write only numeric outputs:

```bash
Rscript sim_script/01_foundation_recovery/pss_net_compare.R
Rscript sim_script/02_scaling_design/pss_net_design.R
Rscript sim_script/02_scaling_design/pss_net_design_nl.R
Rscript sim_script/02_scaling_design/pss_net_design_nl_seq.R
Rscript sim_script/01_foundation_recovery/pss_net_glv_ss.R
```

Analysis scripts read those outputs and create figures/tables:

```bash
Rscript analysis_script/plot_design_curves.R
Rscript analysis_script/summarize_mcc_comparison.R
```

Generated outputs are written under:

- `results/sim_results/`
- `results/figure/`
- `results/table/`

`results/` is currently not version controlled. Regenerate results from scripts
when needed.

## Technical Stack

- Main language: R
- ODE solving: `deSolve`
- Sparse estimation: `ADSIHT`, with `grpreg` group lasso as a baseline
- Plotting: `ggplot2`
- Manuscript: LaTeX; `tectonic main.tex` works in the current environment

## Current Status

- Core perturbed steady-state regression workflow is implemented.
- Analysis and plotting have been separated from formal simulation scripts.
- gLV is treated as a simulation benchmark demonstrating how a multiplicative
  nonlinear system can reduce to an additive linear steady-state equation.
- Real-data analysis remains future work.

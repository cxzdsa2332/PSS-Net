# PSS-Net

## Project Scope

PSS-Net is a project for sparse network inference in complex nonlinear dynamical
systems from perturbed steady-state data.  The main statistical object is a
directed coupling network among state variables, estimated through sparse
additive nonparametric ODE models and double-sparse regression.

Ecological and microbiome systems are important motivating applications.  The
generalized Lotka--Volterra (gLV) model is a key simulation benchmark and
example, not the sole definition of the project.

## Agent Instruction Synchronization

- `CLAUDE.md` and `.github/copilot-instructions.md` are the same project-level
  agent instruction set.
- Whenever either file is modified, update the other file in the same change.
- After editing, run `cmp -s CLAUDE.md .github/copilot-instructions.md` to
  confirm that they are identical.

## Directory Policy

```
PSS-Net/
├── data/             # Dataset notes and future real-data inputs
├── sim_script/       # Formal simulation/inference scripts; numeric outputs only
├── sim_script/manual/# Exploratory or historical scripts; exempt from split rule
├── analysis_script/  # Plotting and summary scripts reading results/sim_results/
├── methods/          # Method notes and design documents
├── manuscript/       # LaTeX manuscript
├── ref/              # References, BibTeX, and literature notes
├── results/          # Local generated outputs; ignored by Git for now
├── CLAUDE.md
└── .github/copilot-instructions.md
```

## Reproducible Analysis Rules

- Formal scripts in `sim_script/` must only run simulations/inference and write
  numeric intermediate data to `results/sim_results/`.
- Formal scripts must not create figures, manuscript tables, or other styled
  outputs directly.
- Plotting scripts in `analysis_script/` must read existing files from
  `results/sim_results/` and write figures to `results/figure/`.
- Table/summary scripts in `analysis_script/` must read existing files from
  `results/sim_results/` and write tables to `results/table/`.
- Exploratory, one-off, or historical scripts belong in `sim_script/manual/`.
  Files in that directory may mix calculation, plotting, and interactive checks,
  but they are not official reproduction entry points.
- If an exploratory script becomes part of the formal workflow, move it out of
  `sim_script/manual/` and split it into computation and analysis scripts.
- `results/` is ignored by Git for now. Do not rely on committed results for
  reproducibility; regenerate outputs from scripts.

## Coding Conventions

- Use R for simulation, inference, and plotting unless there is a clear reason
  to use another language.
- Use `snake_case` for R script names and variables.
- Put a short header at the top of every formal script stating purpose, inputs,
  and outputs.
- Use `ggplot2` for plots.
- Use `deSolve` for R-based ODE simulation.
- Prefer structured CSV outputs for intermediate simulation results.

## Method Conventions

- The primary model is a sparse additive ODE for coupled state variables:
  self-feedback terms are estimated but are not counted as directed cross-node
  edges.
- Work node-wise: build one design matrix per target node and avoid unnecessary
  block-diagonal mega-design matrices in implementation.
- Use no-intercept basis functions so that `f_ji(0)=0` is respected.
- Center and scale the design matrix before sparse regression; recover the
  intercept from the uncentered steady-state equation.
- ADSIHT is the preferred estimator because it supports group sparsity across
  source nodes and within-group sparsity across basis functions.
- Group lasso is a baseline comparator, not the default estimator.
- Local effect signs should be interpreted through the Jacobian at an
  unperturbed reference steady state.

## Manuscript Conventions

- The manuscript's main positioning should be complex systems / nonlinear
  dynamics.  Ecology and microbiome examples should be presented as important
  applications or simulation benchmarks.
- gLV belongs primarily in simulation or examples, not as the sole method model.
- Keep notation consistent across `methods/`, `manuscript/preamble.tex`, and
  `manuscript/sections/02_method.tex`.
- Use `\citep` and `\citet` consistently with `natbib`.
- For biological examples, italicize Latin species names and give abbreviations
  at first use.

## Current Formal Entry Points

Formal simulation scripts:

- `sim_script/pss_net_compare.R`
- `sim_script/pss_net_design.R`
- `sim_script/pss_net_design_nl.R`
- `sim_script/pss_net_design_nl_seq.R`
- `sim_script/pss_net_glv_ss.R`

Formal analysis scripts:

- `analysis_script/plot_design_curves.R`
- `analysis_script/summarize_mcc_comparison.R`

## Build Notes

- The manuscript Makefile assumes `latexmk`.
- In the current environment, `tectonic main.tex` is the available compiler path.
- Generated LaTeX artifacts and `manuscript/main.pdf` are ignored by Git.

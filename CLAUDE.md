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
├── note/             # Project-process thinking: open questions, decisions, caveats
├── manuscript/       # LaTeX manuscript
├── ref/              # References, BibTeX, and literature notes
├── results/          # Local generated outputs; ignored by Git for now
├── CLAUDE.md
└── .github/copilot-instructions.md
```

### Notes (`note/`)

- `note/` holds **project-process thinking**: open questions, design decisions,
  identifiability caveats, dead ends, and "why we chose X" rationale that arises
  during the work but is not a finished method spec.
- Distinct from `methods/` (polished method/design documents) and `ref/`
  (external literature). A note may later graduate into `methods/` or the
  manuscript once settled.
- One topic per Markdown file, `snake_case` filename, with a short purpose line
  at the top and the date. Link related notes/methods by relative path.

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
  - Related prior art: the multi-task block-diagonal construction in
    `dong2026multitask` (Communications Physics 2026,
    https://www.nature.com/articles/s42005-026-02687-4) and the exploratory
    `ref/v0.1.txt` both assemble one large block-diagonal design matrix and solve
    all nodes jointly. That is mathematically equivalent to the node-wise loop;
    if reusing that style, the key caveat is **centering** (remove per-column /
    intercept means so `f_ji(0)=0` and the steady-state equation stay
    identifiable). We default to the node-wise loop for efficiency.
  - Fixed grouping for the joint block-diagonal solve: with the stacked design
    `X = I_p ⊗ Psi_cs`, the ADSIHT group vector must be **`p*p` groups, each
    repeated `M` times** — one group per `(target, source)` block, i.e.
    `group = rep(1:(p*p), each = M)` (the v0.1.txt / `dong2026multitask` scheme).
    Do NOT group a source across all targets (`rep(rep(1:p, each=M), times=p)`):
    that makes each group size `p*M`, and ADSIHT's within-group sparsity then
    floods false positives (empirically MCC collapses to ~0).
  - Empirical note (`sim_script/pss_net_scalefree.R` vs
    `pss_net_joint_smalln.R`; see `note/joint_vs_nodewise_structure.md`): whether
    the joint solve (correct `p*p` grouping) beats node-wise is
    **structure-dependent**. On **scale-free / hub networks** joint is
    consistently better in both edge MCC and hub identification (out-degree rank
    correlation), because global DSIC accumulates the repeated weak signal of hub
    sources shared across many target tasks (e.g. p=50, N=24: MCC 0.20 → 0.23,
    hub Spearman 0.25 → 0.36). On **homogeneous / uniform-degree networks** the
    two tie (joint only trades higher precision for lower recall). Joint is far
    more expensive (dense `I_p ⊗ Psi`; infeasible memory at p=100), so node-wise
    stays the default; prefer joint when the network is heterogeneous (hubs) and
    `N` is small.
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
- Target-journal formatting standard: follow the conventions of *Nature* (and
  Nature-family journals such as *Nature Communications* / *Nature Methods*) and
  of *Bioinformatics* (Oxford). When their conventions differ, default to the
  *Nature* style for the main text and figures, and note the *Bioinformatics*
  alternative where relevant. Concretely:
  - Concise, results-forward abstract; structured numbered references in the
    journal's citation style (author–year here via `natbib`, switch to the
    journal style at submission).
  - Main text kept short with methods/derivations in a Methods or Supplementary
    section, per Nature-family structure.
  - Figures: self-contained captions, panel labels (a, b, c), SI units, sans-serif
    figure fonts; avoid chartjunk.
  - Follow each journal's word/figure limits and abbreviation rules; define every
    abbreviation at first use in both abstract and main text.
  See the `manuscript-style` skill for the enforced checklist.

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

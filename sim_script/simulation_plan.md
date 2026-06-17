# Simulation Plan for PSS-Net

Purpose: organize formal simulation scripts around three publication-level simulation goals. The three folders below correspond to the three main multi-panel figures envisioned for the manuscript. Each goal should eventually include external method comparisons where appropriate.

## Goal 1: Foundation, Identifiability, and Baseline Recovery

Folder: `sim_script/01_foundation_recovery/`

This goal establishes what perturbed steady-state (PSS) data can identify and whether the core estimator recovers directed coupling networks under controlled settings. The central message is that PSS-Net recovers interpretable steady-state coupling functions, while the dynamical mechanism itself has identifiable boundaries.

Included scripts:

- `pss_net_glv_ss.R`: verifies that positive steady states of multiplicative gLV obey the same algebraic steady-state equation used by PSS-Net.
- `pss_net_discriminate.R`: tests which mechanisms can and cannot be distinguished from steady states, especially multiplicative gLV versus additive linear dynamics and linear versus nonlinear interactions.
- `pss_net_compare.R`: compares the current double-sparse ADSIHT estimator with group lasso on the baseline benchmark.

Main figure concept: Figure 1, foundation and identifiability.

Suggested panels:

- Panel 1a: Conceptual PSS-Net workflow from perturbation data to interpretable directed coupling network.
- Panel 1b: Agreement between multiplicative gLV steady states and the algebraic PSS steady-state equation.
- Panel 1c: ADSIHT versus group lasso for edge MCC, precision, recall, and coefficient/Jacobian error.
- Panel 1d: Goodness-of-fit comparison for M=1 versus M=2 under linear, quadratic, and saturating interactions.
- Panel 1e: Representative true versus inferred network, highlighting TP, FP, FN, sign, and core regulators.

## Goal 2: Sample Complexity and Perturbation Design

Folder: `sim_script/02_scaling_design/`

This goal tests whether recovery follows sparse high-dimensional scaling and whether experiment design can reduce the number of perturbation conditions needed for reliable network inference. This is the main simulation axis for the design contribution of the manuscript.

Included scripts:

- `pss_net_highdim.R`: scans p and N to test recovery as a function of the rescaled budget `N / (s log p)`.
- `pss_net_design.R`: compares random, maximin, and sequential D-optimal designs in a linear benchmark.
- `pss_net_design_nl.R`: repeats perturbation-design comparisons under nonlinear quadratic interactions.
- `pss_net_design_nl_seq.R`: stresses the design comparison under stronger nonlinear interactions.

Main figure concept: Figure 2, scaling and design.

Suggested panels:

- Panel 2a: MCC versus `N / (s log p)` for multiple dimensions p, testing sparse recovery scaling.
- Panel 2b: Linear benchmark design curves comparing random, maximin, and D-optimal perturbations.
- Panel 2c: Nonlinear benchmark design curves showing why feature-space information matters beyond uniform spacing in perturbation space.
- Panel 2d: Strong nonlinear benchmark design curves to identify the robustness and limits of linearized D-optimal scoring.
- Panel 2e: Targeted versus non-targeted perturbation trade-off, if the deleted exploratory script is restored: overall MCC versus hub-edge and non-hub-edge recovery.

## Goal 3: Robustness, Network Structure, and Method Benchmarks

Folder: `sim_script/03_robustness_benchmarks/`

This goal evaluates how PSS-Net behaves under more realistic constraints and structural heterogeneity. It should also become the home for external method comparisons against black-box or association-based network inference methods.

Included scripts:

- `pss_net_compositional.R`: quantifies the failure modes caused by relative abundance, CLR transformation, and noisy total-abundance reconstruction.
- `pss_net_scalefree.R`: tests node-wise versus joint block-diagonal estimation under hub-like scale-free structure and evaluates hub recovery.
- `pss_net_joint_smalln.R`: provides a homogeneous-network comparison for interpreting when joint estimation helps or does not help.

Planned external comparisons:

- GENIE3 or GRNBoost2 as tree-based black-box directed network baselines using expression/state data only.
- Correlation, partial correlation, graphical lasso, SPIEC-EASI, or SparCC as association/compositional baselines where appropriate.
- Lasso, elastic net, or STLSQ/SINDy-style library regression as sparse regression baselines using the same steady-state equation.
- Existing group lasso results remain an internal structured-sparsity baseline.

Main figure concept: Figure 3, robustness and benchmarks.

Suggested panels:

- Panel 3a: Node-wise versus joint MCC in scale-free networks.
- Panel 3b: Hub recovery in scale-free networks, using estimated out-degree Spearman correlation and top-k hit rate.
- Panel 3c: Homogeneous or uniform-degree negative control showing that joint estimation is structure-dependent.
- Panel 3d: Absolute abundance versus relative abundance, CLR, and noisy total-abundance reconstruction.
- Panel 3e: PSS-Net versus external network inference methods on matched simulated PSS data.
- Panel 3f: Method capability matrix: uses perturbation input, outputs direction, estimates nonlinear edge functions, requires time series, handles compositional data.

## Current Priority

The immediate publication gap is external method comparison. The existing simulations already cover internal validation, scaling, design, and limitations, but the final manuscript should include fair benchmarks against methods that do not use perturbation input, methods that use sparse regression without double sparsity, and compositional/association methods for microbiome-style data.

# Supplementary Simulation Plan

This directory records planned supplementary simulations for PSS-Net. These studies are important for reviewer confidence but should not dominate the main simulation figures unless they overturn a main claim.

## S1: Model Misspecification

Questions:

- What happens when some true regulators are hidden or unobserved?
- How sensitive is recovery to off-target perturbations or noisy perturbation magnitudes?
- How badly does recovery degrade when the basis dictionary is wrong or over-complete?
- Can PSS-Net detect failure under multiple steady states or non-convergent conditions?

Planned outputs:

- MCC/AUPR versus hidden-node fraction.
- MCC/AUPR versus perturbation-input noise level.
- Edge-sign and Jacobian RMSE under wrong basis dictionaries.
- Failure-rate table for invalid steady states or non-convergence.

## S2: Perturbation Realism

Questions:

- How much correlated perturbation can be tolerated?
- What fraction of nodes must be perturbable for whole-network recovery?
- How does targeted perturbation trade hub recovery against global network recovery?

Planned outputs:

- Recovery curves versus perturbable-node fraction.
- Hub-edge recall and non-hub-edge recall under targeted perturbation.
- Diagnostic warning when a target node has nearly constant response `u_j`.

## S3: Measurement Noise and Compositional Sensitivity

Questions:

- How do additive, multiplicative, count-like, and compositional noise affect recovery?
- How sensitive is absolute-abundance reconstruction to total-biomass measurement error?

Planned outputs:

- MCC/AUPR versus measurement noise.
- Total-biomass CV sweep for relative abundance times noisy total abundance.
- Comparison of absolute, relative, CLR, and total-corrected inputs.

## S4: Monte Carlo Stability

Questions:

- Are main-figure differences stable under larger seed counts?
- How often do individual methods fail or return degenerate networks?

Planned outputs:

- Paired method differences with confidence intervals.
- Failure-rate accounting and NA-handling rules.
- Sensitivity to thresholding versus rank-based AUPR/AUROC.

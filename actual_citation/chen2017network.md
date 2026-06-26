# chen2017network

- **Title:** Network Reconstruction From High-Dimensional Ordinary Differential Equations
- **Authors:** Shizhe Chen, Ali Shojaie, Daniela M. Witten
- **Journal/Year:** Journal of the American Statistical Association, 112(520), 1697–1707, 2017
- **DOI:** 10.1080/01621459.2016.1229197
- **Source consulted:** Europe PMC REST (DOI:10.1080/01621459.2016.1229197)

## Verbatim abstract

"We consider the task of learning a dynamical system from high-dimensional time-course
data. For instance, we might wish to estimate a gene regulatory network from gene
expression data measured at discrete time points. We model the dynamical system
nonparametrically as a system of additive ordinary differential equations. Most existing
methods for parameter estimation in ordinary differential equations estimate the
derivatives from noisy observations. This is known to be challenging and inefficient. We
propose a novel approach that does not involve derivative estimation. We show that the
proposed method can consistently recover the true network structure even in high
dimensions, and we demonstrate empirical improvement over competing approaches."

## Supports in intro

B2/B4 — a rigorous high-dimensional additive-ODE network-reconstruction method with
*consistency guarantees*, and one that explicitly "does not involve derivative
estimation" (it matches integrated trajectories rather than estimating $dx/dt$).
Strengthens the additive-ODE lineage cited alongside \citep{henderson2014network,
wu2014saode}, and is thematically aligned with PSS-Net's avoidance of numerical
differentiation. Key distinction: Chen et al. still work from time-course trajectories
(integral matching), whereas PSS-Net uses perturbed *steady states*, turning the problem
into an algebraic one with no trajectory at all.

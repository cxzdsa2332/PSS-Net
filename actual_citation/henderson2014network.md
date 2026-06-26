# henderson2014network

- **Title:** Network Reconstruction Using Nonparametric Additive ODE Models
- **Authors:** James Henderson, George Michailidis
- **Journal/Year:** PLoS ONE, 9(4), e94003, 2014
- **DOI:** 10.1371/journal.pone.0094003
- **Source consulted:** https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0094003

## Verbatim abstract

"Network representations of biological systems are widespread and reconstructing
unknown networks from data is a focal problem for computational biologists. ...
In this paper, we introduce an approach to reconstructing directed networks based on
dynamic systems models. Our approach generalizes commonly used ODE models based on
linear or nonlinear dynamics by extending the functional class for the functions
involved from parametric to nonparametric models. Concomitantly we limit the
complexity by imposing an additive structure on the estimated slope functions. Thus
the submodel associated with each node is a sum of univariate functions. These
univariate component functions form the basis for a novel coupling metric that we
define in order to quantify the strength of proposed relationships and hence rank
potential edges. ... We compare our method to those that similarly rely on dynamic
systems models and use the results to attempt to disentangle the distinct roles of
linearity, sparsity, and derivative estimation."

## Supports in intro

B2/B4 — the additive-nonparametric ODE prior art: each node's submodel is "a sum of
univariate functions," and edges are ranked by univariate coupling strength. PSS-Net
adopts the same additive decomposition but estimates it at perturbed steady state
(algebraic) rather than from differentiated trajectories, and selects edges by
double-sparse regression instead of a coupling-metric ranking. Their explicit framing
of "the distinct roles of linearity, sparsity, and derivative estimation" is the gap
PSS-Net addresses on the derivative-estimation axis.

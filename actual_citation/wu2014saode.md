# wu2014saode

- **Title:** Sparse Additive Ordinary Differential Equations for Dynamic Gene Regulatory Network Modeling
- **Authors:** Hulin Wu, Tao Lu, Hongqi Xue, Hua Liang
- **Journal/Year:** Journal of the American Statistical Association, 109(506), 700–716, 2014
- **DOI:** 10.1080/01621459.2013.859617
- **Source consulted:** Europe PMC REST (DOI:10.1080/01621459.2013.859617)

## Verbatim abstract

"The gene regulation network (GRN) is a high-dimensional complex system, which can be
represented by various mathematical or statistical models. The ordinary differential
equation (ODE) model is one of the popular dynamic GRN models. High-dimensional linear
ODE models have been proposed to identify GRNs, but with a limitation of the linear
regulation effect assumption. In this article, we propose a sparse additive ODE
(SA-ODE) model, coupled with ODE estimation methods and adaptive group LASSO
techniques, to model dynamic GRNs that could flexibly deal with nonlinear regulation
effects. The asymptotic properties of the proposed method are established and
simulation studies are performed to validate the proposed approach. An application
example for identifying the nonlinear dynamic GRN of T-cell activation is used to
illustrate the usefulness of the proposed method."

## Supports in intro

B2/B4 — the SA-ODE precedent for combining additive nonparametric regulation with
*group* sparsity (adaptive group LASSO). It explicitly motivates leaving the linear
ODE class because of "the limitation of the linear regulation effect assumption."
PSS-Net keeps the sparse-additive structure but (i) works at perturbed steady state
rather than from time-course derivatives and (ii) uses double sparsity (group +
within-group) instead of group LASSO alone.

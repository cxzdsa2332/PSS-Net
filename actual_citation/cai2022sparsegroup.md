# cai2022sparsegroup

- **Title:** Sparse Group Lasso: Optimal Sample Complexity, Convergence Rate, and Statistical Inference
- **Authors:** T. Tony Cai, Anru R. Zhang, Yuchen Zhou
- **Journal/Year:** IEEE Transactions on Information Theory, 68(9), 5975–6002, 2022
- **arXiv:** 1909.09851
- **Source consulted:** https://arxiv.org/abs/1909.09851

## Verbatim abstract

"We study sparse group Lasso for high-dimensional double sparse linear regression,
where the parameter of interest is simultaneously element-wise and group-wise sparse.
This problem is an important instance of the simultaneously structured model -- an
actively studied topic in statistics and machine learning. In the noiseless case,
matching upper and lower bounds on sample complexity are established for the exact
recovery of sparse vectors and for stable estimation of approximately sparse vectors,
respectively. In the noisy case, upper and matching minimax lower bounds for estimation
error are obtained. We also consider the debiased sparse group Lasso and investigate
its asymptotic property for the purpose of statistical inference. Finally, numerical
studies are provided to support the theoretical results."

## Supports in intro

B4 — formal statistical grounding of the "double sparse" (simultaneously element-wise
and group-wise sparse) regression problem that PSS-Net's node-wise inference instantiates;
provides the sample-complexity rationale for exploiting both group and within-group
structure.

# zhang2023minimax

- **Title:** A minimax optimal approach to high-dimensional double sparse linear regression
- **Authors:** Yanhang Zhang, Zhifan Li, Shixiang Liu, Jianxin Yin
- **Journal/Year:** Journal of Machine Learning Research, 25, 2024
- **arXiv:** 2305.04182 (ADSIHT; CRAN package)
- **Source consulted:** https://arxiv.org/abs/2305.04182

## Verbatim abstract

"In this paper, we focus our attention on the high-dimensional double sparse linear
regression, that is, a combination of element-wise and group-wise sparsity. To address
this problem, we propose an IHT-style (iterative hard thresholding) procedure that
dynamically updates the threshold at each step. We establish the matching upper and
lower bounds for parameter estimation, showing the optimality of our proposal in the
minimax sense. More importantly, we introduce a fully adaptive optimal procedure
designed to address unknown sparsity and noise levels. Our adaptive procedure
demonstrates optimal statistical accuracy with fast convergence. Additionally, we
elucidate the significance of the element-wise sparsity level s_0 as the trade-off
between IHT and group IHT, underscoring the superior performance of our method over
both. Leveraging the beta-min condition, we establish that our IHT-style procedure can
attain the oracle estimation rate and achieve almost full recovery of the true support
set at both the element level and group level. Finally, we demonstrate the superiority
of our method by comparing it with several state-of-the-art algorithms on both
synthetic and real-world datasets."

## Supports in intro

B4 — the estimator PSS-Net uses (ADSIHT). Double sparsity = "a combination of
element-wise and group-wise sparsity," which maps onto PSS-Net's structure: group =
source node (select which sources couple to a target), element = basis coefficient
within a source's univariate function. Provides the minimax-optimal, adaptive
selection that motivates choosing it over plain group lasso.

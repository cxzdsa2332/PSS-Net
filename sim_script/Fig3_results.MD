# Fig3 raw results

Generated from current result CSV files. Tables report group means and SDs only.

## Data files

| result | file | rows | repeats |
| --- | --- | --- | --- |
| Fig3b | results/sim_results/Fig3b_external_benchmark_main.csv | 1920 | 30 |
| Fig3d | results/sim_results/Fig3c_structure_dependence.csv | 240 | 30 |
| Fig3e | results/sim_results/Fig3c_compositional_data_limitation.csv | 120 | 30 |

## Fig3b: external benchmark at representative budget

Representative budget: `N = 171`, `N_over_slogp = 25.138`.

| truth | N | N_over_slogp | method | n_seed | MCC_mean | MCC_sd | AUPRC_mean | AUPRC_sd | AUROC_mean | AUROC_sd | Precision_mean | Precision_sd | Recall_mean | Recall_sd | FuncNRMSE_mean | FuncNRMSE_sd | EdgeJacNRMSE_mean | EdgeJacNRMSE_sd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| linear | 171 | 25.138 | aiMeRA | 30 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | NA | NA | 0.023 | 0.003 |
| linear | 171 | 25.138 | Correlation | 30 | 0.165 | 0.084 | 0.289 | 0.046 | 0.775 | 0.033 | 0.507 | 0.234 | 0.071 | 0.049 | NA | NA | NA | NA |
| linear | 171 | 25.138 | ElasticNet | 30 | 0.312 | 0.021 | 1.000 | 0.000 | 1.000 | 0.000 | 0.160 | 0.012 | 1.000 | 0.000 | 0.063 | 0.005 | 0.043 | 0.004 |
| linear | 171 | 25.138 | GroupLasso | 30 | 0.232 | 0.010 | 1.000 | 0.000 | 1.000 | 0.000 | 0.119 | 0.004 | 1.000 | 0.000 | 0.103 | 0.005 | 0.036 | 0.003 |
| linear | 171 | 25.138 | Lasso | 30 | 0.840 | 0.046 | 1.000 | 0.000 | 1.000 | 0.000 | 0.728 | 0.070 | 1.000 | 0.000 | 0.084 | 0.006 | 0.073 | 0.004 |
| linear | 171 | 25.138 | PartialCor | 30 | 0.095 | 0.091 | 0.258 | 0.039 | 0.741 | 0.036 | 0.372 | 0.345 | 0.033 | 0.041 | NA | NA | NA | NA |
| linear | 171 | 25.138 | PSS_Net | 30 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 1.000 | 0.000 | 0.023 | 0.005 | 0.021 | 0.003 |
| linear | 171 | 25.138 | PySINDy_STLSQ | 30 | 0.173 | 0.094 | 0.120 | 0.042 | 0.509 | 0.047 | 0.867 | 0.346 | 0.041 | 0.034 | 0.960 | 0.033 | 0.963 | 0.030 |
| strong_nonlinear | 171 | 25.138 | aiMeRA | 30 | 0.964 | 0.016 | 0.987 | 0.014 | 0.995 | 0.007 | 1.000 | 0.000 | 0.934 | 0.030 | NA | NA | 0.038 | 0.008 |
| strong_nonlinear | 171 | 25.138 | Correlation | 30 | 0.262 | 0.074 | 0.294 | 0.047 | 0.746 | 0.045 | 0.511 | 0.070 | 0.174 | 0.079 | NA | NA | NA | NA |
| strong_nonlinear | 171 | 25.138 | ElasticNet | 30 | 0.302 | 0.057 | 0.996 | 0.009 | 0.998 | 0.005 | 0.157 | 0.042 | 0.998 | 0.007 | 0.128 | 0.034 | 0.042 | 0.018 |
| strong_nonlinear | 171 | 25.138 | GroupLasso | 30 | 0.250 | 0.127 | 0.995 | 0.016 | 0.998 | 0.006 | 0.145 | 0.155 | 0.995 | 0.022 | 0.154 | 0.071 | 0.040 | 0.055 |
| strong_nonlinear | 171 | 25.138 | Lasso | 30 | 0.629 | 0.089 | 0.991 | 0.012 | 0.995 | 0.007 | 0.450 | 0.109 | 0.990 | 0.015 | 0.151 | 0.025 | 0.058 | 0.009 |
| strong_nonlinear | 171 | 25.138 | PartialCor | 30 | 0.196 | 0.096 | 0.252 | 0.055 | 0.711 | 0.052 | 0.491 | 0.190 | 0.103 | 0.070 | NA | NA | NA | NA |
| strong_nonlinear | 171 | 25.138 | PSS_Net | 30 | 0.965 | 0.017 | 0.990 | 0.010 | 0.996 | 0.006 | 1.000 | 0.000 | 0.936 | 0.030 | 0.066 | 0.017 | 0.022 | 0.005 |
| strong_nonlinear | 171 | 25.138 | PySINDy_STLSQ | 30 | 0.454 | 0.044 | 0.307 | 0.040 | 0.596 | 0.040 | 1.000 | 0.000 | 0.219 | 0.041 | 0.803 | 0.037 | 0.610 | 0.066 |

## Fig3b/Fig3c diagnostic: all-budget edge-strength and edge-function recovery

| truth | method | n_seed | EdgeJacNRMSE_mean | EdgeJacNRMSE_sd | FuncNRMSE_mean | FuncNRMSE_sd | EdgeLinearJacRMSE_mean | EdgeLinearJacRMSE_sd | EdgeNonlinearJacRMSE_mean | EdgeNonlinearJacRMSE_sd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| linear | ElasticNet | 30 | 0.072 | 0.048 | 0.097 | 0.056 | 0.059 | 0.023 | 0.046 | 0.014 |
| linear | GroupLasso | 30 | 0.048 | 0.020 | 0.114 | 0.018 | 0.090 | 0.010 | 0.082 | 0.007 |
| linear | Lasso | 30 | 0.087 | 0.025 | 0.100 | 0.030 | 0.042 | 0.016 | 0.029 | 0.013 |
| linear | PSS_Net | 30 | 0.025 | 0.008 | 0.031 | 0.014 | 0.036 | 0.024 | 0.034 | 0.025 |
| linear | PySINDy_STLSQ | 30 | 0.969 | 0.025 | 0.966 | 0.028 | 0.220 | 0.009 | 0.000 | 0.000 |
| strong_nonlinear | ElasticNet | 30 | 0.071 | 0.046 | 0.166 | 0.068 | 0.124 | 0.028 | 0.129 | 0.033 |
| strong_nonlinear | GroupLasso | 30 | 0.115 | 0.141 | 0.248 | 0.181 | 0.152 | 0.031 | 0.158 | 0.042 |
| strong_nonlinear | Lasso | 30 | 0.071 | 0.023 | 0.169 | 0.043 | 0.118 | 0.019 | 0.121 | 0.021 |
| strong_nonlinear | PSS_Net | 30 | 0.030 | 0.019 | 0.079 | 0.031 | 0.122 | 0.036 | 0.118 | 0.035 |
| strong_nonlinear | PySINDy_STLSQ | 30 | 0.621 | 0.067 | 0.810 | 0.040 | 0.231 | 0.013 | 0.227 | 0.022 |

## Fig3d: structure dependence

| topology | N | n_seed | MCC_nodewise_mean | MCC_nodewise_sd | MCC_joint_mean | MCC_joint_sd | rho_nodewise_mean | rho_nodewise_sd | rho_joint_mean | rho_joint_sd | topk_nodewise_mean | topk_nodewise_sd | topk_joint_mean | topk_joint_sd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| er | 16 | 30 | 0.094 | 0.030 | 0.112 | 0.037 | 0.113 | 0.147 | 0.084 | 0.155 | 0.160 | 0.122 | 0.133 | 0.169 |
| er | 24 | 30 | 0.254 | 0.045 | 0.305 | 0.054 | 0.240 | 0.132 | 0.331 | 0.129 | 0.213 | 0.157 | 0.247 | 0.187 |
| er | 32 | 30 | 0.441 | 0.059 | 0.440 | 0.090 | 0.391 | 0.141 | 0.356 | 0.130 | 0.287 | 0.201 | 0.320 | 0.201 |
| er | 40 | 30 | 0.590 | 0.046 | 0.603 | 0.064 | 0.548 | 0.114 | 0.562 | 0.135 | 0.460 | 0.190 | 0.393 | 0.238 |
| scalefree | 16 | 30 | 0.080 | 0.036 | 0.110 | 0.058 | 0.144 | 0.138 | 0.188 | 0.182 | 0.247 | 0.136 | 0.273 | 0.162 |
| scalefree | 24 | 30 | 0.202 | 0.082 | 0.259 | 0.119 | 0.295 | 0.160 | 0.367 | 0.181 | 0.467 | 0.206 | 0.520 | 0.179 |
| scalefree | 32 | 30 | 0.355 | 0.107 | 0.408 | 0.130 | 0.442 | 0.124 | 0.557 | 0.135 | 0.680 | 0.220 | 0.660 | 0.198 |
| scalefree | 40 | 30 | 0.484 | 0.138 | 0.520 | 0.141 | 0.529 | 0.191 | 0.642 | 0.175 | 0.700 | 0.188 | 0.687 | 0.194 |

## Fig3e: compositional data limitation

| input | n_seed | MCC_mean | MCC_sd | Pr_mean | Pr_sd | Re_mean | Re_sd | n_eff_mean | n_eff_sd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| abs | 30 | 0.890 | 0.052 | 0.925 | 0.041 | 0.929 | 0.054 | 199.900 | 0.403 |
| clr | 30 | 0.317 | 0.129 | 0.533 | 0.089 | 0.584 | 0.089 | 199.900 | 0.403 |
| rel | 30 | 0.520 | 0.124 | 0.739 | 0.098 | 0.596 | 0.102 | 199.900 | 0.403 |
| rel_x_T | 30 | -0.047 | 0.134 | 0.318 | 0.049 | 0.581 | 0.082 | 199.900 | 0.403 |

# xiao2017mapping

- **Title:** Mapping the ecological networks of microbial communities
- **Authors:** Yandong Xiao, Marco Tulio Angulo, Jonathan Friedman, Matthew K. Waldor, Scott T. Weiss, Yang-Yu Liu
- **Journal/Year:** Nature Communications, 8(1), 2042, 2017
- **DOI:** 10.1038/s41467-017-02090-2
- **Source consulted:** Europe PMC REST (DOI:10.1038/s41467-017-02090-2)

## Verbatim abstract

"Mapping the ecological networks of microbial communities is a necessary step toward
understanding their assembly rules and predicting their temporal behavior. However,
existing methods require assuming a particular population dynamics model, which is not
known a priori. Moreover, those methods require fitting longitudinal abundance data,
which are often not informative enough for reliable inference. To overcome these
limitations, here we develop a new method based on steady-state abundance data. Our
method can infer the network topology and inter-taxa interaction types without assuming
any particular population dynamics model. Additionally, when the population dynamics is
assumed to follow the classic Generalized Lotka-Volterra model, our method can infer the
inter-taxa interaction strengths and intrinsic growth rates. We systematically validate
our method using simulated data, and then apply it to four experimental data sets. Our
method represents a key step towards reliable modeling of complex, real-world microbial
communities, such as the human gut microbiota."

## Supports in intro

B3 — the closest microbial-ecology prior art for the equilibrium paradigm: infers
network topology and interaction *signs/types* from steady-state abundance data, with no
dynamics-model assumption, and recovers interaction *strengths* only once a fixed gLV
form is imposed. Directly supports the sentence that steady-state data alone can reveal
topology/signs, and motivates PSS-Net's step to full nonlinear interaction functions
rather than signs or fixed-form strengths. (Per CLAUDE.md, the closest related steady-
state inference work; Xiao et al. use steady-state differences for sign/type.)

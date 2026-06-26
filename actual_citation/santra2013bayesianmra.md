# santra2013bayesianmra

- **Title:** Integrating Bayesian variable selection with Modular Response Analysis to infer biochemical network topology
- **Authors:** Tapesh Santra, Walter Kolch, Boris N. Kholodenko
- **Journal/Year:** BMC Systems Biology, 7, 57, 2013
- **DOI:** 10.1186/1752-0509-7-57
- **Source consulted:** Europe PMC REST (DOI:10.1186/1752-0509-7-57)

## Verbatim abstract

"Background: Recent advancements in genetics and proteomics have led to the acquisition
of large quantitative data sets. However, the use of these data to reverse engineer
biochemical networks has remained a challenging problem. Many methods have been proposed
to infer biochemical network topologies from different types of biological data. Here,
we focus on unraveling network topologies from steady state responses of biochemical
networks to successive experimental perturbations.
Results: We propose a computational algorithm which combines a deterministic network
inference method termed Modular Response Analysis (MRA) and a statistical model
selection algorithm called Bayesian Variable Selection, to infer functional interactions
in cellular signaling pathways and gene regulatory networks. It can be used to identify
interactions among individual molecules involved in a biochemical pathway or reveal how
different functional modules of a biological network interact with each other to exchange
information. In cases where not all network components are known, our method reveals
functional interactions which are not direct but correspond to the interaction routes
through unknown elements. Using computer simulated perturbation responses of signaling
pathways and gene regulatory networks from the DREAM challenge, we demonstrate that the
proposed method is robust against noise and scalable to large networks. We also show that
our method can infer network topologies using incomplete perturbation datasets.
Consequently, we have used this algorithm to explore the ERBB regulated G1/S transition
pathway in certain breast cancer cells to understand the molecular mechanisms which cause
these cells to become drug resistant. The algorithm successfully inferred many well
characterized interactions of this pathway by analyzing experimentally obtained
perturbation data. Additionally, it identified some molecular interactions which promote
drug resistance in breast cancer cells.
Conclusions: The proposed algorithm provides a robust, scalable and cost effective
solution for inferring network topologies from biological data. It can potentially be
applied to explore novel pathways which play important roles in life threatening disease
like cancer."

## Supports in intro

B3 — extends the MRA lineage (with Kholodenko as a co-author) by coupling MRA with
*Bayesian variable selection* on steady-state perturbation responses, i.e. an explicit
statistical model-selection step over edges. Strengthens the MRA citation and
foreshadows PSS-Net's use of sparse selection, while remaining within the
local/linearised-response framework that PSS-Net moves beyond.

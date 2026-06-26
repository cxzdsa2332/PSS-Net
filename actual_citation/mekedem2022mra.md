# mekedem2022mra

- **Title:** Application of modular response analysis to medium- to large-size biological systems
- **Authors:** Meriem Mekedem, Patrice Ravel, Jacques Colinge
- **Journal/Year:** PLOS Computational Biology, 18(4), e1009312, 2022
- **DOI:** 10.1371/journal.pcbi.1009312
- **Source consulted:** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009312

## Verbatim abstract

"The development of high-throughput genomic technologies associated with recent genetic
perturbation techniques such as short hairpin RNA (shRNA), gene trapping, or gene
editing (CRISPR/Cas9) has made it possible to obtain large perturbation data sets.
These data sets are invaluable sources of information regarding the function of genes,
and they offer unique opportunities to reverse engineer gene regulatory networks in
specific cell types. Modular response analysis (MRA) is a well-accepted mathematical
modeling method that is precisely aimed at such network inference tasks, but its use
has been limited to rather small biological systems so far. In this study, we show that
MRA can be employed on large systems with almost 1,000 network components. In
particular, we show that MRA performance surpasses general-purpose mutual
information-based algorithms. Part of these competitive results was obtained by the
application of a novel heuristic that pruned MRA-inferred interactions a posteriori. We
also exploited a block structure in MRA linear algebra to parallelize large system
resolutions."

## Supports in intro

B3 — shows the perturbation/steady-state route scaling to ~1,000 components and
outperforming mutual-information methods, but still within the *linear* MRA framework.
Motivates PSS-Net's step to nonparametric nonlinear interaction functions while keeping
the perturbation-response philosophy.

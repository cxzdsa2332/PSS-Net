# trapnell2014pseudotime

- **Title:** The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells
- **Authors:** Cole Trapnell, Davide Cacchiarelli, Jonna Grimsby, Prapti Pokharel, Shuqiang Li, Michael Morse, Niall J. Lennon, Kenneth J. Livak, Tarjei S. Mikkelsen, John L. Rinn
- **Journal/Year:** Nature Biotechnology, 32(4), 381–386, 2014
- **DOI:** 10.1038/nbt.2859
- **Source consulted:** Europe PMC REST (DOI:10.1038/nbt.2859)

## Verbatim abstract

"Defining the transcriptional dynamics of a temporal process such as cell
differentiation is challenging owing to the high variability in gene expression between
individual cells. Time-series gene expression analyses of bulk cells have difficulty
distinguishing early and late phases of a transcriptional cascade or identifying rare
subpopulations of cells, and single-cell proteomic methods rely on a priori knowledge of
key distinguishing markers. Here we describe Monocle, an unsupervised algorithm that
increases the temporal resolution of transcriptome dynamics using single-cell RNA-Seq
data collected at multiple time points. Applied to the differentiation of primary human
myoblasts, Monocle revealed switch-like changes in expression of key regulatory factors,
sequential waves of gene regulation, and expression of regulators that were not known to
act in differentiation. We validated some of these predicted regulators in a
loss-of-function screen. Monocle can in principle be used to recover single-cell gene
expression kinetics from a wide array of cellular processes, including differentiation,
proliferation and oncogenic transformation."

## Supports in intro

B3 — canonical pseudotemporal-ordering method (Monocle): orders single cells along an
inferred trajectory to reconstruct dynamics from data lacking dense/explicit time
stamps. Cited for the "pseudo-time" route that converts snapshot data into a pseudo
time-course before dynamics/network inference.

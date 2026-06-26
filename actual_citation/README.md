# actual_citation/

Verbatim source excerpts for every reference cited in `manuscript/sections/01_introduction.tex`.

Purpose: provenance check. Each file records (i) the bib key and metadata, (ii) the
public source URL consulted, (iii) the **verbatim abstract / key sentence** from the
original paper, and (iv) the specific claim in the Introduction that the citation
supports. This lets any co-author confirm that every `\citep` is grounded in what the
source actually says, not a paraphrase from memory.

Provenance note: abstracts were retrieved 2026-06-26 from open sources — arXiv,
PLoS, Europe PMC REST API (`ebi.ac.uk/europepmc`), Crossref, and Semantic Scholar.
Where a retrieval tool condensed the abstract, the verbatim opening sentence(s) are
quoted exactly and any condensed remainder is explicitly marked `[summary]`.

| key | claim supported (intro block) |
|-----|-------------------------------|
| barzel2013universality | B1/B2 — universality of network dynamics motivates a model-agnostic functional form |
| stein2013ecological | B1/B2 — gLV time-series inference of gut microbiota |
| venturelli2018deciphering | B1 — pairwise interactions drive synthetic gut community dynamics |
| dong2023idopnetwork | B1 — network reconstruction from community abundance data |
| bucci2016mdsine | B2/B3 — dynamical-systems inference from microbiome time series |
| henderson2014network | B2/B4 — nonparametric additive ODE network reconstruction |
| wu2014saode | B2/B4 — sparse additive ODE, nonlinear regulation, group lasso |
| brunton2016sindy | B3 — sparse identification of nonlinear dynamics from data |
| kaptanoglu2022pysindy | B3 — open-source SINDy implementation (software ref; bib-verified) |
| sindy2023microbiota | B3 — SINDy applied to microbiota pairwise interactions |
| kholodenko2002untangling | B3 — steady-state perturbation / modular response analysis |
| mekedem2022mra | B3 — MRA scaled to large perturbation systems |
| meister2013pss | B3/B4 — perturbed steady-state reformulation removes time derivatives |
| zhang2023minimax | B4 — minimax-optimal double-sparse estimator (ADSIHT) |
| cai2022sparsegroup | B4 — sparse-group (double) sparsity sample complexity |

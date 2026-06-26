# dong2026multitask

- **Title:** Multi-task learning of complex networks via nonlinear ordinary differential equations
- **Authors:** Ang Dong, Changjian Fa, Zhifan Li, Shing-Tung Yau, Rongling Wu
- **Journal/Year:** Communications Physics, 2026
- **DOI:** 10.1038/s42005-026-02687-4
- **Source consulted:** OpenAlex (reconstructed abstract, DOI:10.1038/s42005-026-02687-4)

## Abstract (OpenAlex reconstruction from abstract_inverted_index — near-verbatim)

"Complex systems are characterized by many underlying entities and their intricate
interactions. We contextualize ecological niche theory through evolutionary game theory
and a system of nonlinear mixed ordinary differential equations (nMODEs) to reconstruct
informative, dynamic, omnidirectional, and personalized networks (idopNetworks) for
complex systems at any dimension. We implement a multi-task learning (MTL) algorithm
into the matrix form of linearized nMODEs to execute two coupled tasks for group-level
and elementwise sparsity on nonlinear feature representations. Beyond existing
networking practice, MTL-based idopNetworks can capture all-around interacting links,
nonlinearities, and emergent properties of a system, which, to a larger extent,
approximate the complexity of complex systems. We apply our model to learn gene
regulatory networks from transcriptional data for parasite Plasmodium falciparum,
identifying previously-unknown regulatory roles of several genes in mediating malaria
infection. Our model provides insight of machine learning to analyze, model, and
interpret complex systems in non-Euclidean space."

(Reconstructed from OpenAlex inverted index; word order faithful to the source index.)

## Supports in intro

B2/B4 — directly adjacent prior art from the same group: a nonlinear-ODE network model
that performs "group-level and elementwise sparsity on nonlinear feature
representations," i.e. the same double-sparse, nonlinear, additive philosophy as
PSS-Net, but fit jointly via a multi-task / block-diagonal linearized-nMODE solve on
(static/transcriptional) data. PSS-Net differs by using *perturbed steady states*
(algebraic, no derivatives) and a node-wise double-sparse solve. See CLAUDE.md for the
joint-vs-node-wise grouping caveat tied to this reference.

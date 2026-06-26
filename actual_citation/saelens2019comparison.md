# saelens2019comparison

- **Title:** A comparison of single-cell trajectory inference methods
- **Authors:** Wouter Saelens, Robrecht Cannoodt, Helena Todorov, Yvan Saeys
- **Journal/Year:** Nature Biotechnology, 37(5), 547–554, 2019
- **DOI:** 10.1038/s41587-019-0071-9
- **Source consulted:** Europe PMC REST (DOI:10.1038/s41587-019-0071-9)

## Verbatim abstract

"Trajectory inference approaches analyze genome-wide omics data from thousands of single
cells and computationally infer the order of these cells along developmental
trajectories. Although more than 70 trajectory inference tools have already been
developed, it is challenging to compare their performance because the input they require
and output models they produce vary substantially. Here, we benchmark 45 of these
methods on 110 real and 229 synthetic datasets for cellular ordering, topology,
scalability and usability. Our results highlight the complementarity of existing tools,
and that the choice of method should depend mostly on the dataset dimensions and
trajectory topology. Based on these results, we develop a set of guidelines to help users
select the best method for their dataset. Our freely available data and evaluation
pipeline ( https://benchmark.dynverse.org ) will aid in the development of improved tools
designed to analyze increasingly large and complex single-cell datasets."

## Supports in intro

B3 — benchmark of 45 trajectory-inference methods showing that "the input they require
and output models they produce vary substantially" and that the best choice depends on
dataset/topology. Supports the critique that the pseudo-time *ordering* is itself
uncertain and strongly method-dependent, so downstream ODE/network inference inherits
this added uncertainty on top of the derivative-estimation problem.

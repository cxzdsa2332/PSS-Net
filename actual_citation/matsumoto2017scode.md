# matsumoto2017scode

- **Title:** SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation
- **Authors:** Hirotaka Matsumoto, Hisanori Kiryu, Chikara Furusawa, Minoru S. H. Ko, Shigeru B. H. Ko, Norio Gouda, Tetsutaro Hayashi, Itoshi Nikaido
- **Journal/Year:** Bioinformatics, 33(15), 2314–2321, 2017
- **DOI:** 10.1093/bioinformatics/btx194
- **Source consulted:** Europe PMC REST (DOI:10.1093/bioinformatics/btx194)

## Verbatim abstract

"The analysis of RNA-Seq data from individual differentiating cells enables us to
reconstruct the differentiation process and the degree of differentiation (in
pseudo-time) of each cell. Such analyses can reveal detailed expression dynamics and
functional relationships for differentiation. To further elucidate differentiation
processes, more insight into gene regulatory networks is required. The pseudo-time can be
regarded as time information and, therefore, single-cell RNA-Seq data are time-course
data with high time resolution. Although time-course data are useful for inferring
networks, conventional inference algorithms for such data suffer from high time
complexity when the number of samples and genes is large. Therefore, a novel algorithm is
necessary to infer networks from single-cell RNA-Seq during differentiation. In this
study, we developed the novel and efficient algorithm SCODE to infer regulatory networks,
based on ordinary differential equations. We applied SCODE to three single-cell RNA-Seq
datasets and confirmed that SCODE can reconstruct observed expression dynamics. We
evaluated SCODE by comparing its inferred networks with use of a DNaseI-footprint based
network. The performance of SCODE was best for two of the datasets and nearly best for
the remaining dataset. We also compared the runtimes and showed that the runtimes for
SCODE are significantly shorter than for alternatives. Thus, our algorithm provides a
promising approach for further single-cell differentiation analyses."

## Supports in intro

B3 — the key reference tying pseudotime to ODE-based network inference: it explicitly
states "the pseudo-time can be regarded as time information and, therefore, single-cell
RNA-Seq data are time-course data," then infers gene regulatory networks "based on
ordinary differential equations." Directly supports the clause that, after pseudo-time
ordering, networks are inferred with ODE methods treating pseudotime as time — and
illustrates that such methods still inherit the time-course derivative/ordering issues.

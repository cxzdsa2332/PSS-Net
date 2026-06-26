# stein2013ecological

- **Title:** Ecological Modeling from Time-Series Inference: Insight into Dynamics and Stability of Intestinal Microbiota
- **Authors:** Richard R. Stein, Vanni Bucci, Nora C. Toussaint, Charlie G. Buffie, Gunnar Rätsch, Eric G. Pamer, Chris Sander, João B. Xavier
- **Journal/Year:** PLoS Computational Biology, 9(12), e1003388, 2013
- **DOI:** 10.1371/journal.pcbi.1003388
- **Source consulted:** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003388

## Verbatim abstract

"The intestinal microbiota is a microbial ecosystem of crucial importance to human
health. Understanding how the microbiota confers resistance against enteric pathogens
and how antibiotics disrupt that resistance is key to the prevention and cure of
intestinal infections. We present a novel method to infer microbial community ecology
directly from time-resolved metagenomics. This method extends generalized
Lotka–Volterra dynamics to account for external perturbations. Data from recent
experiments on antibiotic-mediated Clostridium difficile infection is analyzed to
quantify microbial interactions, commensal-pathogen interactions, and the effect of the
antibiotic on the community. Stability analysis reveals that the microbiota is
intrinsically stable, explaining how antibiotic perturbations and C. difficile
inoculation can produce catastrophic shifts that persist even after removal of the
perturbations. Importantly, the analysis suggests a subnetwork of bacterial groups
implicated in protection against C. difficile. Due to its generality, our method can be
applied to any high-resolution ecological time-series data to infer community structure
and response to external stimuli."

## Supports in intro

B1/B2 — the canonical gLV-from-time-series approach for the gut microbiota. Note the
requirement of "high-resolution ecological time-series data," which PSS-Net replaces
with perturbed steady-state measurements. Also a source for the gLV motivation and for
the linear-interaction assumption that PSS-Net relaxes.

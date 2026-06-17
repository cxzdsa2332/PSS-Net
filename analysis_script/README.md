# Analysis scripts

This directory contains scripts that turn numeric simulation outputs into
figures, tables, and compact summaries. Formal analysis scripts should read from
`results/sim_results/` and write to `results/figure/` or `results/table/`.

Short-term convention:

- Use `Fig1a_`, `Fig1b_`, ... prefixes for manuscript-facing panel scripts.
- Keep generic helper or legacy summary scripts here until enough panel scripts
  exist to justify figure-level subfolders.
- Use `assemble_Fig1.R`, `assemble_Fig2.R`, and `assemble_Fig3.R` later if the
  final manuscript figures are assembled from separate panel outputs.

# Provenance

The implemented algorithms were inspected and adapted from upstream source;
this is not a clean-room claim. Retained author notices are in `NOTICE`, and the
combined adaptation is distributed under GPL-3.0-or-later.

| Source | Pin / audited functions |
| --- | --- |
| [csaw](https://github.com/bioc/csaw/tree/ce50500871243747fae230c752dc2d4fec5093ee) | 1.46.0; mergeWindows, combineTests |
| [edgeR](https://github.com/bioc/edgeR/tree/27ba7b39d819e68df2ecdea9f32a4b6039b2a0b5) | 4.10.7; calcNormFactors, dispCoxReid, glmQLFit, glmQLFTest |
| [limma](https://bioconductor.org/packages/limma/) | Runtime 3.68.4; fitFDist, squeezeVar, legacy untrended mode |

The NB log-link solver is a Python Fisher-scoring implementation with step
halving. It does not copy edgeR's C optimizer. Numerical comparisons are scoped
to the recorded fixtures and tolerances. Region p-values implement the
unweighted Simes procedure; no metapod source is redistributed.

References:

- Lun ATL, Smyth GK (2016). csaw: a Bioconductor package for differential binding
  analysis of ChIP-seq data using sliding windows. Nucleic Acids Research 44:e45.
  https://doi.org/10.1093/nar/gkv1191
- Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package for
  differential expression analysis of digital gene expression data.
  https://doi.org/10.1093/bioinformatics/btp616
- Robinson MD, Oshlack A (2010). A scaling normalization method for differential
  expression analysis of RNA-seq data. https://doi.org/10.1186/gb-2010-11-3-r25
- Smyth GK (2004). Linear models and empirical Bayes methods for assessing
  differential expression in microarray experiments.
  https://doi.org/10.2202/1544-6115.1027

All reference datasets are synthetic, generated with seed 20261022. Runtime
versions and settings are recorded in `reference/provenance.json`; file hashes
are in `reference/manifest.json`. No patient data or private paths are committed.

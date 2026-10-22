# Validation and limits

`reference/generate.R` generates real reference outputs from csaw 1.46.0,
edgeR 4.10.7, and limma 3.68.4 under R 4.6.1. The checked-in fixtures cover
100-window two-group, covariate, sparse/zero, and estimated-dispersion cases.
They also test TMM, common Cox–Reid dispersion, merging, Simes/BH regional tests,
and BED/GenomicRanges overlap conversion.

Tests compare intermediate deviances, posterior dispersions, prior parameters,
degrees of freedom, F statistics, p-values and adjusted p-values. Test quantities
use relative tolerance 2e-5 plus absolute 2e-7. TMM uses relative 1e-12; common
dispersion uses relative 1e-4 plus absolute 1e-6. Region probabilities use 1e-14.

Unshrunk log2 fold changes on non-boundary fits use relative 2e-5 plus absolute
2e-4: edgeR's generic GLM solver stops at a different convergence criterion than
the Python solver. The measured largest discrepancy on the covariate fixture is
about 1.14e-4 log2 units. Boundary coefficients approximate infinite effects and
are intentionally not tested as finite, identifiable quantities; their final
test statistics are compared. No universal numerical tolerance is asserted.

`python reference/compare.py` refreshes the measured JSON and Markdown report,
including decisions at FDR 0.05. Passing these fixtures establishes bounded
numerical agreement, not general biological validity. Broad null/power studies,
real ChIP-seq datasets, donor/design robustness, and independent scientific
review remain needed before a stable research release.

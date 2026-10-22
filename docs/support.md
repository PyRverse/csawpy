# Supported surface

| Operation | Status and boundary |
| --- | --- |
| Count fragment/window overlaps | Implemented, integer counts; BED half-open intervals; GenomicRanges reference checked |
| Generate sliding windows | Implemented; complete windows only, no trailing short window |
| TMM normalization | Weighted TMM, edgeR reference checked; estimate from background bins |
| Abundance filtering | Simple CPM/sample-count rule; not csaw filterWindowsGlobal |
| NB dispersion | Supplied scalar/vector, or unweighted common Cox–Reid estimate on all eligible features |
| GLM | Full-rank numeric design with positive residual df, one numeric contrast |
| QL tests | Legacy, nonrobust, untrended, no coefficient shrinkage; Poisson bound |
| Region merging | Positional single-linkage, no strand or maximum-width split |
| Region tests | Unweighted Simes then BH; no direction/representative-test summaries |
| Current edgeR v4 default QL | Deferred; never silently substituted |
| Robust moderation, tagwise/trended dispersion, weights | Deferred |
| BAM I/O, annotation, plots | Deferred; preprocessed fragments/counts accepted |

Invalid counts, rank-deficient designs, incompatible dimensions, nonpositive
library sizes, and invalid p-values raise errors. Numerical nonconvergence raises
an error rather than returning a successful-looking test. All-zero windows within
an informative dataset have p=1; an entirely zero dataset cannot estimate the
moderation prior and raises an error. Nonzero windows with no residual df produce
a warning and must be inspected.

The common-dispersion optimizer uses all rows with at least five counts, with no
systematic subset. Default R csaw workflows often use more elaborate filtering,
normalization, dispersion, and robust QL options. A result under this explicitly
supported model is not a claim of equivalence to every csaw workflow.

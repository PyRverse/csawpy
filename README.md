# csawpy

**Window-based differential binding in Python.** Pronounced “seesaw pie”.

Part of [PyRverse](https://github.com/PyRverse). This is an **alpha, partial native
port** of a specified csaw/edgeR workflow. Python analyses require no R runtime.
The supported workflow has genuine R-generated reference tests; it has not had
independent scientific review or broad real-data validation.

```sh
python -m pip install .
```

Run the command from a local checkout. This version is not published on PyPI.

```python
import numpy as np
import pandas as pd
import csawpy as cs

counts = np.array(
    [
        [50, 60, 120, 130],
        [40, 45, 90, 95],
        [100, 90, 95, 100],
        [70, 75, 68, 72],
        [30, 35, 32, 34],
        [80, 85, 78, 83],
    ]
)
windows = pd.DataFrame(
    {
        "chrom": ["chr1"] * 6,
        "start": [0, 50, 500, 550, 1000, 1050],
        "end": [150, 200, 650, 700, 1150, 1200],
    }
)
design = np.column_stack([np.ones(4), [0, 0, 1, 1]])
result = cs.analyze(
    counts,
    windows,
    design,
    contrast=[0, 1],
    library_sizes=[1_000_000] * 4,
    dispersion=0.1,
)
print(result.regions[["chrom", "start", "end", "PValue", "FDR"]])
```

The example uses a supplied NB dispersion of 0.1 for illustration. Omitting it
estimates a common Cox–Reid dispersion. Select an appropriate design, dispersion
strategy, and abundance filter for the experiment before interpreting results.

## Supported workflow

1. Supply mapped, filtered fragments as BED-style `chrom/start/end` DataFrames
   and count overlaps with `count_windows`, or supply prepared window counts.
2. Compute `tmm_factors` on background bins and pass the resulting factors with
   mapped-fragment library sizes. Counts in overlapping windows do not sum to a
   library size. `filter_windows` provides a simple CPM filter.
3. Call `test_windows` or `analyze` with a numeric design and one contrast.
4. Merge **all tested windows** by location and combine their evidence using
   Simes, then BH across regions. Clustering only significant windows is invalid.

Coordinates are zero-based, half-open. Counts have windows in rows and samples
in columns. `region_ids` correspond to input window order, including unsorted
input. Both design and contrast are positional; align sample order explicitly.

## Statistical contract

The QL implementation targets **edgeR's legacy mode**, with `robust=FALSE`,
`abundance.trend=FALSE`, `prior.count=0`, and the Poisson lower bound enabled.
This is a supported older model in current edgeR, not its v4 default. Log fold
changes are unshrunk. Scalar or per-window NB dispersions are accepted; automatic
estimation currently supplies a common dispersion only.

Current limitations include BAM parsing, paired-mate reconstruction, csaw global
background filtering, robust/abundance-trended moderation, edgeR v4 adjusted
deviances, tagwise dispersion estimation, multi-df tests, weighted regional
tests, strand-aware merging, and max-width clustering. See
[the support matrix](docs/support.md) and [validation](docs/validation.md).

## Development

```sh
python -m pip install -e '.[dev]'
python -m pytest
python -m ruff check .
python -m build
python examples/quickstart.py
python reference/compare.py
```

Regenerate reference data with `Rscript reference/generate.R reference` in the
pinned R environment described in `reference/provenance.json`. R is required
only for regenerating reference outputs, not for ordinary testing or use.

GNU GPL-3.0-or-later. Original method authors are credited in [NOTICE](NOTICE)
and [provenance](docs/provenance.md). PyRverse is independent of Bioconductor and
the original package maintainers.

"""Small executable example; replace synthetic counts with experimental data."""

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
    counts, windows, design, [0, 1], library_sizes=[1_000_000] * 4, dispersion=0.1
)
print(result.regions.to_string(index=False))

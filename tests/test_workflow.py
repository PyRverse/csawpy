import json
from pathlib import Path

import numpy as np
import pandas as pd

import csawpy as cs


def test_region_results_stay_aligned_for_unsorted_input():
    fixture = Path(__file__).parents[1] / "reference" / "two-group-input.json"
    data = json.loads(fixture.read_text())
    counts = np.asarray(data["counts"])
    windows = pd.DataFrame(
        {
            "chrom": ["chr1"] * len(counts),
            "start": np.arange(len(counts)) * 200,
            "end": np.arange(len(counts)) * 200 + 150,
        }
    )
    order = np.arange(len(counts))[::-1]
    result = cs.analyze(
        counts[order],
        windows.iloc[order],
        data["design"],
        data["contrast"],
        library_sizes=data["library_sizes"],
        dispersion=0.15,
        merge_gap=0,
    )
    # Singleton regions sorted by genomic position must retain the corresponding p-value.
    np.testing.assert_allclose(
        result.regions.PValue, result.windows.PValue.to_numpy()[::-1]
    )
    np.testing.assert_array_equal(result.regions.start, windows.start)

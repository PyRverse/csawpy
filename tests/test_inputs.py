import numpy as np
import pandas as pd
import pytest

import csawpy as cs


def test_half_open_boundaries_and_order():
    w = pd.DataFrame({"chrom": ["chr1"] * 3, "start": [10, 0, 20], "end": [20, 10, 30]})
    f = pd.DataFrame({"chrom": ["chr1"] * 2, "start": [10, 20], "end": [20, 30]})
    np.testing.assert_array_equal(cs.count_windows(w, {"x": f}).x, [1, 0, 1])
    c = cs.merge_windows(w, tol=0)
    np.testing.assert_array_equal(c.ids, [1, 1, 1])
    assert c.regions.start.iloc[0] == 0
    assert c.regions.end.iloc[0] == 30


@pytest.mark.parametrize(
    "counts",
    [
        [[1, -1], [2, 3]],
        [[1, np.nan], [2, 3]],
        [[0.1, 2], [3, 4]],
        [[2**54, 1], [1, 2]],
        [],
    ],
)
def test_invalid_counts(counts):
    with pytest.raises(ValueError):
        cs.tmm_factors(counts)


def test_integer_precision_is_checked_before_float_conversion():
    with pytest.raises(ValueError, match="exact float64"):
        cs.tmm_factors(np.array([[2**53 + 1, 10], [10, 10]], dtype=np.int64))


def test_coordinate_precision_is_checked_before_float_conversion():
    intervals = pd.DataFrame(
        {"chrom": ["chr1"], "start": [2**53 + 1], "end": [2**53 + 10]}
    )
    with pytest.raises(ValueError, match="exact float64"):
        cs.merge_windows(intervals)


def test_integer_valued_float_window_width():
    result = cs.sliding_windows({"chr1": 301}, width=150.0, spacing=100.0)
    np.testing.assert_array_equal(result.start, [0, 100])


@pytest.mark.parametrize("pvalues", [[-0.1], [1.1], [np.nan], [np.inf]])
def test_invalid_pvalues(pvalues):
    with pytest.raises(ValueError):
        cs.combine_tests([1], pd.DataFrame({"PValue": pvalues}))


def test_singleton_and_simes():
    r = cs.combine_tests([1, 1, 2], pd.DataFrame({"PValue": [0.01, 0.04, 0.2]}))
    np.testing.assert_allclose(r.PValue, [0.02, 0.2])
    np.testing.assert_allclose(r.FDR, [0.04, 0.2])


def test_design_validation():
    y = np.array([[10, 12, 14, 16], [30, 28, 27, 31], [5, 8, 9, 8]])
    with pytest.raises(ValueError, match="full column rank"):
        cs.test_windows(
            y, np.ones((4, 2)), [0, 1], library_sizes=[100] * 4, dispersion=0.1
        )
    with pytest.raises(ValueError, match="contrast"):
        cs.test_windows(
            y, np.ones((4, 1)), [0], library_sizes=[100] * 4, dispersion=0.1
        )


def test_windows_generation_and_filter():
    w = cs.sliding_windows({"chr1": 301}, width=150, spacing=100)
    np.testing.assert_array_equal(w.start, [0, 100])
    np.testing.assert_array_equal(w.end, [150, 250])
    keep = cs.filter_windows(
        [[10, 10, 0], [0, 0, 0]], [1e6] * 3, min_cpm=5, min_samples=2
    )
    np.testing.assert_array_equal(keep, [True, False])


def test_empty_interval_results():
    w = pd.DataFrame(columns=["chrom", "start", "end"])
    assert cs.merge_windows(w).regions.empty
    assert cs.count_windows(w, {"a": w}).empty


def test_scaling_libraries_preserves_contrast():
    y = np.array([[10, 12, 24, 25], [30, 28, 27, 31], [5, 8, 9, 8]])
    x = np.column_stack([np.ones(4), [0, 0, 1, 1]])
    a = cs.test_windows(y, x, [0, 1], library_sizes=[100] * 4, dispersion=0.1)
    b = cs.test_windows(y, x, [0, 1], library_sizes=[1e6] * 4, dispersion=0.1)
    np.testing.assert_allclose(
        a[["logFC", "PValue"]], b[["logFC", "PValue"]], atol=1e-7
    )

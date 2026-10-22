"""Native window-based differential binding with explicit statistical scope."""

from ._intervals import WindowClusters, count_windows, merge_windows, sliding_windows
from ._statistics import (
    BindingResult,
    analyze,
    combine_tests,
    estimate_dispersion,
    filter_windows,
    test_windows,
    tmm_factors,
)

__version__ = "0.1.0a1"
__all__ = [
    "BindingResult",
    "WindowClusters",
    "analyze",
    "combine_tests",
    "count_windows",
    "estimate_dispersion",
    "filter_windows",
    "merge_windows",
    "sliding_windows",
    "test_windows",
    "tmm_factors",
]

"""Genomic intervals use BED's zero-based, half-open coordinates throughout.

Clustering follows csaw::mergeWindows (Aaron Lun), without max.width or strands.
See NOTICE and docs/provenance.md for attribution and supported settings.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd


def _intervals(value: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(value, pd.DataFrame):
        raise TypeError("intervals must be a DataFrame with chrom, start, end")
    if not {"chrom", "start", "end"}.issubset(value):
        raise ValueError("intervals require chrom, start, end columns")
    out = value[["chrom", "start", "end"]].copy()
    if out.isna().any().any():
        raise ValueError("intervals cannot contain missing values")
    if not out.chrom.map(lambda x: isinstance(x, str) and bool(x)).all():
        raise ValueError("chromosome names must be nonempty strings")
    for col in ("start", "end"):
        raw = out[col].to_numpy()
        if len(raw) and raw.dtype.kind not in "iuf":
            raise ValueError(
                "coordinates must be numeric integers, not strings or booleans"
            )
        if raw.dtype.kind in "iu" and np.any(raw > 2**53):
            raise ValueError("coordinates exceed the exact float64 integer range")
        a = out[col].to_numpy(dtype=float)
        if (
            np.any(~np.isfinite(a))
            or np.any(a != np.floor(a))
            or np.any(a < 0)
            or np.any(a > 2**53)
        ):
            raise ValueError(
                "coordinates must be nonnegative exactly representable integers"
            )
        out[col] = a.astype(np.int64)
    if (out.end <= out.start).any():
        raise ValueError("interval ends must exceed starts")
    return out


def sliding_windows(
    chrom_sizes: Mapping[str, int], width: int = 150, spacing: int = 50
) -> pd.DataFrame:
    """Generate full-width windows; incomplete terminal windows are omitted."""
    if (
        isinstance(width, bool)
        or isinstance(spacing, bool)
        or int(width) != width
        or int(spacing) != spacing
        or min(width, spacing) <= 0
    ):
        raise ValueError("width and spacing must be positive integers")
    width, spacing = int(width), int(spacing)
    rows = []
    for chrom, size in chrom_sizes.items():
        if isinstance(size, bool) or int(size) != size or size < 0:
            raise ValueError("chromosome lengths must be nonnegative integers")
        rows.extend(
            (chrom, s, s + width) for s in range(0, int(size) - width + 1, int(spacing))
        )
    return _intervals(pd.DataFrame(rows, columns=["chrom", "start", "end"]))


def count_windows(
    windows: pd.DataFrame, fragments: Mapping[str, pd.DataFrame]
) -> pd.DataFrame:
    """Count each fragment once in each overlapping window.

    Supply already filtered/extended fragments. Paired-end mates must be collapsed
    to one fragment before this call. No MAPQ/duplicate/strand inference is done.
    """
    w = _intervals(windows)
    if not fragments:
        raise ValueError("at least one sample of fragments is required")
    output = np.zeros((len(w), len(fragments)), dtype=np.int64)
    for j, f in enumerate(fragments.values()):
        f = _intervals(f)
        for chrom in pd.unique(w.chrom):
            mask = (w.chrom == chrom).to_numpy()
            frag = f[f.chrom == chrom]
            starts = np.sort(frag.start.to_numpy())
            ends = np.sort(frag.end.to_numpy())
            output[mask, j] = np.searchsorted(
                starts, w.loc[mask, "end"], side="left"
            ) - np.searchsorted(ends, w.loc[mask, "start"], side="right")
    return pd.DataFrame(output, index=windows.index, columns=list(fragments))


@dataclass(frozen=True)
class WindowClusters:
    ids: np.ndarray
    regions: pd.DataFrame


def merge_windows(windows: pd.DataFrame, tol: int = 100) -> WindowClusters:
    """Single-linkage merge by genomic distance, independent of test p-values.

    IDs are one-based and aligned to the original input order. Negative tolerance
    requires overlap. Chromosomes retain their order of first appearance.
    """
    w = _intervals(windows)
    if isinstance(tol, bool) or not np.isfinite(tol) or int(tol) != tol:
        raise ValueError("tol must be a finite integer")
    ids = np.zeros(len(w), dtype=int)
    regions = []
    for chrom in pd.unique(w.chrom):
        idx = np.flatnonzero(w.chrom.to_numpy() == chrom)
        order = idx[np.lexsort((w.end.to_numpy()[idx], w.start.to_numpy()[idx]))]
        left = right = None
        for i in order:
            start, end = int(w.iloc[i].start), int(w.iloc[i].end)
            if right is None or start - right > tol:
                regions.append([chrom, start, end])
                left, right = start, end
            else:
                right = max(right, end)
                regions[-1] = [chrom, left, right]
            ids[i] = len(regions)
    return WindowClusters(ids, pd.DataFrame(regions, columns=["chrom", "start", "end"]))

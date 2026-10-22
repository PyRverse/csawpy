"""Statistical core; adaptations of csaw/edgeR/limma are credited in NOTICE.

The supported QL mode corresponds to edgeR legacy=True, robust=False,
abundance.trend=False, prior.count=0. Current edgeR v4 defaults differ.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import linalg, optimize, special, stats

from ._intervals import merge_windows


def _counts(counts):
    original = np.asarray(counts)
    if original.dtype.kind not in "iuf":
        raise ValueError("counts must be numeric integers, not strings or booleans")
    if original.dtype.kind in "iu" and np.any(original > 2**53):
        raise ValueError("counts exceed the exact float64 integer range")
    y = np.asarray(counts, dtype=float)
    if y.ndim != 2 or min(y.shape) == 0 or not np.all(np.isfinite(y)) or np.any(y < 0):
        raise ValueError(
            "counts must be a nonempty finite nonnegative feature-by-sample matrix"
        )
    if np.any(y != np.floor(y)) or np.any(y > 2**53):
        raise ValueError("counts must be integers no greater than 2**53")
    return y


def _positive_vector(values, n, label):
    a = np.asarray(values, dtype=float)
    if a.shape != (n,) or np.any(~np.isfinite(a)) or np.any(a <= 0):
        raise ValueError(f"{label} must contain {n} finite positive values")
    return a


def _design(design, n):
    x = np.asarray(design, dtype=float)
    if x.ndim != 2 or x.shape[0] != n or x.shape[1] == 0 or not np.all(np.isfinite(x)):
        raise ValueError("design must be a finite samples-by-coefficients matrix")
    if np.linalg.matrix_rank(x) != x.shape[1] or x.shape[1] >= n:
        raise ValueError(
            "design must have full column rank and positive residual degrees of freedom"
        )
    return x


def _bh(p):
    p = np.asarray(p, dtype=float)
    valid = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    v = p[valid]
    if not len(v):
        return out
    order = np.argsort(v, kind="stable")
    adjusted = np.minimum.accumulate(
        (v[order] * len(v) / np.arange(1, len(v) + 1))[::-1]
    )[::-1]
    back = np.empty_like(v)
    back[order] = np.minimum(1, adjusted)
    out[valid] = back
    return out


def tmm_factors(
    counts, library_sizes=None, *, reference=None, logratio_trim=0.3, sum_trim=0.05
):
    """Weighted TMM factors with product one, matching edgeR's TMM settings.

    For a binding workflow, estimate these on background bins with mapped-fragment
    library sizes, then apply the factors to the selected windows.
    """
    y = _counts(counts)
    lib = _positive_vector(
        y.sum(axis=0) if library_sizes is None else library_sizes,
        y.shape[1],
        "library_sizes",
    )
    if not (0 <= logratio_trim < 0.5 and 0 <= sum_trim < 0.5):
        raise ValueError("trimming fractions must lie in [0, 0.5)")
    y = y[np.any(y > 0, axis=1)]
    if len(y) == 0 or y.shape[1] == 1:
        return np.ones(len(lib))
    if reference is None:
        q75 = np.quantile(y, 0.75, axis=0) / lib
        reference = (
            np.argmax(np.sqrt(y).sum(axis=0))
            if np.median(q75) < 1e-20
            else np.argmin(np.abs(q75 - q75.mean()))
        )
    if (
        isinstance(reference, bool)
        or int(reference) != reference
        or not 0 <= reference < y.shape[1]
    ):
        raise ValueError("reference must be a zero-based sample column index")
    factors = np.ones(y.shape[1])
    r = y[:, int(reference)]
    for j, o in enumerate(y.T):
        good = (o > 0) & (r > 0)
        obs, ref = o[good], r[good]
        if len(obs) == 0:
            continue
        m = np.log2((obs / lib[j]) / (ref / lib[int(reference)]))
        a = (np.log2(obs / lib[j]) + np.log2(ref / lib[int(reference)])) / 2
        v = (lib[j] - obs) / lib[j] / obs + (lib[int(reference)] - ref) / lib[
            int(reference)
        ] / ref
        if np.max(np.abs(m)) < 1e-6:
            continue
        if np.any(v <= 0):
            raise ValueError(
                "library_sizes must exceed individual feature counts for weighted TMM"
            )
        n = len(m)
        lo_m, lo_a = np.floor(n * logratio_trim) + 1, np.floor(n * sum_trim) + 1
        mr, ar = stats.rankdata(m), stats.rankdata(a)
        keep = (mr >= lo_m) & (mr <= n + 1 - lo_m) & (ar >= lo_a) & (ar <= n + 1 - lo_a)
        if np.any(keep):
            factors[j] = 2 ** (np.sum(m[keep] / v[keep]) / np.sum(1 / v[keep]))
    return factors / np.exp(np.mean(np.log(factors)))


def filter_windows(counts, library_sizes, *, min_cpm=1.0, min_samples=2):
    """Simple abundance filter, explicitly distinct from csaw filterWindowsGlobal."""
    y = _counts(counts)
    lib = _positive_vector(library_sizes, y.shape[1], "library_sizes")
    if (
        not np.isfinite(min_cpm)
        or min_cpm < 0
        or isinstance(min_samples, bool)
        or int(min_samples) != min_samples
        or not 1 <= min_samples <= y.shape[1]
    ):
        raise ValueError("invalid CPM or sample threshold")
    return np.sum(y / lib * 1e6 >= min_cpm, axis=1) >= min_samples


def _loglik(y, mu, phi):
    if phi <= 1e-12:
        return np.sum(special.xlogy(y, mu) - mu - special.gammaln(y + 1))
    r = 1 / phi
    return np.sum(
        special.gammaln(y + r)
        - special.gammaln(r)
        - special.gammaln(y + 1)
        + special.xlogy(y, mu / (mu + r))
        - r * np.log1p(mu / r)
    )


def _fit(y, x, offset, phi):
    if x.shape[1] == 0:
        mu = np.exp(offset)
        return np.zeros(0), mu, _loglik(y, mu, phi)
    beta = np.linalg.lstsq(x, np.log(y + 0.5) - offset, rcond=None)[0]
    mu = np.exp(np.clip(offset + x @ beta, -40, 40))
    ll = _loglik(y, mu, phi)
    for _ in range(150):
        score = x.T @ ((y - mu) / (1 + phi * mu))
        if np.max(np.abs(score)) < 1e-8:
            break
        w = mu / (1 + phi * mu)
        step = np.linalg.lstsq(x.T @ (w[:, None] * x), score, rcond=None)[0]
        scale = 1.0
        accepted = False
        for _ in range(30):
            candidate = beta + scale * step
            new_mu = np.exp(np.clip(offset + x @ candidate, -40, 40))
            new_ll = _loglik(y, new_mu, phi)
            if new_ll >= ll - 1e-10:
                beta, mu, ll = candidate, new_mu, new_ll
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            if np.max(np.abs(score)) < 1e-5:
                break
            raise RuntimeError("negative-binomial fit failed to converge")
    else:
        raise RuntimeError("negative-binomial fit exceeded iteration limit")
    return beta, mu, ll


def _deviance(y, mu, phi):
    if phi <= 1e-12:
        return max(0.0, 2 * np.sum(special.xlogy(y, y / mu) - y + mu))
    return max(
        0.0,
        2
        * np.sum(
            special.xlogy(y, y / mu)
            - (y + 1 / phi) * np.log((1 + phi * y) / (1 + phi * mu))
        ),
    )


def _moderate(var, df):
    # limma squeezeVar(...legacy=TRUE,robust=FALSE,covariate=NULL).
    if len(var) < 3:
        return var.copy(), np.zeros(len(var)), var.copy()
    valid = (df > 0) & np.isfinite(var)
    if not np.any(valid):
        raise ValueError("no windows with residual degrees of freedom")
    v, d = var[valid], df[valid]
    if len(v) == 1:
        return var.copy(), np.zeros(len(var)), var.copy()
    med = np.median(v)
    v = np.maximum(v, 1e-5 * (med if med > 0 else 1))
    e = np.log(v) + np.log(d / 2) - special.digamma(d / 2)
    ev = np.var(e, ddof=1) - np.mean(special.polygamma(1, d / 2))
    if ev <= 0:
        prior_df = np.inf
        prior_var = np.mean(v)
        post = np.full(len(var), prior_var)
    else:
        root = optimize.brentq(lambda z: special.polygamma(1, np.exp(z)) - ev, -30, 40)
        prior_df = 2 * np.exp(root)
        prior_var = np.exp(
            e.mean() + special.digamma(prior_df / 2) - np.log(prior_df / 2)
        )
        post = (df * var + prior_df * prior_var) / (df + prior_df)
    return post, np.full(len(var), prior_df), np.full(len(var), prior_var)


def estimate_dispersion(
    counts, design, library_sizes, *, norm_factors=None, min_row_sum=5
):
    """Common NB dispersion by Cox-Reid adjusted profile likelihood.

    Matches the unweighted, unsubsampled edgeR dispCoxReid objective on [0,4].
    This is a common dispersion, not edgeR's trended/tagwise dispersion pipeline.
    """
    y = _counts(counts)
    x = _design(design, y.shape[1])
    lib = _positive_vector(library_sizes, y.shape[1], "library_sizes")
    norm = (
        np.ones(y.shape[1])
        if norm_factors is None
        else _positive_vector(norm_factors, y.shape[1], "norm_factors")
    )
    if not np.isfinite(min_row_sum) or min_row_sum < 0:
        raise ValueError("min_row_sum must be finite and nonnegative")
    y = y[y.sum(axis=1) >= max(min_row_sum, 1)]
    if not len(y):
        raise ValueError("no windows pass min_row_sum")
    offset = np.log(lib * norm)

    def objective(fourth_root):
        phi = fourth_root**4
        value = 0.0
        for row in y:
            _, mu, ll = _fit(row, x, offset, phi)
            info = x.T @ ((mu / (1 + phi * mu))[:, None] * x)
            sign, logdet = np.linalg.slogdet(info)
            if sign <= 0:
                return np.inf
            value -= ll - logdet / 2
        return value

    result = optimize.minimize_scalar(
        objective, bounds=(1e-4, 4**0.25), method="bounded", options={"xatol": 1e-5}
    )
    if not result.success:
        raise RuntimeError("dispersion optimization failed")
    return float(result.x**4)


def test_windows(
    counts,
    design,
    contrast,
    *,
    library_sizes,
    dispersion=None,
    norm_factors=None,
    poisson_bound=True,
):
    """Legacy nonrobust, untrended QL F-tests for a single numeric contrast.

    With dispersion=None, estimate a common Cox-Reid dispersion. This does not
    implement current edgeR v4 adjusted deviances or robust/trended moderation.
    All-zero windows receive p=1; failed numerical fits raise an error.
    """
    y = _counts(counts)
    x = _design(design, y.shape[1])
    c = np.asarray(contrast, dtype=float)
    if c.shape != (x.shape[1],) or not np.all(np.isfinite(c)) or not np.any(c != 0):
        raise ValueError("contrast must be a finite nonzero coefficient vector")
    lib = _positive_vector(library_sizes, y.shape[1], "library_sizes")
    norm = (
        np.ones(y.shape[1])
        if norm_factors is None
        else _positive_vector(norm_factors, y.shape[1], "norm_factors")
    )
    if dispersion is None:
        dispersion = estimate_dispersion(y, x, lib, norm_factors=norm)
    phi = np.asarray(dispersion, dtype=float)
    if phi.ndim == 0:
        phi = np.full(len(y), phi)
    if phi.shape != (len(y),) or np.any(~np.isfinite(phi)) or np.any(phi < 0):
        raise ValueError("dispersion must be a nonnegative scalar or per-window vector")
    offset = np.log(lib * norm)
    reduced = x @ linalg.null_space(c[None, :])
    df = np.full(len(y), y.shape[1] - x.shape[1], dtype=float)
    dev, lr, logfc = np.zeros((3, len(y)))
    fits = []
    for i, row in enumerate(y):
        if not np.any(row):
            df[i] = 0
            fits.append(np.zeros(y.shape[1]))
            continue
        beta, mu, ll = _fit(row, x, offset, phi[i])
        _, _, ll0 = _fit(row, reduced, offset, phi[i])
        dev[i] = _deviance(row, mu, phi[i])
        lr[i] = max(0, 2 * (ll - ll0))
        logfc[i] = c @ beta / np.log(2)
        nonzero = ~((mu < 1e-4) & (row < 1e-4))
        df[i] = nonzero.sum() - (
            np.linalg.matrix_rank(x[nonzero]) if np.any(nonzero) else 0
        )
        fits.append(mu)
    raw = np.divide(dev, df, out=np.zeros(len(y)), where=df > 0)
    post, prior_df, prior_var = _moderate(raw, df)
    total_df = np.minimum(prior_df + df, len(y) * (y.shape[1] - x.shape[1]))
    fstat = np.divide(lr, post, out=np.where(lr > 0, np.inf, 0.0), where=post > 0)
    pval = stats.f.sf(fstat, 1, total_df)
    zero = y.sum(axis=1) == 0
    pval[zero] = 1
    if poisson_bound:
        for i, row in enumerate(y):
            if not zero[i] and np.any(post[i] * (1 + phi[i] * fits[i]) < 1):
                _, _, ll = _fit(row, x, offset, 0)
                _, _, ll0 = _fit(row, reduced, offset, 0)
                pval[i] = max(pval[i], stats.chi2.sf(max(0, 2 * (ll - ll0)), 1))
    if np.any((df == 0) & ~zero):
        warnings.warn(
            "Some nonzero windows have no residual degrees of freedom; inspect df_residual",
            RuntimeWarning,
            stacklevel=2,
        )
    index = counts.index if isinstance(counts, pd.DataFrame) else None
    table = pd.DataFrame(
        {
            "logFC": logfc,
            "F": fstat,
            "PValue": pval,
            "FDR": _bh(pval),
            "dispersion": phi,
            "deviance": dev,
            "df_residual": df,
            "ql_dispersion": raw,
            "ql_posterior": post,
            "prior_df": prior_df,
            "prior_variance": prior_var,
            "df_total": total_df,
        },
        index=index,
    )
    table.attrs["method"] = "edgeR legacy nonrobust untrended QL; unshrunk logFC"
    return table


def combine_tests(ids, table):
    """Unweighted Simes p-values per region, then BH across regions.

    Clusters must be defined independently of p-values. Do not cluster only
    significant windows. Direction summaries and weighted tests are deferred.
    """
    ids = np.asarray(ids)
    if ids.ndim != 1 or len(ids) != len(table) or pd.isna(ids).any():
        raise ValueError("one nonmissing region ID is required per test")
    if "PValue" not in table:
        raise ValueError("table needs a PValue column")
    p = np.asarray(table.PValue, dtype=float)
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be finite and between zero and one")
    rows = []
    for key in pd.unique(ids):
        pos = np.flatnonzero(ids == key)
        ordered = np.sort(p[pos])
        simes = min(1.0, float(np.min(ordered * len(pos) / np.arange(1, len(pos) + 1))))
        rows.append((key, len(pos), simes))
    out = pd.DataFrame(rows, columns=["region_id", "num_tests", "PValue"])
    out["FDR"] = _bh(out.PValue)
    return out


@dataclass(frozen=True)
class BindingResult:
    windows: pd.DataFrame
    regions: pd.DataFrame
    region_ids: np.ndarray


def analyze(
    counts,
    windows,
    design,
    contrast,
    *,
    library_sizes,
    dispersion=None,
    norm_factors=None,
    merge_gap=100,
):
    """Test already abundance-filtered windows and aggregate every tested window."""
    if len(counts) != len(windows):
        raise ValueError("counts and windows must have the same number of rows")
    if isinstance(counts, pd.DataFrame) and not counts.index.equals(windows.index):
        raise ValueError("counts and windows DataFrame indices must match in order")
    clusters = merge_windows(windows, tol=merge_gap)
    tests = test_windows(
        counts,
        design,
        contrast,
        library_sizes=library_sizes,
        dispersion=dispersion,
        norm_factors=norm_factors,
    )
    combined = combine_tests(clusters.ids, tests)
    regions = clusters.regions.assign(
        region_id=np.arange(1, len(clusters.regions) + 1)
    ).merge(combined, on="region_id", validate="one_to_one", sort=False)
    return BindingResult(tests, regions, clusters.ids)

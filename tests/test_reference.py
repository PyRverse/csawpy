import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import csawpy as cs

REF = Path(__file__).resolve().parents[1] / "reference"


@pytest.mark.parametrize("case", ["two-group", "covariate", "zeros", "estimated"])
def test_edgeR_legacy_ql_reference(case):
    data = json.loads((REF / f"{case}-input.json").read_text())
    expected = pd.read_csv(REF / f"{case}-expected.csv")
    result = cs.test_windows(**data)
    for col in [
        "F",
        "PValue",
        "FDR",
        "deviance",
        "ql_posterior",
        "prior_variance",
        "prior_df",
        "df_residual",
        "df_total",
    ]:
        np.testing.assert_allclose(
            result[col], expected[col], rtol=2e-5, atol=2e-7, err_msg=col
        )
    # Boundary fits have non-identifiable finite approximations to infinite logFC.
    valid = np.all(np.asarray(data["counts"]) > 0, axis=1)
    np.testing.assert_allclose(
        result.logFC[valid], expected.logFC[valid], rtol=2e-5, atol=2e-4
    )


def test_reference_normalization_and_common_dispersion():
    data = json.loads((REF / "two-group-input.json").read_text())
    expected = json.loads((REF / "normalization.json").read_text())
    np.testing.assert_allclose(
        cs.tmm_factors(data["counts"], data["library_sizes"]),
        expected["tmm"],
        rtol=1e-12,
    )
    phi = cs.estimate_dispersion(data["counts"], data["design"], data["library_sizes"])
    np.testing.assert_allclose(phi, expected["common_dispersion"], rtol=1e-4, atol=1e-6)


def test_csaw_region_reference():
    data = json.loads((REF / "regions.json").read_text())
    clusters = cs.merge_windows(pd.DataFrame(data["windows"]), tol=50)
    np.testing.assert_array_equal(clusters.ids, data["ids"])
    pd.testing.assert_frame_equal(
        clusters.regions, pd.DataFrame(data["regions"]), check_dtype=False
    )
    combined = cs.combine_tests(
        clusters.ids, pd.DataFrame({"PValue": data["pvalues"]})
    ).sort_values("region_id")
    expected = pd.DataFrame(data["combined"])
    np.testing.assert_array_equal(combined.num_tests, expected["num.tests"])
    np.testing.assert_allclose(
        combined[["PValue", "FDR"]], expected[["PValue", "FDR"]], rtol=1e-14
    )


def test_genomicranges_fragment_overlap_reference():
    data = json.loads((REF / "overlaps.json").read_text())
    result = cs.count_windows(
        pd.DataFrame(data["windows"]), {"sample": pd.DataFrame(data["fragments"])}
    )
    np.testing.assert_array_equal(result["sample"], data["counts"])

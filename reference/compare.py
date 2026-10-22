"""Recompute measured parity reports; requires an installed csawpy."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import csawpy

root = Path(__file__).resolve().parent
for fixture in root.iterdir():
    if fixture.suffix in {".R", ".csv", ".json"}:
        fixture.write_bytes(fixture.read_bytes().replace(b"\r\n", b"\n"))
report = {
    "package_version": csawpy.__version__,
    "reference": json.loads((root / "provenance.json").read_text()),
    "cases": {},
}
for case in ("two-group", "covariate", "zeros", "estimated"):
    data = json.loads((root / f"{case}-input.json").read_text())
    reference = pd.read_csv(root / f"{case}-expected.csv")
    actual = csawpy.test_windows(**data)
    cols = [
        "F",
        "PValue",
        "FDR",
        "deviance",
        "ql_posterior",
        "prior_variance",
        "df_residual",
        "df_total",
    ]
    errors = {col: float(np.max(np.abs(actual[col] - reference[col]))) for col in cols}
    decisions = int(np.sum((actual.FDR < 0.05) != (reference.FDR < 0.05)))
    report["cases"][case] = {
        "windows": len(actual),
        "max_absolute_errors": errors,
        "fdr_0_05_decision_disagreements": decisions,
    }
output = root.parent / "reports"
output.mkdir(exist_ok=True)
(output / "parity.json").write_text(
    json.dumps(report, indent=2, allow_nan=False) + "\n"
)
lines = [
    "# Measured R reference comparisons",
    "",
    "| Case | Windows | Max p error | Max FDR error | FDR .05 disagreements |",
    "| --- | ---: | ---: | ---: | ---: |",
]
for case, record in report["cases"].items():
    errors = record["max_absolute_errors"]
    lines.append(
        f"| {case} | {record['windows']} | {errors['PValue']:.3g} | {errors['FDR']:.3g} | {record['fdr_0_05_decision_disagreements']} |"
    )
lines += [
    "",
    "See docs/validation.md for supported settings and acceptance criteria. These synthetic fixture results do not establish general statistical calibration.",
]
(output / "parity.md").write_text("\n".join(lines) + "\n")
files = sorted(
    p
    for p in root.iterdir()
    if p.suffix in {".R", ".csv", ".json"} and p.name != "manifest.json"
)
manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print("\n".join(lines))

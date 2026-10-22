import hashlib
import json
from pathlib import Path


def test_reference_integrity():
    root = Path(__file__).parents[1] / "reference"
    expected = json.loads((root / "manifest.json").read_text())
    for name, digest in expected.items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest, name

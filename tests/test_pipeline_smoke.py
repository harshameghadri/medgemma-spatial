"""Smoke tests for data loading utilities — no GPU required."""
import json
import os
import tempfile

from src.report_generation.prompt_builder import load_features_json


def test_load_features_json_roundtrip():
    """JSON written to disk should be readable back with correct values."""
    data = {
        "annotation": {"n_spots": 100},
        "spatial_heterogeneity": {},
        "uncertainty": {},
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as f:
        json.dump(data, f)
        path = f.name
    try:
        loaded = load_features_json(path)
        assert loaded["annotation"]["n_spots"] == 100
    finally:
        os.unlink(path)


def test_load_features_json_missing_file():
    """Missing file should raise FileNotFoundError."""
    import pytest

    with pytest.raises((FileNotFoundError, OSError)):
        load_features_json("/nonexistent/path/features.json")

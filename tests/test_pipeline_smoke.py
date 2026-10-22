"""End-to-end smoke test of the spatial pipeline on synthetic data (no model)."""

import json

import pytest

from src.report_generation.prompt_builder import generate_medgemma_prompt
from src.streamlit_adapter import (
    annotate_spatial_regions,
    calculate_spatial_heterogeneity,
)


@pytest.fixture(scope="module")
def pipeline_output():
    from tests.conftest import make_synthetic_tissue

    adata, annotation = annotate_spatial_regions(make_synthetic_tissue())
    spatial = calculate_spatial_heterogeneity(adata)
    return adata, annotation, spatial


@pytest.mark.slow
def test_pipeline_runs_and_returns_expected_keys(pipeline_output):
    adata, annotation, spatial = pipeline_output
    assert "cell_type" in adata.obs
    for key in ("n_spots", "n_clusters", "n_cell_types", "cell_type_counts"):
        assert key in annotation
    assert annotation["n_spots"] == adata.n_obs


@pytest.mark.slow
def test_features_are_json_serialisable(pipeline_output):
    _, annotation, spatial = pipeline_output
    features = {"annotation": annotation, "spatial_heterogeneity": spatial}
    json.dumps(features, default=str)


@pytest.mark.slow
def test_prompt_renders_from_pipeline_output(pipeline_output):
    _, annotation, spatial = pipeline_output
    features = {"annotation": annotation, "spatial_heterogeneity": spatial}
    prompt = generate_medgemma_prompt(features)
    assert "Report:" in prompt

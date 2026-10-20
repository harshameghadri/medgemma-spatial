"""Regression tests for the pipeline -> prompt data contract (docs/REVIEW.md S1-S7)."""

import json

import numpy as np
import pytest

from src.report_generation.prompt_builder import (
    generate_medgemma_prompt,
    is_immune_label,
    summarise_composition,
)
from src.streamlit_adapter import (
    annotate_spatial_regions,
    build_features,
    calculate_spatial_heterogeneity,
    looks_like_raw_counts,
)
from src.spatial_analysis.uncertainty_spatial_analysis import (
    compute_morans_i_with_uncertainty,
)
from tests.conftest import EPITHELIAL_GENES, T_CELL_GENES, make_synthetic_tissue


@pytest.fixture(scope="module")
def run():
    adata, annotation = annotate_spatial_regions(make_synthetic_tissue(), use_celltypist=False)
    spatial = calculate_spatial_heterogeneity(adata, n_permutations=50)
    return adata, build_features(adata, annotation, spatial)


def test_no_errors_on_clean_input(run):
    _, features = run
    assert features["errors"] == []


def test_feature_schema_is_flat_and_complete(run):
    _, features = run
    spatial = features["spatial_heterogeneity"]
    for key in ("morans_i_mean", "n_enriched_pairs", "entropy_mean"):
        assert isinstance(spatial[key], (int, float)), key
    json.dumps(features)


def test_t_cell_region_is_annotated(run):
    adata, _ = run
    t_region = adata.obs.loc[adata.obs["is_t_region"], "cell_type"]
    assert (t_region == "T_cells").mean() > 0.8


def test_s1_prompt_receives_nonzero_spatial_metrics(run):
    _, features = run
    spatial = features["spatial_heterogeneity"]
    prompt = generate_medgemma_prompt(features)
    assert "Moran's I): " in prompt and f"({spatial['morans_i_mean']:.2f})" in prompt
    assert spatial["n_enriched_pairs"] >= 2
    assert f"co-enriched cell-type pairs: {spatial['n_enriched_pairs']}" in prompt


def test_s2_t_cell_rich_tissue_is_not_cold(run):
    _, features = run
    comp = summarise_composition(features)
    assert comp["immune_fraction"] > 0.4
    assert "cold immune desert" not in generate_medgemma_prompt(features)


@pytest.mark.parametrize("label", ["T_cells", "NK_cells", "Plasma_cells", "Macrophages",
                                   "Tcm/Naive helper T cells", "Mast_cells", "Monocytes"])
def test_s2_immune_labels_recognised(label):
    assert is_immune_label(label)


@pytest.mark.parametrize("label", ["Epithelial", "CAF_generic", "Endothelial", "Stem cells", "Unknown"])
def test_s2_non_immune_labels_rejected(label):
    assert not is_immune_label(label)


def test_s3_confidence_is_measured_not_defaulted(run):
    adata, features = run
    assert features["annotation"]["mean_confidence"] == pytest.approx(
        float(adata.obs["annotation_confidence"].mean()))
    assert features["uncertainty"]["mean_prediction_entropy"] == pytest.approx(
        float(adata.obs["prediction_entropy"].mean()))


def test_s3_missing_annotation_yields_none_not_default():
    adata, annotation = annotate_spatial_regions(make_synthetic_tissue(), use_markers=False)
    features = build_features(adata, annotation, {"errors": []})
    assert features["annotation"]["mean_confidence"] is None
    assert features["uncertainty"]["mean_prediction_entropy"] is None


def test_s4_base_graph_not_overwritten(run):
    adata, _ = run
    degree = np.asarray((adata.obsp["spatial_connectivities"] > 0).sum(axis=1)).ravel()
    assert degree.max() == 6
    assert not any(k.startswith("nhood_r") for k in adata.obsp)


def test_s5_morans_i_detects_patterned_genes(run):
    adata, _ = run
    res = compute_morans_i_with_uncertainty(adata, genes=T_CELL_GENES[:5] + EPITHELIAL_GENES)
    assert min(res["morans_i"]) > 0.3
    assert all(p < 0.001 for p in res["p_value"])


def test_s6_raw_counts_detected_regardless_of_gene_count():
    adata = make_synthetic_tissue()[:, :600].copy()
    assert looks_like_raw_counts(adata.X)
    out, _ = annotate_spatial_regions(adata, use_celltypist=False)
    assert "log1p_norm" in out.layers
    assert not looks_like_raw_counts(out.layers["log1p_norm"])


def test_s7_missing_coordinates_reported_not_zeroed():
    adata, _ = annotate_spatial_regions(make_synthetic_tissue(), use_celltypist=False)
    del adata.obsm["spatial"]
    spatial = calculate_spatial_heterogeneity(adata)
    assert spatial["morans_i_mean"] is None
    assert any("spatial" in e.lower() for e in spatial["errors"])
    assert "not available" in generate_medgemma_prompt(
        build_features(adata, {"errors": []}, spatial))


def test_s7_nan_coordinates_reported():
    adata, _ = annotate_spatial_regions(make_synthetic_tissue(), use_celltypist=False)
    adata.obsm["spatial"][0, 0] = np.nan
    spatial = calculate_spatial_heterogeneity(adata)
    assert any("NaN" in e for e in spatial["errors"])


def test_s6_scaled_hvg_input_uses_full_raw_matrix():
    import scanpy as sc

    adata = make_synthetic_tissue()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=500)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata)
    out, metrics = annotate_spatial_regions(adata, use_celltypist=False)
    assert metrics["n_genes"] == adata.raw.n_vars
    assert metrics["errors"] == []
    t_region = out.obs.loc[out.obs["is_t_region"], "cell_type"]
    assert (t_region == "T_cells").mean() > 0.8

"""Unit tests for prompt_builder — no GPU required."""
import json
import pytest

from src.report_generation.prompt_builder import (
    generate_medgemma_prompt,
    evaluate_report_quality,
)

SAMPLE_FEATURES = {
    "annotation": {
        "cell_type_counts": {
            "T cells": 300,
            "Macrophages": 200,
            "Epithelial cells": 500,
            "Fibroblasts": 100,
        }
    },
    "spatial_heterogeneity": {"morans_i_mean": 0.45, "n_enriched_pairs": 3},
    "uncertainty": {"mean_prediction_entropy": 1.2},
}


def test_prompt_moa_focus_returns_string():
    """MoA-focused prompt should be a non-empty string."""
    p = generate_medgemma_prompt(SAMPLE_FEATURES, moa_focus=True)
    assert isinstance(p, str)
    assert len(p) > 200


def test_prompt_standard_returns_string():
    """Standard (non-MoA) prompt should contain generation instruction."""
    p = generate_medgemma_prompt(SAMPLE_FEATURES, moa_focus=False)
    assert isinstance(p, str)
    assert "Generate" in p


def test_evaluate_quality_returns_required_keys():
    """Quality evaluation must return all expected metric keys."""
    result = evaluate_report_quality(
        "This suggests an immune hot tumour.", SAMPLE_FEATURES
    )
    for key in ["word_count", "has_moa", "parroting_risk", "has_interpretation"]:
        assert key in result, f"Missing key: {key}"


def test_low_parroting_clean_report():
    """A well-written interpretive sentence should not be HIGH parroting risk."""
    clean = (
        "The tissue architecture suggests active immune infiltration "
        "mediated by IFN-γ signalling, consistent with an inflamed tumour phenotype."
    )
    result = evaluate_report_quality(clean, SAMPLE_FEATURES)
    assert result["parroting_risk"] in ("LOW", "MODERATE"), (
        f"Expected LOW/MODERATE, got {result['parroting_risk']}"
    )


def test_moa_prompt_contains_phenotype_hint():
    """MoA prompt should classify tissue phenotype."""
    hot_features = {
        "annotation": {
            "cell_type_counts": {
                "T cells": 800,
                "NK cells": 200,
                "Macrophages": 300,
                "Epithelial cells": 100,
            }
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.55, "n_enriched_pairs": 5},
        "uncertainty": {"mean_prediction_entropy": 0.8},
    }
    p = generate_medgemma_prompt(hot_features, moa_focus=True)
    assert "immune" in p.lower() or "hot" in p.lower() or "phenotype" in p.lower()

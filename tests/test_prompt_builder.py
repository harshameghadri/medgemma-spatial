"""Prompt builder must tolerate missing fields without crashing."""

import pytest

from src.report_generation.prompt_builder import generate_medgemma_prompt


@pytest.mark.parametrize("moa_focus", [True, False])
def test_prompt_handles_empty_features(moa_focus):
    prompt = generate_medgemma_prompt({}, moa_focus=moa_focus)
    assert isinstance(prompt, str) and len(prompt) > 100


def test_prompt_does_not_leak_raw_counts():
    features = {
        "annotation": {"cell_type_counts": {"T_cells": 12345, "Epithelial": 67890}},
        "spatial_heterogeneity": {"morans_i_mean": 0.5, "n_enriched_pairs": 10},
        "uncertainty": {"mean_prediction_entropy": 1.0},
    }
    prompt = generate_medgemma_prompt(features)
    assert "12345" not in prompt and "67890" not in prompt

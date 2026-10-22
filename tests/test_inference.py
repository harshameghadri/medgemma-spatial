"""MedGemma inference path exercised end-to-end with a tiny random Gemma 3 checkpoint."""

import numpy as np
import pytest

from src.report_generation import medgemma_multimodal as mm
from src.report_generation.prompt_builder import generate_medgemma_prompt

TINY_GEMMA3 = "trl-internal-testing/tiny-Gemma3ForConditionalGeneration"

FEATURES = {
    "annotation": {"cell_type_counts": {"T_cells": 300, "Epithelial": 500, "CAF_generic": 200}},
    "spatial_heterogeneity": {"morans_i_mean": 0.42, "n_enriched_pairs": 7, "entropy_mean": 0.9},
    "uncertainty": {"mean_prediction_entropy": 0.6},
}


def test_prompt_does_not_request_fabricated_citations():
    prompt = generate_medgemma_prompt(FEATURES)
    assert "never invent a reference" in prompt
    assert "cite author surname" not in prompt


def test_multimodal_prompt_tolerates_missing_metrics():
    prompt = mm.create_multimodal_prompt({})
    assert "not available" in prompt


def test_float_hires_image_converted_to_rgb():
    class FakeAdata:
        uns = {"spatial": {"lib": {"images": {"hires": np.random.rand(40, 60, 4).astype(np.float32)}}}}

    image = mm.load_hires_image_from_adata(FakeAdata())
    assert image.mode == "RGB" and image.size == (60, 40)


@pytest.mark.network
def test_textonly_generation_returns_only_new_tokens():
    report, meta = mm.generate_textonly_fallback(
        FEATURES, model_id=TINY_GEMMA3, device="cpu", max_new_tokens=8)
    assert isinstance(report, str)
    assert "SPATIAL TRANSCRIPTOMICS DATA SUMMARY" not in report
    assert meta["mode"] == "text_only" and meta["device"] == "cpu"


@pytest.mark.network
def test_multimodal_generation_with_hne_image():
    class FakeAdata:
        uns = {"spatial": {"lib": {"images": {"hires": (np.random.rand(64, 64, 3) * 255).astype(np.uint8)}}}}

    report, meta = mm.generate_multimodal_report(
        FakeAdata(), FEATURES, model_id=TINY_GEMMA3, device="cpu", max_new_tokens=8)
    assert isinstance(report, str)
    assert meta["image_size"] == [896, 896]


@pytest.mark.network
def test_model_is_loaded_once():
    first = mm.load_medgemma(TINY_GEMMA3, "cpu")
    second = mm.load_medgemma(TINY_GEMMA3, "cpu")
    assert first[1] is second[1]

"""Import-level checks that catch dependency drift before deploy."""

import importlib

import pytest

MODULES = [
    "src.streamlit_adapter",
    "src.report_generation.prompt_builder",
    "src.report_generation.medgemma_multimodal",
    "src.spatial_analysis.uncertainty_spatial_analysis",
    "src.spatial_analysis.gene_signature_scoring",
]


@pytest.mark.parametrize("module", MODULES)
def test_module_imports(module):
    importlib.import_module(module)


def test_transformers_supports_gemma3():
    """MedGemma 4B is a Gemma 3 checkpoint; transformers < 4.50 cannot load it."""
    importlib.import_module("transformers.models.gemma3")

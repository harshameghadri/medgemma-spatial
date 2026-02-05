#!/usr/bin/env python3
"""
MedGemma Multimodal Report Generation with H&E Images.

Integrates spatial transcriptomics features with H&E tissue images for
comprehensive pathology reports using MedGemma 1.5 multimodal capabilities.
"""

from PIL import Image
from transformers import pipeline
import torch
import numpy as np
from typing import Dict, Optional, Tuple
import sys
from pathlib import Path


def load_hires_image_from_adata(adata, library_id: Optional[str] = None) -> Image.Image:
    """
    Extract high-resolution H&E image from AnnData spatial slot.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial information
    library_id : str, optional
        Library ID for spatial data. If None, uses first available.

    Returns
    -------
    PIL.Image.Image
        High-resolution tissue image

    Raises
    ------
    KeyError
        If spatial data not found in AnnData
    """
    if 'spatial' not in adata.uns:
        raise KeyError("No spatial data found in adata.uns['spatial']")

    if library_id is None:
        library_id = list(adata.uns['spatial'].keys())[0]
        print(f"  Using library_id: {library_id}")

    try:
        img_array = adata.uns['spatial'][library_id]['images']['hires']
    except KeyError as e:
        available = list(adata.uns['spatial'][library_id]['images'].keys())
        raise KeyError(f"'hires' image not found. Available: {available}") from e

    if img_array.dtype == np.float64 or img_array.dtype == np.float32:
        img_array = (img_array * 255).astype(np.uint8)

    if img_array.ndim == 2:
        img_array = np.stack([img_array] * 3, axis=-1)

    image = Image.fromarray(img_array)
    print(f"  Loaded image: {image.size} pixels, mode: {image.mode}")

    return image


def preprocess_for_medgemma(
    image: Image.Image,
    target_size: int = 896,
    maintain_aspect: bool = True
) -> Image.Image:
    """
    Resize and preprocess H&E image for MedGemma input.

    MedGemma 1.5 uses SigLIP encoder with 896x896 input size.

    Parameters
    ----------
    image : PIL.Image.Image
        Input H&E tissue image
    target_size : int, default=896
        Target dimension (MedGemma expects 896x896)
    maintain_aspect : bool, default=True
        If True, maintain aspect ratio with padding. If False, stretch.

    Returns
    -------
    PIL.Image.Image
        Preprocessed image ready for MedGemma
    """
    if maintain_aspect:
        image_copy = image.copy()
        image_copy.thumbnail((target_size, target_size), Image.Resampling.LANCZOS)

        canvas = Image.new('RGB', (target_size, target_size), (255, 255, 255))
        offset_x = (target_size - image_copy.width) // 2
        offset_y = (target_size - image_copy.height) // 2
        canvas.paste(image_copy, (offset_x, offset_y))

        return canvas
    else:
        return image.resize((target_size, target_size), Image.Resampling.LANCZOS)


def create_multimodal_prompt(features: Dict) -> str:
    """
    Create prompt for multimodal MedGemma combining H&E image + spatial features.

    Uses anti-parroting strategy: exposes aggregated metrics only.

    Parameters
    ----------
    features : dict
        Spatial transcriptomics features (from uncertainty_spatial_analysis.py)

    Returns
    -------
    str
        Prompt for MedGemma multimodal input
    """
    annotation = features.get('annotation', {})
    cell_types = annotation.get('cell_type_counts', {})
    spatial = features.get('spatial_heterogeneity', {})
    uncertainty = features.get('uncertainty', {})

    n_spots = sum(cell_types.values())
    major_populations = {k: v for k, v in cell_types.items() if v / n_spots > 0.10}
    immune_types = {k: v for k, v in cell_types.items()
                    if any(x in k.lower() for x in ['t cell', 'b cell', 'macrophage',
                                                      'nk cell', 'dendritic', 'neutrophil'])}

    has_spatial_het = spatial.get('morans_i_mean', 0) > 0.3

    prompt = f"""You are a computational pathologist analyzing spatial transcriptomics data with corresponding H&E histology.

AVAILABLE INFORMATION:

H&E IMAGE:
- High-resolution tissue section provided above
- Use morphological features to interpret spatial patterns

MOLECULAR DATA:
- Total tissue spots: {n_spots}
- Major cell populations (>10% abundance): {len(major_populations)} types
- Immune cell diversity: {len(immune_types)} distinct immune populations
- Spatial organization: {'heterogeneous' if has_spatial_het else 'uniform'} (Moran's I = {spatial.get('morans_i_mean', 'N/A'):.3f})
- Prediction uncertainty: {uncertainty.get('mean_prediction_entropy', 'N/A'):.3f}

TASK:
Generate a 200-word clinical pathology report that:
1. Describes tissue morphology visible in the H&E image
2. Correlates morphological findings with molecular spatial patterns
3. Interprets biological significance of spatial organization
4. Provides clinical context

CRITICAL:
- DO NOT repeat exact numbers from the data
- SYNTHESIZE insights from both image and molecular data
- Focus on biological interpretation, not data enumeration

Report:"""

    return prompt


def generate_multimodal_report(
    adata,
    features: Dict,
    model_id: str = "google/medgemma-1.5-4b-it",
    device: str = "mps",
    max_new_tokens: int = 500,
    library_id: Optional[str] = None
) -> Tuple[str, Dict]:
    """
    Generate pathology report using BOTH H&E image and spatial features.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial images in adata.uns['spatial']
    features : dict
        Spatial transcriptomics features
    model_id : str, default="google/medgemma-1.5-4b-it"
        HuggingFace model identifier
    device : str, default="mps"
        Device for inference ('mps', 'cuda', or 'cpu')
    max_new_tokens : int, default=500
        Maximum tokens to generate
    library_id : str, optional
        Library ID for spatial data

    Returns
    -------
    tuple of (str, dict)
        Generated report and metadata

    Raises
    ------
    ImportError
        If transformers or PIL not available
    RuntimeError
        If model fails to load or generate
    """
    print("\n[MULTIMODAL REPORT GENERATION]")
    print("="*80)

    print("\n[1/5] Loading H&E image from AnnData...")
    try:
        he_image = load_hires_image_from_adata(adata, library_id=library_id)
    except Exception as e:
        print(f"  ✗ Failed to load image: {e}")
        raise

    print("\n[2/5] Preprocessing image for MedGemma...")
    he_image_processed = preprocess_for_medgemma(he_image, target_size=896)
    print(f"  ✓ Resized to: {he_image_processed.size}")

    print("\n[3/5] Creating multimodal prompt...")
    prompt = create_multimodal_prompt(features)
    print(f"  ✓ Prompt length: {len(prompt)} characters")

    print(f"\n[4/5] Loading MedGemma 1.5 multimodal model ({model_id})...")
    try:
        pipe = pipeline(
            "image-text-to-text",
            model=model_id,
            torch_dtype=torch.bfloat16,
            device=device
        )
        print(f"  ✓ Model loaded on device: {device}")
    except Exception as e:
        print(f"  ✗ Model loading failed: {e}")
        print("\nTROUBLESHOOTING:")
        print("  1. Check model ID is correct (google/medgemma-1.5-4b-it)")
        print("  2. Ensure transformers>=4.45.0 installed")
        print("  3. Try device='cpu' if MPS fails")
        raise

    print("\n[5/5] Generating multimodal report...")
    messages = [{
        "role": "user",
        "content": [
            {"type": "image", "image": he_image_processed},
            {"type": "text", "text": prompt}
        ]
    }]

    try:
        output = pipe(text=messages, max_new_tokens=max_new_tokens)
        report = output[0]['generated_text']
        print(f"  ✓ Generated {len(report)} characters")
    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        raise

    metadata = {
        'mode': 'multimodal',
        'model_id': model_id,
        'device': device,
        'image_size': he_image_processed.size,
        'report_length': len(report)
    }

    print("\n" + "="*80)
    print("[MULTIMODAL REPORT COMPLETE]")

    return report, metadata


def generate_textonly_fallback(
    features: Dict,
    model_id: str = "google/medgemma-4b-it",
    device: str = "mps",
    max_new_tokens: int = 500
) -> Tuple[str, Dict]:
    """
    Fallback text-only report generation if multimodal fails.

    Parameters
    ----------
    features : dict
        Spatial transcriptomics features
    model_id : str
        Text-only model ID
    device : str
        Device for inference
    max_new_tokens : int
        Maximum tokens to generate

    Returns
    -------
    tuple of (str, dict)
        Generated report and metadata
    """
    from transformers import AutoTokenizer, AutoModelForCausalLM

    print("\n[TEXT-ONLY FALLBACK MODE]")
    print("="*80)

    print(f"\n[1/3] Loading {model_id}...")
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForCausalLM.from_pretrained(
        model_id,
        torch_dtype=torch.bfloat16,
        device_map="auto",
        low_cpu_mem_usage=True
    )
    print(f"  ✓ Model loaded")

    print("\n[2/3] Creating text-only prompt...")
    prompt = create_multimodal_prompt(features)
    prompt = prompt.replace("H&E IMAGE:\n- High-resolution tissue section provided above\n- Use morphological features to interpret spatial patterns\n\n", "")
    print(f"  ✓ Prompt ready")

    print("\n[3/3] Generating report...")
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
    outputs = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        do_sample=False,
        pad_token_id=tokenizer.eos_token_id
    )

    full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
    report = full_output[len(prompt):].strip()
    print(f"  ✓ Generated {len(report)} characters")

    metadata = {
        'mode': 'text_only_fallback',
        'model_id': model_id,
        'device': device,
        'report_length': len(report)
    }

    print("\n" + "="*80)

    return report, metadata


if __name__ == "__main__":
    print("MedGemma Multimodal Module - Test")
    print("="*80)
    print("\nUsage:")
    print("  from src.report_generation.medgemma_multimodal import generate_multimodal_report")
    print("  report, meta = generate_multimodal_report(adata, features)")
    print("\nRequirements:")
    print("  - transformers>=4.45.0")
    print("  - PIL (pillow)")
    print("  - torch>=2.4.0")
    print("  - MedGemma 1.5 model access on HuggingFace")
    print("\nFor testing:")
    print("  python scripts/test_multimodal_report.py --h5ad <path>")

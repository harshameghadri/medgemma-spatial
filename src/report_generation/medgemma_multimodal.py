"""
MedGemma report generation (text-only or H&E image + spatial features).

Loads ``google/medgemma-4b-it`` (a Gemma 3 image-text-to-text checkpoint) with
its processor and chat template. Device and dtype are chosen for the hardware:
CUDA uses bfloat16 (or 4-bit NF4 when ``quantize_4bit=True`` and bitsandbytes
is installed); Apple MPS and CPU use bfloat16 unquantised (~8.6 GB of weights).
An optional LoRA adapter is applied when ``MEDGEMMA_ADAPTER`` (or ``adapter_id``)
names a PEFT repo or directory.
"""

import functools
import os
from typing import Dict, Optional, Tuple

import numpy as np
from PIL import Image

DEFAULT_MODEL_ID = os.environ.get("MEDGEMMA_MODEL_ID", "google/medgemma-4b-it")
SYSTEM_PROMPT = (
    "You are an expert computational pathologist. Interpret spatial transcriptomics "
    "summaries for research purposes. Do not invent citations; if you mention prior "
    "work, describe it generically unless you are certain of the reference."
)


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

    if np.issubdtype(img_array.dtype, np.floating):
        img_array = (np.clip(img_array, 0, 1) * 255).astype(np.uint8)

    if img_array.ndim == 2:
        img_array = np.stack([img_array] * 3, axis=-1)

    image = Image.fromarray(img_array).convert('RGB')
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
    """Text-only report prompt (the one the LoRA adapter trains on) plus an instruction to use the H&E image."""
    from src.report_generation.prompt_builder import generate_medgemma_prompt

    prompt = generate_medgemma_prompt(features)
    return prompt.replace(
        "SPATIAL TRANSCRIPTOMICS DATA SUMMARY:",
        "H&E IMAGE: the tissue section above. Describe only morphology you can see in it; "
        "if it is unclear, say so rather than guess.\n\nSPATIAL TRANSCRIPTOMICS DATA SUMMARY:", 1)


def pick_device(device: Optional[str] = None) -> str:
    """Resolve 'auto'/None to cuda, mps or cpu."""
    import torch

    if device and device != "auto":
        return device
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


@functools.lru_cache(maxsize=2)
def load_medgemma(model_id: str = DEFAULT_MODEL_ID, device: Optional[str] = None,
                  adapter_id: Optional[str] = None, quantize_4bit: bool = False):
    """Load (processor, model) once per process; cached across calls."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor

    device = pick_device(device)
    token = os.environ.get("HF_TOKEN")
    kwargs = {"dtype": torch.bfloat16, "token": token, "low_cpu_mem_usage": True}

    if quantize_4bit and device == "cuda":
        from transformers import BitsAndBytesConfig

        kwargs["quantization_config"] = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_compute_dtype=torch.bfloat16,
        )
        kwargs["device_map"] = "auto"

    processor = AutoProcessor.from_pretrained(model_id, token=token)
    model = AutoModelForImageTextToText.from_pretrained(model_id, **kwargs)
    if "device_map" not in kwargs:
        model = model.to(device)

    adapter_id = adapter_id or os.environ.get("MEDGEMMA_ADAPTER") or None
    if adapter_id:
        from peft import PeftModel

        model = PeftModel.from_pretrained(model, adapter_id, token=token)
    model.eval()
    return processor, model


def generate_from_messages(messages, model_id: str = DEFAULT_MODEL_ID, device: Optional[str] = None,
                           max_new_tokens: int = 400, adapter_id: Optional[str] = None) -> str:
    """Apply the chat template, generate greedily and decode only the new tokens."""
    import torch

    processor, model = load_medgemma(model_id, device, adapter_id)
    inputs = processor.apply_chat_template(
        messages, add_generation_prompt=True, tokenize=True,
        return_dict=True, return_tensors="pt",
    ).to(model.device, dtype=model.dtype)
    input_len = inputs["input_ids"].shape[-1]
    with torch.inference_mode():
        output = model.generate(**inputs, max_new_tokens=max_new_tokens, do_sample=False)
    return processor.decode(output[0][input_len:], skip_special_tokens=True).strip()


def _system_message():
    return {"role": "system", "content": [{"type": "text", "text": SYSTEM_PROMPT}]}


def generate_multimodal_report(
    adata,
    features: Dict,
    model_id: str = DEFAULT_MODEL_ID,
    device: Optional[str] = None,
    max_new_tokens: int = 400,
    library_id: Optional[str] = None,
    adapter_id: Optional[str] = None,
) -> Tuple[str, Dict]:
    """Report from the H&E image in adata.uns['spatial'] plus spatial features."""
    image = preprocess_for_medgemma(load_hires_image_from_adata(adata, library_id=library_id))
    prompt = create_multimodal_prompt(features)
    messages = [
        _system_message(),
        {"role": "user", "content": [
            {"type": "image", "image": image},
            {"type": "text", "text": prompt},
        ]},
    ]
    report = generate_from_messages(messages, model_id, device, max_new_tokens, adapter_id)
    metadata = {
        "mode": "multimodal",
        "model_id": model_id,
        "adapter_id": adapter_id or os.environ.get("MEDGEMMA_ADAPTER"),
        "device": pick_device(device),
        "image_size": list(image.size),
        "report_length": len(report),
        "prompt": prompt,
    }
    return report, metadata


def generate_textonly_fallback(
    features: Dict,
    model_id: str = DEFAULT_MODEL_ID,
    device: Optional[str] = None,
    max_new_tokens: int = 400,
    adapter_id: Optional[str] = None,
) -> Tuple[str, Dict]:
    """Report from spatial features only, using the same prompt the LoRA adapter trains on."""
    from src.report_generation.prompt_builder import generate_medgemma_prompt

    prompt = generate_medgemma_prompt(features)
    messages = [
        _system_message(),
        {"role": "user", "content": [{"type": "text", "text": prompt}]},
    ]
    report = generate_from_messages(messages, model_id, device, max_new_tokens, adapter_id)
    metadata = {
        "mode": "text_only",
        "model_id": model_id,
        "adapter_id": adapter_id or os.environ.get("MEDGEMMA_ADAPTER"),
        "device": pick_device(device),
        "report_length": len(report),
        "prompt": prompt,
    }
    return report, metadata

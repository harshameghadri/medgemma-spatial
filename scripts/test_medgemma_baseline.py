#!/usr/bin/env python3
"""
Test pre-trained MedGemma on spatial transcriptomics (baseline).

This establishes whether MedGemma knows anything about spatial
transcriptomics, which determines if our fine-tuning is a "novel task".
"""

import json


def create_spatial_prompt():
    """Create a prompt with spatial transcriptomics features."""

    features = {
        "tissue_type": "colon",
        "cell_types": {
            "Epithelial cells": 19410,
            "Endothelial cells": 6490,
            "NK cells": 4380,
            "Smooth muscle cells": 3960,
            "Macrophages": 2740
        },
        "spatial_metrics": {
            "morans_i": 0.65,
            "spatial_entropy": 1.82
        }
    }

    prompt = f"""You are a computational pathologist. Analyze the following Visium HD spatial transcriptomics data:

SPATIAL TRANSCRIPTOMICS FEATURES:
{json.dumps(features, indent=2)}

Generate a histopathology report that correlates the molecular findings with expected tissue morphology on H&E staining.

Report:"""

    return prompt, features


def test_pretrained_medgemma():
    """Test if pre-trained MedGemma understands spatial transcriptomics."""

    print("="*80)
    print("MEDGEMMA BASELINE TEST - Spatial Transcriptomics Understanding")
    print("="*80)

    prompt, features = create_spatial_prompt()

    print("\n[1/3] Loading MedGemma 4B (instruction-tuned)...")

    try:
        from transformers import AutoTokenizer, AutoModelForCausalLM
        import torch

        model_id = "google/medgemma-4b-it"

        print(f"  Model: {model_id}")
        print(f"  Loading...")

        tokenizer = AutoTokenizer.from_pretrained(model_id)
        model = AutoModelForCausalLM.from_pretrained(
            model_id,
            torch_dtype=torch.bfloat16,
            device_map="auto",
            low_cpu_mem_usage=True
        )

        print(f"  ✓ Model loaded successfully")

    except Exception as e:
        print(f"  ✗ Failed to load: {e}")
        print("\n" + "="*80)
        print("RECOMMENDATION:")
        print("="*80)
        print("\n1. Install HuggingFace CLI:")
        print("   pip install --upgrade huggingface_hub")
        print("\n2. Login to HuggingFace:")
        print("   huggingface-cli login")
        print("\n3. Accept MedGemma license:")
        print("   https://huggingface.co/google/medgemma-4b-it")
        print("\n4. Re-run this script")
        return

    print("\n[2/3] Generating report from spatial features...")

    try:
        inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

        print(f"  Prompt length: {len(prompt)} chars")
        print(f"  Generating (max 300 tokens)...")

        outputs = model.generate(
            **inputs,
            max_new_tokens=300,
            do_sample=False,
            pad_token_id=tokenizer.eos_token_id
        )

        full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
        report = full_output[len(prompt):].strip()

        print(f"  ✓ Generated {len(report)} characters")

    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        import traceback
        traceback.print_exc()
        return

    print("\n[3/3] Analyzing baseline performance...")

    # Check if report shows understanding of spatial transcriptomics
    understands_spatial = any(term in report.lower() for term in [
        'spatial', 'moran', 'autocorrelation', 'spatial organization',
        'spatial pattern', 'spatial distribution', 'clustering'
    ])

    understands_correlation = any(term in report.lower() for term in [
        'h&e', 'histology', 'morphology', 'staining', 'microscopy',
        'tissue architecture', 'cellular arrangement'
    ])

    mentions_cell_types = any(ct.lower() in report.lower()
                              for ct in features['cell_types'].keys())

    mentions_metrics = any(str(round(val, 2)) in report or str(round(val, 1)) in report
                          for val in features['spatial_metrics'].values())

    # Classification
    if understands_spatial and understands_correlation:
        level = "ADVANCED"
        color = "✓"
    elif understands_correlation:
        level = "INTERMEDIATE"
        color = "⚠️"
    else:
        level = "NOVICE"
        color = "✗"

    print(f"\n{'='*80}")
    print(f"BASELINE REPORT")
    print(f"{'='*80}\n")
    print(report)
    print(f"\n{'='*80}")
    print(f"ANALYSIS")
    print(f"{'='*80}")
    print(f"\n  Spatial understanding: {color if understands_spatial else '✗'} "
          f"{'YES' if understands_spatial else 'NO'}")
    print(f"  Morphological correlation: {color if understands_correlation else '✗'} "
          f"{'YES' if understands_correlation else 'NO'}")
    print(f"  Mentions cell types: {'✓' if mentions_cell_types else '✗'} "
          f"{'YES' if mentions_cell_types else 'NO'}")
    print(f"  Cites spatial metrics: {'⚠️' if mentions_metrics else '✗'} "
          f"{'YES (parroting!)' if mentions_metrics else 'NO'}")

    print(f"\n  Overall level: {level}")

    # Determine if novel task
    print(f"\n{'='*80}")
    print(f"NOVEL TASK DETERMINATION")
    print(f"{'='*80}")

    if level == "NOVICE":
        print("\n✅ NOVEL TASK CONFIRMED")
        print("  Pre-trained MedGemma does NOT understand spatial transcriptomics")
        print("  Fine-tuning to add this capability qualifies as novel task")
    elif level == "INTERMEDIATE":
        print("\n⚠️ PARTIAL NOVELTY")
        print("  MedGemma has some histology knowledge but weak on spatial")
        print("  Fine-tuning to improve correlation may still qualify")
    else:
        print("\n❌ NOT NOVEL")
        print("  MedGemma already understands spatial transcriptomics")
        print("  Need to find a different novel task")

    print(f"\n{'='*80}")

    # Save results
    import os
    os.makedirs("outputs/baseline_test", exist_ok=True)

    with open("outputs/baseline_test/medgemma_baseline_report.txt", "w") as f:
        f.write(report)

    with open("outputs/baseline_test/medgemma_baseline_analysis.json", "w") as f:
        json.dump({
            "understands_spatial": understands_spatial,
            "understands_correlation": understands_correlation,
            "mentions_cell_types": mentions_cell_types,
            "mentions_metrics": mentions_metrics,
            "level": level,
            "novel_task": level == "NOVICE",
            "report_length": len(report)
        }, f, indent=2)

    print(f"✓ Results saved to: outputs/baseline_test/")

    return {
        "report": report,
        "level": level,
        "novel_task": level == "NOVICE"
    }


if __name__ == "__main__":
    test_pretrained_medgemma()

#!/usr/bin/env python3
"""
Test pre-trained MedGemma on spatial transcriptomics (RIGOROUS VERSION).

This version does NOT expose raw data - tests if model can interpret
spatial patterns from abstract descriptions only.
"""

import json


def create_abstract_spatial_prompt():
    """Create prompt WITHOUT raw numbers - tests true understanding."""

    prompt = """You are a computational pathologist analyzing Visium HD spatial transcriptomics data from human colon tissue.

AVAILABLE INFORMATION:
- Tissue composition shows 2 major cell populations (>10% each)
- Immune infiltration detected with 5 immune cell types present
- Spatial organization: heterogeneous (high spatial autocorrelation)
- Prediction uncertainty: moderate (15% spots with high entropy)

TASK:
Generate a histopathology report that:
1. Interprets what this spatial heterogeneity means biologically
2. Correlates molecular findings with expected H&E morphology
3. Explains clinical significance

DO NOT list cell types or percentages. Provide biological interpretation only.

Report:"""

    return prompt


def create_baseline_prompt_with_data():
    """Original prompt WITH raw data - allows parroting."""

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


def test_medgemma_understanding():
    """Test both abstract and data-rich prompts."""

    print("="*80)
    print("MEDGEMMA BASELINE TEST V2 - Rigorous Spatial Understanding")
    print("="*80)

    print("\n[1/4] Loading MedGemma 4B (instruction-tuned)...")

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

    # TEST 1: Abstract prompt (no raw data)
    print("\n" + "="*80)
    print("[2/4] TEST 1: Abstract Spatial Prompt (No Raw Data)")
    print("="*80)

    abstract_prompt = create_abstract_spatial_prompt()

    print(f"\n  Prompt length: {len(abstract_prompt)} chars")
    print(f"  Exposes raw numbers: NO")
    print(f"  Generating...")

    try:
        inputs = tokenizer(abstract_prompt, return_tensors="pt").to(model.device)

        outputs = model.generate(
            **inputs,
            max_new_tokens=300,
            do_sample=False,
            pad_token_id=tokenizer.eos_token_id
        )

        full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
        abstract_report = full_output[len(abstract_prompt):].strip()

        print(f"  ✓ Generated {len(abstract_report)} characters")

    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        abstract_report = ""

    # TEST 2: Baseline prompt (with raw data)
    print("\n" + "="*80)
    print("[3/4] TEST 2: Baseline Prompt (With Raw Data)")
    print("="*80)

    baseline_prompt, features = create_baseline_prompt_with_data()

    print(f"\n  Prompt length: {len(baseline_prompt)} chars")
    print(f"  Exposes raw numbers: YES")
    print(f"  Generating...")

    try:
        inputs = tokenizer(baseline_prompt, return_tensors="pt").to(model.device)

        outputs = model.generate(
            **inputs,
            max_new_tokens=300,
            do_sample=False,
            pad_token_id=tokenizer.eos_token_id
        )

        full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)
        baseline_report = full_output[len(baseline_prompt):].strip()

        print(f"  ✓ Generated {len(baseline_report)} characters")

    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        baseline_report = ""

    # ANALYSIS
    print("\n" + "="*80)
    print("[4/4] Comparative Analysis")
    print("="*80)

    # Analyze abstract report (no data to parrot)
    abstract_lower = abstract_report.lower()
    abstract_understands_spatial = any(term in abstract_lower for term in [
        'spatial', 'heterogeneity', 'organization', 'pattern', 'distribution',
        'autocorrelation', 'clustering', 'neighborhood'
    ])
    abstract_interprets = any(term in abstract_lower for term in [
        'suggests', 'indicates', 'consistent with', 'reflects', 'implies',
        'characteristic', 'typical', 'associated'
    ])
    abstract_mentions_morphology = any(term in abstract_lower for term in [
        'h&e', 'histology', 'morphology', 'staining', 'architecture',
        'microscopy', 'cellular arrangement', 'tissue structure'
    ])

    # Analyze baseline report (can parrot)
    baseline_lower = baseline_report.lower()
    baseline_parrots = any(str(count) in baseline_report for count in features['cell_types'].values())
    baseline_mentions_cell_types = sum(1 for ct in features['cell_types'].keys()
                                       if ct.lower() in baseline_lower) >= 3

    print("\n  TEST 1: Abstract Prompt (True Understanding Required)")
    print(f"    Understands spatial concepts: {'✓ YES' if abstract_understands_spatial else '✗ NO'}")
    print(f"    Provides interpretation: {'✓ YES' if abstract_interprets else '✗ NO'}")
    print(f"    Mentions morphology/H&E: {'✓ YES' if abstract_mentions_morphology else '✗ NO'}")

    print("\n  TEST 2: Baseline Prompt (Parroting Possible)")
    print(f"    Parrots exact counts: {'⚠️  YES' if baseline_parrots else '✓ NO'}")
    print(f"    Lists cell types verbatim: {'⚠️  YES' if baseline_mentions_cell_types else '✓ NO'}")

    # VERDICT
    print("\n" + "="*80)
    print("NOVEL TASK DETERMINATION")
    print("="*80)

    if abstract_understands_spatial and abstract_interprets and abstract_mentions_morphology:
        verdict = "NOT NOVEL"
        print("\n❌ NOT NOVEL")
        print("  MedGemma understands spatial transcriptomics WITHOUT seeing data")
        print("  It can interpret abstract spatial patterns → morphology")
        print("  Fine-tuning on this task would NOT qualify as novel")
        novel_task = False
    elif baseline_parrots or baseline_mentions_cell_types:
        verdict = "PARROTING (NOT UNDERSTANDING)"
        print("\n⚠️  PARROTING, NOT UNDERSTANDING")
        print("  MedGemma appears to understand when given raw data")
        print("  But this is just regurgitating input, not interpreting")
        print("  RECOMMENDATION: Use abstract prompts to force interpretation")
        novel_task = False
    else:
        verdict = "NOVEL TASK CONFIRMED"
        print("\n✅ NOVEL TASK CONFIRMED")
        print("  MedGemma does NOT understand spatial transcriptomics")
        print("  Fine-tuning to add this capability qualifies as novel")
        novel_task = True

    print("\n" + "="*80)
    print("REPORTS")
    print("="*80)

    print("\n  TEST 1 REPORT (Abstract):")
    print("  " + "-"*76)
    preview = abstract_report[:300] + "..." if len(abstract_report) > 300 else abstract_report
    for line in preview.split('\n'):
        print(f"  {line}")

    print("\n  TEST 2 REPORT (Baseline):")
    print("  " + "-"*76)
    preview = baseline_report[:300] + "..." if len(baseline_report) > 300 else baseline_report
    for line in preview.split('\n'):
        print(f"  {line}")

    # Save results
    import os
    os.makedirs("outputs/baseline_test_v2", exist_ok=True)

    with open("outputs/baseline_test_v2/abstract_report.txt", "w") as f:
        f.write(abstract_report)

    with open("outputs/baseline_test_v2/baseline_report.txt", "w") as f:
        f.write(baseline_report)

    with open("outputs/baseline_test_v2/analysis.json", "w") as f:
        json.dump({
            "abstract": {
                "understands_spatial": abstract_understands_spatial,
                "interprets": abstract_interprets,
                "mentions_morphology": abstract_mentions_morphology,
                "report_length": len(abstract_report)
            },
            "baseline": {
                "parrots_counts": baseline_parrots,
                "lists_cell_types": baseline_mentions_cell_types,
                "report_length": len(baseline_report)
            },
            "verdict": verdict,
            "novel_task": novel_task
        }, f, indent=2)

    print(f"\n  ✓ Results saved to: outputs/baseline_test_v2/")

    print("\n" + "="*80)

    return {
        "verdict": verdict,
        "novel_task": novel_task,
        "abstract_report": abstract_report,
        "baseline_report": baseline_report
    }


if __name__ == "__main__":
    test_medgemma_understanding()

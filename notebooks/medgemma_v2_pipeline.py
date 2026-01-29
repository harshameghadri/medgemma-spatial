"""
MedGemma V2: Anti-Parrot Uncertainty-Aware Clinical Report Generator

Architecture:
- Stage 0: Annotation quality assessment (doublet detection)
- Stage 1: Uncertainty-aware spatial pattern extraction
- Stage 2: Comparative analysis (vs reference cohorts)
- Stage 3: Conditional IF-THEN reasoning
- Stage 4: Hypothesis generation (testable claims)
- Stage 5: Verification (Llama-3.1 or Qwen2.5)
- Stage 6: Self-audit (programmatic novelty validation)
"""

import torch
import os
import json
from transformers import AutoTokenizer, AutoModelForCausalLM
from dotenv import load_dotenv
import numpy as np
from typing import Dict, Optional
import subprocess
import sys

# Import self-audit module
sys.path.append(os.path.dirname(__file__))
from medgemma_self_audit import run_self_audit


# ==============================================================================
# DEVICE DETECTION (UPGRADE #10: CPU/GPU ARCHITECTURE)
# ==============================================================================

def get_optimal_device():
    """
    Smart device selection for cross-platform compatibility (Kaggle/Local).
    """
    if torch.cuda.is_available():
        device = "cuda"
        dtype = torch.float16
        print(f"Using CUDA GPU: {torch.cuda.get_device_name(0)}")
    elif torch.backends.mps.is_available():
        print("WARNING: MPS backend has generation bugs in PyTorch 2.10+")
        print("Falling back to CPU for stability")
        device = "cpu"
        dtype = torch.float32
    else:
        device = "cpu"
        dtype = torch.float32
        print("No GPU detected, using CPU")

    return device, dtype


# ==============================================================================
# STAGE 2: COMPARATIVE ANALYSIS (FORCED COMPARISON)
# ==============================================================================

def load_reference_cohort(tissue_type='breast_cancer'):
    """
    Load reference phenotypes for comparative analysis.

    Note: In production, this would load from a database.
    For now, using synthetic reference ranges.
    """
    reference_cohorts = {
        'breast_cancer': {
            'hot_tumor': {
                'tumor_immune_interface_pct': (40, 60),
                'spatial_entropy': (0.6, 0.9),
                'immune_markers': ['CD8A', 'GZMA', 'PRF1'],
                'description': 'High immune infiltration, active T-cell response'
            },
            'cold_tumor': {
                'tumor_immune_interface_pct': (5, 20),
                'spatial_entropy': (0.1, 0.3),
                'immune_markers': ['CD68', 'CD163'],
                'description': 'Immune excluded, macrophage-dominated'
            },
            'hybrid_immune_inflamed': {
                'tumor_immune_interface_pct': (25, 40),
                'spatial_entropy': (0.4, 0.6),
                'immune_markers': ['CD79A', 'IGHG1', 'CD68'],
                'description': 'Mixed B-cell and macrophage infiltration, intermediate response'
            }
        }
    }

    return reference_cohorts.get(tissue_type, {})


def generate_comparative_prompt(uncertainty_features: Dict, reference_cohort: Dict) -> str:
    """
    Generate prompt that FORCES comparative analysis.
    """
    # Extract key metrics
    stage1 = uncertainty_features.get('stage_1_spatial_patterns', {})
    morans = stage1.get('morans_i', {})
    entropy = stage1.get('spatial_entropy', {})

    # Get signal strength distribution with gene names
    if 'signal_strength' in morans:
        strong_genes = [g for g, s in zip(morans.get('genes', []), morans['signal_strength']) if s == 'STRONG']
        moderate_genes = [g for g, s in zip(morans.get('genes', []), morans['signal_strength']) if s == 'MODERATE']
        weak_genes = [g for g, s in zip(morans.get('genes', []), morans['signal_strength']) if s == 'WEAK']
        strong_count = len(strong_genes)
        moderate_count = len(moderate_genes)
        weak_count = len(weak_genes)
    else:
        strong_genes, moderate_genes, weak_genes = [], [], []
        strong_count, moderate_count, weak_count = 0, 0, 0

    prompt = f"""You are a spatial pathology expert analyzing spatial transcriptomics data.

CRITICAL INSTRUCTIONS:
1. Do NOT repeat numeric values from the data
2. Do NOT use generic phrases like "spatial segregation observed"
3. MUST compare this sample to reference phenotypes
4. MUST report uncertainty (p-values, confidence intervals)
5. MUST acknowledge weak signals explicitly

=== SPATIAL SIGNAL QUALITY ===
- STRONG spatial signals: {strong_count} genes (p < 0.001, |I| > 0.3)
  Top genes: {', '.join(strong_genes[:10]) if strong_genes else 'None'}
- MODERATE signals: {moderate_count} genes (p < 0.05, |I| > 0.1)
  Top genes: {', '.join(moderate_genes[:10]) if moderate_genes else 'None'}
- WEAK/NONE signals: {weak_count} genes

Spatial entropy: {entropy.get('mean_entropy', 'N/A'):.3f} (95% CI: [{entropy.get('ci_lower', 'N/A'):.3f}, {entropy.get('ci_upper', 'N/A'):.3f}])

CRITICAL: You MUST cite these exact counts in your report:
- {strong_count} STRONG signals
- {moderate_count} MODERATE signals

=== REFERENCE PHENOTYPES (for comparison) ===
"""

    # Add reference phenotypes
    for phenotype_name, phenotype_data in reference_cohort.items():
        prompt += f"\n{phenotype_name.upper().replace('_', ' ')}:\n"
        prompt += f"  - Interface: {phenotype_data['tumor_immune_interface_pct'][0]}-{phenotype_data['tumor_immune_interface_pct'][1]}%\n"
        prompt += f"  - Entropy: {phenotype_data['spatial_entropy'][0]:.2f}-{phenotype_data['spatial_entropy'][1]:.2f}\n"
        prompt += f"  - Key markers: {', '.join(phenotype_data['immune_markers'])}\n"
        prompt += f"  - Pattern: {phenotype_data['description']}\n"

    prompt += f"""

=== YOUR TASK ===
Answer these specific questions with EVIDENCE-BASED comparisons:

1. PHENOTYPE CLASSIFICATION (with uncertainty):
   - Which reference phenotype does this sample MOST resemble?
   - What is your confidence level (LOW/MODERATE/HIGH)?
   - What specific metrics support this classification?
   - Does it show hybrid features? If so, which?

2. SIGNAL STRENGTH ASSESSMENT:
   - Given {strong_count} STRONG and {moderate_count} MODERATE signals, is spatial analysis reliable?
   - Should we proceed with deeper mechanistic inference, or STOP here?
   - What additional validation would increase confidence?

3. STATISTICAL UNCERTAINTY:
   - Report p-values for top spatial patterns
   - Mention confidence intervals where available
   - Flag any metrics with high uncertainty

4. COMPARATIVE INSIGHTS (NOT generic descriptions):
   - How does this sample DIFFER from typical {list(reference_cohort.keys())[0].replace('_', ' ')}?
   - Are there unexpected spatial features given the phenotype?
   - What spatial patterns are ABSENT that you'd expect?

5. TESTABLE HYPOTHESES:
   - Propose ONE specific hypothesis that could be tested with additional data
   - What experiment or analysis would validate/refute this?

FORMAT RULES:
- Use "p < 0.001" notation, not "significant"
- Say "entropy = X (95% CI: [a, b])" not just "entropy = X"
- Compare to reference ranges: "Interface at Y% (typical cold tumor: 5-20%)"
- If signal is WEAK, say "WEAK SIGNAL - further analysis not warranted"
- Maximum 250 words, every sentence must cite data or comparisons

BEGIN ANALYSIS:
"""

    return prompt


# ==============================================================================
# STAGE 3: CONDITIONAL IF-THEN REASONING
# ==============================================================================

def generate_conditional_prompt(uncertainty_features: Dict, stage2_output: str) -> str:
    """
    Force conditional reasoning: IF X THEN Y, with explicit logic.
    """
    stage_stopping = uncertainty_features.get('stage_stopping', {})
    decision = stage_stopping.get('decision', 'PROCEED')

    prompt = f"""You are reviewing a spatial pathology preliminary report.

PREVIOUS ANALYSIS:
{stage2_output}

SIGNAL QUALITY DECISION: {decision}

=== CONDITIONAL REASONING TASK ===

Apply IF-THEN logic to generate mechanistic hypotheses:

IF [condition from data] THEN [biological implication] BECAUSE [mechanism]

RULES:
1. Conditions must cite specific metrics (e.g., "IF entropy < 0.3")
2. Implications must be testable
3. Mechanisms must reference known biology (NOT generic pathways)
4. Use format: "IF ... (p < 0.05) THEN ... BECAUSE ..."

Example GOOD reasoning:
"IF spatial entropy = 0.25 (< 0.3, p < 0.001) THEN tumor is spatially homogeneous,
 BECAUSE low entropy indicates lack of microenvironment diversity, consistent with
 clonal expansion without immune infiltration."

Example BAD reasoning (FORBIDDEN):
"Spatial patterns suggest tumor microenvironment heterogeneity." (NO LOGIC, NO CONDITION)

Generate 2-3 conditional statements for this sample.
Maximum 150 words.
"""

    return prompt


# ==============================================================================
# MEDGEMMA V2 PIPELINE
# ==============================================================================

def load_medgemma_model(device="cpu", dtype=torch.float32):
    """Load MedGemma with device auto-detection."""
    load_dotenv()

    if not os.environ.get('HF_TOKEN'):
        raise ValueError("HF_TOKEN not found in .env file")

    model_id = "google/medgemma-4b-it"
    print(f"Loading MedGemma on {device} with {dtype}...")

    tokenizer = AutoTokenizer.from_pretrained(
        model_id,
        token=os.environ.get('HF_TOKEN')
    )

    model = AutoModelForCausalLM.from_pretrained(
        model_id,
        torch_dtype=dtype,
        device_map=device,
        trust_remote_code=True,
        token=os.environ.get('HF_TOKEN'),
        low_cpu_mem_usage=True
    )

    return model, tokenizer


def generate_report_stage(model, tokenizer, prompt: str, max_tokens=512) -> str:
    """Generate report for one stage."""
    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_tokens,
        temperature=0.7,
        top_p=0.9,
        do_sample=True
    )

    response = tokenizer.decode(outputs[0], skip_special_tokens=True)

    # Extract only the generated portion (after prompt)
    if prompt in response:
        response = response.replace(prompt, "").strip()

    return response


def run_medgemma_v2_pipeline(
    uncertainty_features_path: str,
    output_path: str,
    tissue_type: str = 'breast_cancer'
):
    """
    Execute full MedGemma V2 anti-parrot pipeline.
    """
    print("=" * 80)
    print("MEDGEMMA V2: UNCERTAINTY-AWARE CLINICAL REPORT GENERATOR")
    print("=" * 80)

    # Load uncertainty features
    print("\n[1/8] Loading uncertainty-aware spatial features...")
    with open(uncertainty_features_path, 'r') as f:
        uncertainty_features = json.load(f)

    # Check stage stopping decision
    decision = uncertainty_features.get('stage_stopping', {}).get('decision', 'PROCEED')
    print(f"  Stage stopping decision: {decision}")

    if decision == 'STOP_WEAK_SIGNAL':
        print("\n⚠️  WEAK SIGNAL DETECTED - Stopping analysis")
        report = {
            'decision': 'STOPPED',
            'reason': uncertainty_features['stage_stopping']['rationale'],
            'recommendation': 'Signal quality insufficient for clinical interpretation. '
                            'Consider: (1) deeper sequencing, (2) additional samples, '
                            '(3) orthogonal validation (IHC, IF).'
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"\nStopping report saved to: {output_path}")
        return report

    # Load reference cohort
    print("\n[2/8] Loading reference phenotypes...")
    reference_cohort = load_reference_cohort(tissue_type)
    print(f"  Loaded {len(reference_cohort)} reference phenotypes")

    # Load MedGemma
    print("\n[3/8] Loading MedGemma model...")
    device, dtype = get_optimal_device()
    model, tokenizer = load_medgemma_model(device, dtype)
    print("  Model loaded successfully")

    # Stage 2: Comparative analysis
    print("\n[4/8] STAGE 2: Comparative Analysis")
    stage2_prompt = generate_comparative_prompt(uncertainty_features, reference_cohort)
    stage2_output = generate_report_stage(model, tokenizer, stage2_prompt, max_tokens=400)
    print(f"  Generated {len(stage2_output)} characters")

    # Stage 3: Conditional reasoning
    print("\n[5/8] STAGE 3: Conditional IF-THEN Reasoning")
    stage3_prompt = generate_conditional_prompt(uncertainty_features, stage2_output)
    stage3_output = generate_report_stage(model, tokenizer, stage3_prompt, max_tokens=300)
    print(f"  Generated {len(stage3_output)} characters")

    # Combine outputs with de-duplication
    # Split into sentences and remove exact duplicates while preserving order
    def deduplicate_text(text: str) -> str:
        """Remove duplicate sentences while preserving order."""
        sentences = []
        seen = set()
        for line in text.split('\n'):
            line = line.strip()
            if line and line not in seen:
                sentences.append(line)
                seen.add(line)
        return '\n'.join(sentences)

    combined_output = f"{stage2_output}\n\n{stage3_output}"
    combined_output = deduplicate_text(combined_output)

    # Stage 6: Self-audit
    print("\n[6/8] STAGE 6: Self-Audit (Programmatic Validation)")
    audit_results = run_self_audit(
        medgemma_output=combined_output,
        input_features=uncertainty_features,
        uncertainty_data=uncertainty_features
    )

    # Check if audit passed
    if not audit_results['overall']['pass']:
        print("\n❌ SELF-AUDIT FAILED")
        print(f"Critical failures: {', '.join(audit_results['overall']['critical_failures'])}")

        # Optionally: regenerate with stricter prompt
        print("\n[7/8] Regenerating with anti-parrot constraints...")
        retry_prompt = stage2_prompt + "\n\nWARNING: Previous attempt failed audit for: " + \
                      ", ".join(audit_results['overall']['critical_failures']) + \
                      "\nREGENERATE avoiding these issues."

        stage2_output_retry = generate_report_stage(model, tokenizer, retry_prompt, max_tokens=400)
        combined_output = f"{stage2_output_retry}\n\n{stage3_output}"

        # Re-audit
        audit_results = run_self_audit(combined_output, uncertainty_features, uncertainty_features)

    # Compile final report
    final_report = {
        'report_text': combined_output,
        'audit_results': audit_results,
        'metadata': {
            'medgemma_version': 'v2_uncertainty_aware',
            'tissue_type': tissue_type,
            'stage_stopping_decision': decision,
            'device': device,
            'dtype': str(dtype)
        }
    }

    # Save report
    print("\n[8/8] Saving final report...")
    with open(output_path, 'w') as f:
        json.dump(final_report, f, indent=2)

    # Also save human-readable version
    txt_path = output_path.replace('.json', '.txt')
    with open(txt_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MEDGEMMA V2 SPATIAL PATHOLOGY REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(combined_output)
        f.write("\n\n" + "=" * 80 + "\n")
        f.write(f"Audit Status: {'PASS' if audit_results['overall']['pass'] else 'FAIL'}\n")
        if not audit_results['overall']['pass']:
            f.write(f"Failures: {', '.join(audit_results['overall']['critical_failures'])}\n")
        f.write("=" * 80 + "\n")

    print(f"\nReport saved to: {output_path}")
    print(f"Human-readable: {txt_path}")
    print("=" * 80)

    return final_report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="MedGemma V2 Clinical Report Generator")
    parser.add_argument("--features", required=True, help="Path to uncertainty_spatial_features.json")
    parser.add_argument("--output", required=True, help="Output path for clinical report")
    parser.add_argument("--tissue", default="breast_cancer", help="Tissue type for reference cohort")

    args = parser.parse_args()

    run_medgemma_v2_pipeline(
        uncertainty_features_path=args.features,
        output_path=args.output,
        tissue_type=args.tissue
    )

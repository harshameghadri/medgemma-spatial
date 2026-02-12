#!/usr/bin/env python3
"""
MedGemma Report Generator with Anti-Parroting Safeguards.

Generates standardized clinical pathology reports from spatial transcriptomics
features while avoiding regurgitation of input data.
"""

import json
import os
from typing import Dict, List, Optional
import warnings
warnings.filterwarnings('ignore')


def load_features_json(json_path: str) -> Dict:
    """Load features JSON from spatial analysis."""
    with open(json_path, 'r') as f:
        return json.load(f)


def create_anti_parroting_prompt(features: Dict) -> str:
    """
    Create prompt that encourages synthesis, not regurgitation.

    Key strategies:
    1. Don't list raw numbers - ask for interpretation
    2. Request biological significance, not feature enumeration
    3. Ask for clinical implications, not data summary
    4. Require synthesis across multiple features
    5. Demand domain knowledge application
    """

    # Extract key information (but don't expose raw numbers in prompt)
    annotation = features.get('annotation', {})
    # Support both key names from different pipeline versions
    cell_types = (
        annotation.get('cell_type_counts') or
        annotation.get('cluster_distribution') or
        {}
    )
    n_spots = sum(cell_types.values()) if cell_types else annotation.get('n_spots', 1)

    # Identify major cell populations (>10% of spots)
    major_populations = {k: v for k, v in cell_types.items()
                        if v / n_spots > 0.10 and k != 'Unknown'}

    # Identify immune populations
    immune_cells = {k: v for k, v in cell_types.items()
                    if any(term in k.lower() for term in
                          ['t cells', 'b cells', 'nk cells', 'macrophage',
                           'plasma cells', 'dendritic', 'neutrophil'])}

    # Get spatial metrics (but don't list them)
    has_spatial_heterogeneity = features.get('spatial_heterogeneity', {}).get('morans_i_mean', 0) > 0.3
    has_high_uncertainty = features.get('uncertainty', {}).get('mean_prediction_entropy', 0) > 1.5

    # Build interpretive prompt (NO RAW DATA)
    prompt = f"""You are a computational pathologist analyzing Visium HD spatial transcriptomics data from human tissue.

TASK: Generate a concise clinical pathology report (150-200 words) that synthesizes the biological findings.

AVAILABLE INFORMATION:
- Tissue composition shows {len(major_populations)} major cell populations
- Immune infiltration detected with {len(immune_cells)} immune cell types present
- Spatial organization: {'heterogeneous' if has_spatial_heterogeneity else 'relatively uniform'}
- Tissue complexity: {'high diversity' if has_high_uncertainty else 'defined regions'}

REPORTING REQUIREMENTS:
1. Describe the predominant tissue architecture and cellular composition
2. Characterize the immune microenvironment and infiltration patterns
3. Comment on spatial organization and tissue heterogeneity
4. Provide brief clinical or biological interpretation
5. Note any unexpected or significant findings

CRITICAL:
- DO NOT list cell type percentages or raw counts
- DO NOT enumerate spatial metrics or statistics
- SYNTHESIZE biological meaning, don't summarize data
- Use domain knowledge to interpret patterns
- Focus on clinical/biological significance

Generate the report now:"""

    return prompt


def create_medgemma_prompt_with_features(features: Dict) -> str:
    """
    Alternative: Include features but structure to discourage parroting.

    Uses guided questions that require interpretation, not regurgitation.
    """

    annotation = features.get('annotation', {})
    cell_types = annotation.get('cell_type_counts', {})
    spatial = features.get('spatial_heterogeneity', {})
    uncertainty = features.get('uncertainty', {})

    n_spots = sum(cell_types.values())
    major_populations = {k: v for k, v in cell_types.items() if v / n_spots > 0.10}
    immune_types = {k: v for k, v in cell_types.items() if any(x in k.lower() for x in ['t cell', 'b cell', 'macrophage', 'nk cell', 'dendritic'])}

    prompt = f"""You are analyzing Visium HD spatial transcriptomics data from human tissue.

CELLULAR COMPOSITION:
- Total spots analyzed: {n_spots}
- Major cell populations (>10%): {len(major_populations)} types
- Immune cell diversity: {len(immune_types)} distinct immune populations

SPATIAL METRICS:
- Moran's I (autocorrelation): {spatial.get('morans_i_mean', 'N/A')}
- Spatial entropy: {uncertainty.get('mean_prediction_entropy', 'N/A')}

INTERPRETIVE QUESTIONS (answer in narrative form):

1. TISSUE ARCHITECTURE: What does the cellular composition tell us about the tissue structure? Consider epithelial vs stromal vs immune compartments.

2. TUMOR MICROENVIRONMENT: How would you characterize the immune infiltration? What are the implications for anti-tumor immunity?

3. SPATIAL ORGANIZATION: What does the spatial autocorrelation suggest about tissue organization? Are cells clustered or dispersed?

4. BIOLOGICAL SIGNIFICANCE: What is the most notable biological finding? What might this indicate about disease state or treatment response?

5. CLINICAL CONTEXT: How might these findings inform diagnosis, prognosis, or treatment decisions?

Generate a cohesive 150-200 word clinical pathology report that addresses these questions. Focus on biological interpretation and clinical relevance, not data enumeration.

Report:"""

    return prompt


def create_structured_template_prompt(features: Dict) -> str:
    """
    Use structured template to guide interpretation over enumeration.
    """

    annotation = features.get('annotation', {})
    cell_types = annotation.get('cell_type_counts', {})
    n_spots = sum(cell_types.values())

    # Compute derived insights (not raw data)
    epithelial_fraction = sum(v for k, v in cell_types.items()
                              if 'epithelial' in k.lower()) / n_spots
    immune_fraction = sum(v for k, v in cell_types.items()
                          if any(term in k.lower() for term in
                                ['t cells', 'b cells', 'macrophage', 'plasma', 'nk'])) / n_spots
    stromal_fraction = sum(v for k, v in cell_types.items()
                           if any(term in k.lower() for term in
                                 ['fibroblast', 'endothelial', 'smooth muscle'])) / n_spots

    # Classify patterns
    epithelial_status = "intact" if epithelial_fraction > 0.3 else "disrupted"
    immune_status = "high" if immune_fraction > 0.15 else "moderate" if immune_fraction > 0.08 else "low"
    stromal_status = "prominent" if stromal_fraction > 0.20 else "minimal"

    prompt = f"""Generate a clinical pathology report for Visium HD spatial transcriptomics analysis of human tissue.

TISSUE CHARACTERISTICS (interpret these patterns):
- Epithelial architecture: {epithelial_status}
- Immune infiltration: {immune_status}
- Stromal component: {stromal_status}

REPORT TEMPLATE (fill in with biological interpretation):

HISTOLOGICAL COMPOSITION:
[Describe the overall tissue architecture based on epithelial status. What does this suggest about tumor structure?]

IMMUNE MICROENVIRONMENT:
[Characterize the immune infiltration level. What types of immune responses are present? Implications for anti-tumor immunity?]

STROMAL FEATURES:
[Describe the stromal compartment. Is there desmoplasia or vascular remodeling?]

SPATIAL ORGANIZATION:
[Comment on tissue heterogeneity and cellular organization patterns]

CLINICAL SIGNIFICANCE:
[Provide brief interpretation of key findings relevant to diagnosis or prognosis]

Generate the report (150-200 words, continuous narrative format, no bullets):"""

    return prompt


def generate_report_with_medgemma(features_json_path: str,
                                  output_path: str,
                                  prompt_strategy: str = 'anti_parroting',
                                  model_name: str = 'google/medgemma-2b-it') -> Dict:
    """
    Generate clinical report using MedGemma.

    Args:
        features_json_path: Path to spatial features JSON
        output_path: Where to save report
        prompt_strategy: 'anti_parroting', 'guided_questions', or 'structured_template'
        model_name: MedGemma model identifier

    Returns:
        Dict with report text and metadata
    """

    print("="*80)
    print("MEDGEMMA REPORT GENERATION")
    print("="*80)

    # Load features
    print(f"\n[1/5] Loading features from: {features_json_path}")
    features = load_features_json(features_json_path)

    # Select prompt strategy
    print(f"\n[2/5] Creating prompt (strategy: {prompt_strategy})...")
    if prompt_strategy == 'anti_parroting':
        prompt = create_anti_parroting_prompt(features)
    elif prompt_strategy == 'guided_questions':
        prompt = create_medgemma_prompt_with_features(features)
    elif prompt_strategy == 'structured_template':
        prompt = create_structured_template_prompt(features)
    else:
        raise ValueError(f"Unknown prompt strategy: {prompt_strategy}")

    print(f"  Prompt length: {len(prompt)} characters")

    # Load MedGemma model
    print(f"\n[3/5] Loading MedGemma model: {model_name}...")
    try:
        from transformers import AutoTokenizer, AutoModelForCausalLM
        import torch

        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForCausalLM.from_pretrained(
            model_name,
            torch_dtype=torch.float16,
            device_map='auto',
            low_cpu_mem_usage=True
        )
        print(f"  ✓ Model loaded successfully")

    except Exception as e:
        print(f"  ✗ Failed to load model: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'prompt': prompt
        }

    # Generate report
    print(f"\n[4/5] Generating report...")
    try:
        inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

        outputs = model.generate(
            **inputs,
            max_new_tokens=300,
            temperature=0.7,
            top_p=0.9,
            do_sample=True,
            pad_token_id=tokenizer.eos_token_id
        )

        report = tokenizer.decode(outputs[0], skip_special_tokens=True)

        # Extract just the generated part (after prompt)
        report_text = report[len(prompt):].strip()

        print(f"  ✓ Generated {len(report_text)} characters")

    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'prompt': prompt
        }

    # Save report
    print(f"\n[5/5] Saving report to: {output_path}")

    result = {
        'status': 'success',
        'report': report_text,
        'prompt_strategy': prompt_strategy,
        'model': model_name,
        'features_source': features_json_path,
        'prompt': prompt,
        'word_count': len(report_text.split())
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Save as JSON
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)

    # Also save text-only version
    text_path = output_path.replace('.json', '.txt')
    with open(text_path, 'w') as f:
        f.write(report_text)

    print(f"  ✓ Saved to: {output_path}")
    print(f"  ✓ Text version: {text_path}")

    # Analyze report for parroting
    analysis = analyze_report_quality(report_text, features)

    print(f"\n{'='*80}")
    print(f"REPORT QUALITY ANALYSIS")
    print(f"{'='*80}")
    print(f"  Word count: {analysis['word_count']}")
    print(f"  Contains raw numbers: {'⚠️  YES' if analysis['has_raw_numbers'] else '✓ NO'}")
    print(f"  Lists cell types: {'⚠️  YES' if analysis['lists_cell_types'] else '✓ NO'}")
    print(f"  Interpretation present: {'✓ YES' if analysis['has_interpretation'] else '⚠️  NO'}")
    print(f"  Clinical context: {'✓ YES' if analysis['has_clinical_context'] else '⚠️  NO'}")
    print(f"  Parroting risk: {analysis['parroting_risk']}")
    print(f"{'='*80}")

    result['quality_analysis'] = analysis

    return result


def analyze_report_quality(report: str, features: Dict) -> Dict:
    """
    Analyze report for parroting and quality indicators.
    """

    report_lower = report.lower()

    # Check for raw numbers from features
    cell_types = features.get('annotation', {}).get('cell_type_counts', {})
    has_raw_numbers = any(str(count) in report for count in cell_types.values())

    # Check if listing cell types verbatim
    cell_type_mentions = sum(1 for ct in cell_types.keys() if ct.lower() in report_lower)
    lists_cell_types = cell_type_mentions > 3  # More than 3 = likely listing

    # Check for interpretation keywords
    interpretation_terms = ['suggests', 'indicates', 'consistent with', 'characteristic of',
                           'typical of', 'may indicate', 'implies', 'reflects']
    has_interpretation = any(term in report_lower for term in interpretation_terms)

    # Check for clinical context
    clinical_terms = ['prognosis', 'diagnosis', 'treatment', 'clinical', 'therapeutic',
                     'patient', 'outcome', 'response']
    has_clinical_context = any(term in report_lower for term in clinical_terms)

    # Calculate parroting risk
    risk_score = 0
    if has_raw_numbers:
        risk_score += 3
    if lists_cell_types:
        risk_score += 2
    if not has_interpretation:
        risk_score += 2
    if not has_clinical_context:
        risk_score += 1

    if risk_score >= 5:
        parroting_risk = "HIGH"
    elif risk_score >= 3:
        parroting_risk = "MODERATE"
    else:
        parroting_risk = "LOW"

    return {
        'word_count': len(report.split()),
        'has_raw_numbers': has_raw_numbers,
        'lists_cell_types': lists_cell_types,
        'has_interpretation': has_interpretation,
        'has_clinical_context': has_clinical_context,
        'parroting_risk': parroting_risk,
        'risk_score': risk_score
    }


if __name__ == "__main__":
    print("MedGemma Report Generator with Anti-Parroting Safeguards")
    print("Strategies: anti_parroting, guided_questions, structured_template")

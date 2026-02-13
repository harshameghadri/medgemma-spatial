"""
MedGemma prompt construction for spatial transcriptomics reports.

Generates prompts that encourage clinical synthesis rather than
data regurgitation from the spatial analysis pipeline.
"""

import json
import os
from typing import Dict


def load_features_json(json_path: str) -> Dict:
    """Load spatial features JSON from disk."""
    with open(json_path, 'r') as f:
        return json.load(f)


def generate_medgemma_prompt(features: Dict) -> str:
    """
    Build a clinical pathology prompt from spatial features.

    Structures the request so MedGemma synthesizes biological meaning
    rather than listing raw numbers back from the input.
    """
    annotation = features.get('annotation', {})
    cell_types = (
        annotation.get('cell_type_counts') or
        annotation.get('cluster_distribution') or
        {}
    )
    n_spots = sum(cell_types.values()) if cell_types else annotation.get('n_spots', 1)

    major_populations = {k: v for k, v in cell_types.items()
                         if v / max(n_spots, 1) > 0.10 and k != 'Unknown'}

    immune_cells = {k: v for k, v in cell_types.items()
                    if any(term in k.lower() for term in
                           ['t cells', 'b cells', 'nk cells', 'macrophage',
                            'plasma cells', 'dendritic', 'neutrophil'])}

    has_spatial_heterogeneity = features.get('spatial_heterogeneity', {}).get('morans_i_mean', 0) > 0.3
    has_high_uncertainty = features.get('uncertainty', {}).get('mean_prediction_entropy', 0) > 1.5

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


def generate_medgemma_prompt_detailed(features: Dict) -> str:
    """
    Alternative prompt: includes feature context with guided interpretation questions.

    Uses structured questions to direct synthesis over enumeration.
    """
    annotation = features.get('annotation', {})
    cell_types = annotation.get('cell_type_counts', {})
    spatial = features.get('spatial_heterogeneity', {})
    uncertainty = features.get('uncertainty', {})

    n_spots = sum(cell_types.values()) if cell_types else 1
    major_populations = {k: v for k, v in cell_types.items() if v / max(n_spots, 1) > 0.10}
    immune_types = {k: v for k, v in cell_types.items()
                    if any(x in k.lower() for x in ['t cell', 'b cell', 'macrophage', 'nk cell', 'dendritic'])}

    prompt = f"""You are analyzing Visium HD spatial transcriptomics data from human tissue.

CELLULAR COMPOSITION:
- Total spots analyzed: {n_spots}
- Major cell populations (>10%): {len(major_populations)} types
- Immune cell diversity: {len(immune_types)} distinct immune populations

SPATIAL METRICS:
- Moran's I (autocorrelation): {spatial.get('morans_i_mean', 'N/A')}
- Spatial entropy: {uncertainty.get('mean_prediction_entropy', 'N/A')}

INTERPRETIVE QUESTIONS (answer in narrative form):

1. TISSUE ARCHITECTURE: What does the cellular composition tell us about the tissue structure?
2. TUMOR MICROENVIRONMENT: How would you characterize the immune infiltration?
3. SPATIAL ORGANIZATION: What does the spatial autocorrelation suggest about tissue organization?
4. BIOLOGICAL SIGNIFICANCE: What is the most notable biological finding?
5. CLINICAL CONTEXT: How might these findings inform diagnosis, prognosis, or treatment?

Generate a cohesive 150-200 word clinical pathology report addressing these questions.

Report:"""

    return prompt


def evaluate_report_quality(report: str, features: Dict) -> Dict:
    """
    Check a generated report for data regurgitation and clinical quality.

    Returns a dict with quality flags and a risk score.
    """
    report_lower = report.lower()

    cell_types = features.get('annotation', {}).get('cell_type_counts', {})
    has_raw_numbers = any(str(count) in report for count in cell_types.values())

    cell_type_mentions = sum(1 for ct in cell_types.keys() if ct.lower() in report_lower)
    lists_cell_types = cell_type_mentions > 3

    interpretation_terms = ['suggests', 'indicates', 'consistent with', 'characteristic of',
                            'typical of', 'may indicate', 'implies', 'reflects']
    has_interpretation = any(term in report_lower for term in interpretation_terms)

    clinical_terms = ['prognosis', 'diagnosis', 'treatment', 'clinical', 'therapeutic',
                      'patient', 'outcome', 'response']
    has_clinical_context = any(term in report_lower for term in clinical_terms)

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
        'risk_score': risk_score,
    }


def generate_report_with_medgemma(features_json_path: str,
                                   output_path: str,
                                   prompt_strategy: str = 'default',
                                   model_name: str = 'google/medgemma-2b-it') -> Dict:
    """
    Generate a clinical report using MedGemma.

    Args:
        features_json_path: Path to spatial features JSON
        output_path: Where to save the result
        prompt_strategy: 'default' or 'detailed'
        model_name: HuggingFace model ID

    Returns:
        Dict with report text and quality metadata
    """
    print("=" * 80)
    print("MEDGEMMA REPORT GENERATION")
    print("=" * 80)

    features = load_features_json(features_json_path)

    if prompt_strategy == 'detailed':
        prompt = generate_medgemma_prompt_detailed(features)
    else:
        prompt = generate_medgemma_prompt(features)

    print(f"  Prompt length: {len(prompt)} characters")

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
        report_text = report[len(prompt):].strip()

    except Exception as e:
        return {'status': 'error', 'error': str(e), 'prompt': prompt}

    result = {
        'status': 'success',
        'report': report_text,
        'prompt_strategy': prompt_strategy,
        'model': model_name,
        'word_count': len(report_text.split()),
        'quality_analysis': evaluate_report_quality(report_text, features),
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)

    return result

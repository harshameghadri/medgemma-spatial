#!/usr/bin/env python3
"""
Adversarial test for data leakage in MedGemma prompts.

Tests that tissue type, file names, and identifying information are NOT exposed.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from notebooks.medgemma_report_generator import (
    create_anti_parroting_prompt,
    create_medgemma_prompt_with_features,
    create_structured_template_prompt
)


def create_mock_features(tissue_type='colon'):
    """Create mock features simulating real spatial data."""
    return {
        'annotation': {
            'cell_type_counts': {
                'Epithelial cells': 19410,
                'Endothelial cells': 6490,
                'NK cells': 4380,
                'Smooth muscle cells': 3960,
                'Macrophages': 2740,
                'T cells': 1230,
                'B cells': 890,
                'Fibroblasts': 1450
            }
        },
        'spatial_heterogeneity': {
            'morans_i_mean': 0.65,
            'spatial_entropy': 1.82
        },
        'uncertainty': {
            'mean_prediction_entropy': 1.45
        }
    }


def check_leakage(prompt: str, features: dict, tissue_type: str = 'colon') -> dict:
    """
    Check prompt for data leakage.

    Returns dict with leakage indicators and severity score.
    """
    prompt_lower = prompt.lower()

    leakage_found = []

    # Check 1: Tissue type mentioned
    tissue_keywords = ['colon', 'breast', 'lung', 'brain', 'pancreas', 'prostate']
    tissue_mentions = [kw for kw in tissue_keywords if kw in prompt_lower]
    if tissue_mentions:
        leakage_found.append({
            'type': 'TISSUE_TYPE',
            'severity': 'CRITICAL',
            'details': f"Tissue types mentioned: {tissue_mentions}"
        })

    # Check 2: Exact cell counts exposed
    cell_counts = list(features['annotation']['cell_type_counts'].values())
    exact_counts = [str(count) for count in cell_counts if str(count) in prompt]
    if exact_counts:
        leakage_found.append({
            'type': 'RAW_COUNTS',
            'severity': 'HIGH',
            'details': f"Exact counts exposed: {exact_counts}"
        })

    # Check 3: Cancer mentioned without qualifier
    if 'cancer' in prompt_lower and 'human tissue' not in prompt_lower:
        leakage_found.append({
            'type': 'CANCER_SPECIFIC',
            'severity': 'MEDIUM',
            'details': "Cancer mentioned without 'human tissue' generalization"
        })

    # Check 4: Specimen type too specific
    specific_terms = ['biopsy', 'resection', 'needle biopsy', 'core biopsy']
    specific_found = [term for term in specific_terms if term in prompt_lower]
    if specific_found:
        leakage_found.append({
            'type': 'SPECIFIC_SPECIMEN',
            'severity': 'LOW',
            'details': f"Specific specimen types: {specific_found}"
        })

    # Check 5: Cell type names in full (listing many)
    cell_type_names = list(features['annotation']['cell_type_counts'].keys())
    cell_type_mentions = sum(1 for ct in cell_type_names if ct.lower() in prompt_lower)
    if cell_type_mentions > 3:
        leakage_found.append({
            'type': 'CELL_TYPE_ENUMERATION',
            'severity': 'MEDIUM',
            'details': f"Listed {cell_type_mentions}/{len(cell_type_names)} cell types by name"
        })

    # Calculate risk score
    severity_weights = {'CRITICAL': 10, 'HIGH': 5, 'MEDIUM': 2, 'LOW': 1}
    risk_score = sum(severity_weights[item['severity']] for item in leakage_found)

    # Determine risk level
    if risk_score >= 10:
        risk_level = 'CRITICAL'
    elif risk_score >= 5:
        risk_level = 'HIGH'
    elif risk_score >= 2:
        risk_level = 'MEDIUM'
    else:
        risk_level = 'LOW'

    return {
        'risk_level': risk_level,
        'risk_score': risk_score,
        'leakage_count': len(leakage_found),
        'leakage_details': leakage_found,
        'passed': len(leakage_found) == 0
    }


def test_all_prompt_strategies():
    """Test all prompt generation strategies for leakage."""

    print("="*80)
    print("DATA LEAKAGE ADVERSARIAL TEST")
    print("="*80)

    features = create_mock_features(tissue_type='colon')

    strategies = {
        'anti_parroting': create_anti_parroting_prompt,
        'guided_questions': create_medgemma_prompt_with_features,
        'structured_template': create_structured_template_prompt
    }

    results = {}

    for name, func in strategies.items():
        print(f"\n{'='*80}")
        print(f"Testing Strategy: {name}")
        print(f"{'='*80}")

        try:
            prompt = func(features)
            analysis = check_leakage(prompt, features, tissue_type='colon')

            results[name] = {
                'prompt_length': len(prompt),
                'analysis': analysis
            }

            print(f"\n  Risk Level: {analysis['risk_level']}")
            print(f"  Risk Score: {analysis['risk_score']}/10+")
            print(f"  Leakage Count: {analysis['leakage_count']}")
            print(f"  Passed: {'✓ YES' if analysis['passed'] else '✗ NO'}")

            if analysis['leakage_details']:
                print(f"\n  Leakage Details:")
                for leak in analysis['leakage_details']:
                    print(f"    - [{leak['severity']}] {leak['type']}: {leak['details']}")

            # Show prompt preview
            print(f"\n  Prompt Preview (first 300 chars):")
            preview = prompt[:300].replace('\n', '\n    ')
            print(f"    {preview}...")

        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            results[name] = {'error': str(e)}

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    for name, result in results.items():
        if 'error' in result:
            print(f"  {name}: ERROR - {result['error']}")
        else:
            analysis = result['analysis']
            status = '✓ PASSED' if analysis['passed'] else '✗ FAILED'
            print(f"  {name}: {status} (Risk: {analysis['risk_level']}, Score: {analysis['risk_score']})")

    # Recommendation
    print(f"\n{'='*80}")
    print("RECOMMENDATION")
    print(f"{'='*80}")

    safe_strategies = [name for name, result in results.items()
                       if 'analysis' in result and result['analysis']['risk_level'] in ['LOW', 'MEDIUM']]

    if safe_strategies:
        print(f"\n  ✓ Safe strategies: {', '.join(safe_strategies)}")
        print(f"  → Use '{safe_strategies[0]}' as default")
    else:
        print(f"\n  ✗ No safe strategies found!")
        print(f"  → All prompts leak sensitive information")
        print(f"  → CRITICAL: Fix prompts before deployment")

    # Save results
    output_path = Path('outputs/leakage_test_results.json')
    output_path.parent.mkdir(exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n  Results saved to: {output_path}")

    print(f"\n{'='*80}")

    return results


if __name__ == "__main__":
    results = test_all_prompt_strategies()

    # Exit code based on results
    all_passed = all(
        'analysis' in r and r['analysis']['passed']
        for r in results.values()
    )

    sys.exit(0 if all_passed else 1)

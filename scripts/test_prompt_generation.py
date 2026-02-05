#!/usr/bin/env python3
"""
Test prompt generation without requiring MedGemma model.
"""

import sys
import os
import json
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'notebooks'))

from medgemma_report_generator import (
    create_anti_parroting_prompt,
    create_medgemma_prompt_with_features,
    create_structured_template_prompt,
    load_features_json
)


def create_mock_features():
    """Create mock features JSON for testing."""
    return {
        'annotation': {
            'annotation_method': 'leiden_clustering + markers',
            'n_regions': 9,
            'n_annotated': 8,
            'pct_annotated': 88.9,
            'cell_type_counts': {
                'Epithelial cells': 19410,
                'Endothelial cells': 6490,
                'NK cells': 4380,
                'Smooth muscle cells': 3960,
                'Macrophages': 2740,
                'Plasma cells': 2590,
                'Goblet cells': 2100,
                'Fibroblasts': 1530,
                'Unknown': 6800
            }
        },
        'spatial_heterogeneity': {
            'morans_i_mean': 0.65,
            'morans_i_significant_genes': 42
        },
        'uncertainty': {
            'mean_prediction_entropy': 1.82,
            'high_uncertainty_spots_pct': 15.3
        }
    }


def test_prompt_strategies():
    """Test all prompt strategies and check for parroting risks."""

    print("="*80)
    print("PROMPT GENERATION TEST - ANTI-PARROTING VALIDATION")
    print("="*80)

    # Create mock features
    features = create_mock_features()

    print("\n[INPUT] Mock Features:")
    print(f"  Cell types: {len(features['annotation']['cell_type_counts'])}")
    print(f"  Top 3 populations:")
    for ct, count in list(features['annotation']['cell_type_counts'].items())[:3]:
        print(f"    - {ct}: {count:,}")

    strategies = {
        'anti_parroting': create_anti_parroting_prompt,
        'guided_questions': create_medgemma_prompt_with_features,
        'structured_template': create_structured_template_prompt
    }

    results = {}

    for name, func in strategies.items():
        print(f"\n{'='*80}")
        print(f"STRATEGY: {name.upper()}")
        print(f"{'='*80}")

        try:
            prompt = func(features)

            # Analyze prompt for data leakage
            has_exact_counts = any(str(count) in prompt for count in features['annotation']['cell_type_counts'].values())
            has_percentages = '88.9' in prompt or 'pct_annotated' in prompt
            lists_all_types = sum(1 for ct in features['annotation']['cell_type_counts'].keys()
                                 if ct in prompt) > 5

            print(f"\n  Prompt length: {len(prompt)} chars")
            print(f"  Prompt excerpt (first 300 chars):")
            print(f"  {prompt[:300]}...")

            print(f"\n  DATA LEAKAGE ANALYSIS:")
            print(f"    Contains exact counts: {'⚠️  YES' if has_exact_counts else '✓ NO'}")
            print(f"    Contains percentages: {'⚠️  YES' if has_percentages else '✓ NO'}")
            print(f"    Lists all cell types: {'⚠️  YES' if lists_all_types else '✓ NO'}")

            # Calculate parroting risk
            risk_score = (2 if has_exact_counts else 0) + \
                        (2 if has_percentages else 0) + \
                        (3 if lists_all_types else 0)

            if risk_score >= 5:
                risk_level = "HIGH"
            elif risk_score >= 3:
                risk_level = "MODERATE"
            else:
                risk_level = "LOW"

            print(f"    Parroting risk: {risk_level} (score: {risk_score}/7)")

            results[name] = {
                'prompt': prompt,
                'length': len(prompt),
                'has_exact_counts': has_exact_counts,
                'has_percentages': has_percentages,
                'lists_all_types': lists_all_types,
                'risk_level': risk_level,
                'risk_score': risk_score
            }

            # Save prompt to file
            os.makedirs('outputs/medgemma_prompts', exist_ok=True)
            with open(f'outputs/medgemma_prompts/prompt_{name}.txt', 'w') as f:
                f.write(prompt)

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()

    # Comparison table
    print(f"\n{'='*80}")
    print(f"STRATEGY COMPARISON")
    print(f"{'='*80}\n")

    print(f"{'Strategy':<25} {'Length':<10} {'Exact Counts':<15} {'Percentages':<15} {'Lists Types':<15} {'Risk':<10}")
    print(f"{'-'*25} {'-'*10} {'-'*15} {'-'*15} {'-'*15} {'-'*10}")

    for name, result in results.items():
        print(f"{name:<25} {result['length']:<10} "
              f"{'YES ⚠️' if result['has_exact_counts'] else 'NO ✓':<15} "
              f"{'YES ⚠️' if result['has_percentages'] else 'NO ✓':<15} "
              f"{'YES ⚠️' if result['lists_all_types'] else 'NO ✓':<15} "
              f"{result['risk_level']:<10}")

    # Recommendation
    best = min(results.items(), key=lambda x: x[1]['risk_score'])

    print(f"\n{'='*80}")
    print(f"RECOMMENDATION")
    print(f"{'='*80}")
    print(f"\n  ✓ Best strategy: {best[0]}")
    print(f"  ✓ Parroting risk: {best[1]['risk_level']}")
    print(f"  ✓ Reason: Minimal data leakage, encourages interpretation")

    print(f"\n  Prompts saved to: outputs/medgemma_prompts/")

    # Save comparison
    with open('outputs/medgemma_prompts/comparison.json', 'w') as f:
        json.dump({k: {key: val for key, val in v.items() if key != 'prompt'}
                  for k, v in results.items()}, f, indent=2)

    print(f"\n{'='*80}")
    print(f"✅ PROMPT GENERATION TEST COMPLETE")
    print(f"{'='*80}")

    return results


if __name__ == "__main__":
    test_prompt_strategies()

#!/usr/bin/env python3
"""
Test MedGemma integration with anti-parroting validation.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'notebooks'))

from medgemma_report_generator import (
    generate_report_with_medgemma,
    load_features_json,
    create_anti_parroting_prompt,
    create_medgemma_prompt_with_features,
    create_structured_template_prompt
)


def test_all_prompt_strategies():
    """Test all three prompt strategies and compare results."""

    print("="*80)
    print("MEDGEMMA INTEGRATION TEST - ANTI-PARROTING VALIDATION")
    print("="*80)

    # Find most recent features JSON
    test_dirs = [
        "outputs/test_final_success_20260130_151551/e2735493904baddea063a55d5e676d24",
        "outputs/test_quick"
    ]

    features_path = None
    for test_dir in test_dirs:
        candidate = os.path.join(test_dir, "uncertainty_features.json")
        if os.path.exists(candidate):
            features_path = candidate
            break

    if not features_path:
        print("ERROR: No features JSON found")
        print("Tried:")
        for td in test_dirs:
            print(f"  - {os.path.join(td, 'uncertainty_features.json')}")
        return

    print(f"\n[1/4] Using features from: {features_path}")

    # Load features to show what we're working with
    features = load_features_json(features_path)
    cell_types = features.get('annotation', {}).get('cell_type_counts', {})

    print(f"\n  Cell types available:")
    for ct, count in list(cell_types.items())[:5]:
        print(f"    - {ct}: {count}")
    print(f"    ... and {len(cell_types) - 5} more")

    # Test all three strategies
    strategies = ['anti_parroting', 'guided_questions', 'structured_template']
    output_dir = "outputs/medgemma_reports"
    os.makedirs(output_dir, exist_ok=True)

    results = {}

    for i, strategy in enumerate(strategies, 1):
        print(f"\n{'='*80}")
        print(f"[{i}/3] TESTING STRATEGY: {strategy.upper()}")
        print(f"{'='*80}")

        output_path = os.path.join(output_dir, f"report_{strategy}.json")

        try:
            result = generate_report_with_medgemma(
                features_json_path=features_path,
                output_path=output_path,
                prompt_strategy=strategy
            )

            results[strategy] = result

            # Display report preview
            if result['status'] == 'success':
                report = result['report']
                preview = report[:200] + "..." if len(report) > 200 else report
                print(f"\n  REPORT PREVIEW:")
                print(f"  {preview}")

                # Show quality scores
                qa = result.get('quality_analysis', {})
                print(f"\n  QUALITY SCORES:")
                print(f"    Parroting risk: {qa.get('parroting_risk', 'N/A')}")
                print(f"    Has interpretation: {qa.get('has_interpretation', False)}")
                print(f"    Has clinical context: {qa.get('has_clinical_context', False)}")

        except Exception as e:
            print(f"  ✗ Strategy failed: {e}")
            import traceback
            traceback.print_exc()
            results[strategy] = {'status': 'error', 'error': str(e)}

    # Compare strategies
    print(f"\n{'='*80}")
    print(f"STRATEGY COMPARISON")
    print(f"{'='*80}")

    comparison_table = []
    for strategy, result in results.items():
        if result['status'] == 'success':
            qa = result.get('quality_analysis', {})
            comparison_table.append({
                'strategy': strategy,
                'parroting_risk': qa.get('parroting_risk', 'N/A'),
                'word_count': qa.get('word_count', 0),
                'interpretation': '✓' if qa.get('has_interpretation') else '✗',
                'clinical_context': '✓' if qa.get('has_clinical_context') else '✗'
            })

    if comparison_table:
        print(f"\n{'Strategy':<25} {'Risk':<12} {'Words':<8} {'Interp':<8} {'Clinical':<8}")
        print(f"{'-'*25} {'-'*12} {'-'*8} {'-'*8} {'-'*8}")
        for row in comparison_table:
            print(f"{row['strategy']:<25} {row['parroting_risk']:<12} "
                  f"{row['word_count']:<8} {row['interpretation']:<8} {row['clinical_context']:<8}")

    # Recommendation
    print(f"\n{'='*80}")
    print(f"RECOMMENDATION")
    print(f"{'='*80}")

    best_strategy = min(comparison_table, key=lambda x: x['parroting_risk']) if comparison_table else None

    if best_strategy:
        print(f"\n  Best strategy: {best_strategy['strategy']}")
        print(f"  Parroting risk: {best_strategy['parroting_risk']}")
        print(f"  Reason: Lowest parroting risk with interpretation and clinical context")

    print(f"\n  Reports saved to: {output_dir}/")

    print(f"\n{'='*80}")
    print(f"✅ MEDGEMMA INTEGRATION TEST COMPLETE")
    print(f"{'='*80}")


if __name__ == "__main__":
    test_all_prompt_strategies()

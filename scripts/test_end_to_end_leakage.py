#!/usr/bin/env python3
"""
End-to-end leakage test: Generate reports and check for tissue identification.

Simulates real report generation and validates no identifying information leaks.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from notebooks.medgemma_report_generator import create_anti_parroting_prompt


def create_tissue_features(tissue_type: str) -> dict:
    """Create tissue-specific features for testing."""

    base_features = {
        'spatial_heterogeneity': {'morans_i_mean': 0.65},
        'uncertainty': {'mean_prediction_entropy': 1.45}
    }

    # Tissue-specific cell distributions (realistic)
    tissue_profiles = {
        'colon': {
            'Epithelial cells': 19410,
            'Goblet cells': 5200,
            'Enterocytes': 3800,
            'Macrophages': 2740,
            'T cells': 1230,
            'B cells': 890,
            'Fibroblasts': 1450,
            'Endothelial cells': 2300
        },
        'brain': {
            'Neurons': 12000,
            'Astrocytes': 8500,
            'Oligodendrocytes': 6000,
            'Microglia': 3200,
            'Endothelial cells': 1800,
            'Pericytes': 900
        },
        'breast': {
            'Luminal epithelial cells': 15000,
            'Basal epithelial cells': 4000,
            'Fibroblasts': 6500,
            'Adipocytes': 5000,
            'Macrophages': 2500,
            'T cells': 1500,
            'Endothelial cells': 2000
        }
    }

    base_features['annotation'] = {
        'cell_type_counts': tissue_profiles.get(tissue_type, tissue_profiles['colon'])
    }

    return base_features


def check_tissue_identification(prompt: str, expected_tissue: str) -> dict:
    """
    Check if prompt allows identification of tissue type.

    Returns dict with identification risk assessment.
    """
    prompt_lower = prompt.lower()

    # Explicit tissue mentions
    tissue_keywords = {
        'colon': ['colon', 'colorectal', 'intestine', 'bowel', 'enterocyte', 'goblet'],
        'brain': ['brain', 'neuron', 'astrocyte', 'oligodendrocyte', 'neural', 'cerebral'],
        'breast': ['breast', 'mammary', 'luminal', 'ductal', 'lobular']
    }

    explicit_matches = []
    for tissue, keywords in tissue_keywords.items():
        matches = [kw for kw in keywords if kw in prompt_lower]
        if matches:
            explicit_matches.append({
                'tissue': tissue,
                'matched_keywords': matches,
                'is_correct': tissue == expected_tissue
            })

    # Check for cell type names that reveal tissue
    tissue_specific_cells = {
        'colon': ['goblet', 'enterocyte'],
        'brain': ['neuron', 'astrocyte', 'oligodendrocyte', 'microglia'],
        'breast': ['luminal', 'basal epithelial', 'adipocyte']
    }

    cell_type_leaks = []
    for tissue, cells in tissue_specific_cells.items():
        matches = [cell for cell in cells if cell in prompt_lower]
        if matches:
            cell_type_leaks.append({
                'tissue': tissue,
                'revealing_cells': matches,
                'is_correct': tissue == expected_tissue
            })

    identifiable = len(explicit_matches) > 0 or len(cell_type_leaks) > 0

    return {
        'identifiable': identifiable,
        'expected_tissue': expected_tissue,
        'explicit_matches': explicit_matches,
        'cell_type_leaks': cell_type_leaks,
        'passed': not identifiable
    }


def test_multi_tissue():
    """Test report generation for multiple tissue types."""

    print("="*80)
    print("END-TO-END TISSUE IDENTIFICATION TEST")
    print("="*80)
    print("\nTesting if prompts remain tissue-blind across different samples...")

    tissues = ['colon', 'brain', 'breast']
    results = {}

    for tissue in tissues:
        print(f"\n{'='*80}")
        print(f"Testing Tissue: {tissue.upper()}")
        print(f"{'='*80}")

        features = create_tissue_features(tissue)

        print(f"\n  Cell types in sample: {len(features['annotation']['cell_type_counts'])}")

        # Generate prompt
        prompt = create_anti_parroting_prompt(features)

        # Check for identification
        analysis = check_tissue_identification(prompt, tissue)

        results[tissue] = {
            'prompt_length': len(prompt),
            'analysis': analysis
        }

        print(f"\n  Identifiable: {'✗ YES (LEAKED)' if analysis['identifiable'] else '✓ NO (SAFE)'}")
        print(f"  Passed: {'✓ YES' if analysis['passed'] else '✗ NO'}")

        if analysis['explicit_matches']:
            print(f"\n  Explicit tissue mentions:")
            for match in analysis['explicit_matches']:
                status = '✓' if match['is_correct'] else '✗'
                print(f"    {status} {match['tissue']}: {match['matched_keywords']}")

        if analysis['cell_type_leaks']:
            print(f"\n  Cell type reveals tissue:")
            for leak in analysis['cell_type_leaks']:
                status = '✓' if leak['is_correct'] else '✗'
                print(f"    {status} {leak['tissue']}: {leak['revealing_cells']}")

        # Show prompt sample
        print(f"\n  Prompt sample (first 200 chars):")
        preview = prompt[:200].replace('\n', '\n    ')
        print(f"    {preview}...")

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    all_passed = all(r['analysis']['passed'] for r in results.values())

    for tissue, result in results.items():
        status = '✓ PASSED' if result['analysis']['passed'] else '✗ FAILED'
        print(f"  {tissue.capitalize()}: {status}")

    print(f"\n{'='*80}")
    print("VERDICT")
    print(f"{'='*80}")

    if all_passed:
        print(f"\n  ✓ ALL TESTS PASSED")
        print(f"  → Prompts are tissue-blind")
        print(f"  → Model must interpret spatial patterns without prior knowledge")
    else:
        print(f"\n  ✗ SOME TESTS FAILED")
        print(f"  → Tissue-specific information is leaking")
        print(f"  → Model can identify tissue without analyzing data")

    # Save results
    output_path = Path('outputs/end_to_end_leakage_test.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n  Results saved to: {output_path}")
    print(f"\n{'='*80}")

    return all_passed


if __name__ == "__main__":
    passed = test_multi_tissue()
    sys.exit(0 if passed else 1)

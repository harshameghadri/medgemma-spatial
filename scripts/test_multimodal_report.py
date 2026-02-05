#!/usr/bin/env python3
"""
Test script for MedGemma multimodal report generation.

Tests H&E image + spatial features integration.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc
from src.report_generation.medgemma_multimodal import (
    generate_multimodal_report,
    generate_textonly_fallback
)


def extract_features_from_adata(adata):
    """Extract mock features for testing."""
    n_obs = adata.n_obs if hasattr(adata, 'n_obs') else 10000

    features = {
        'annotation': {
            'cell_type_counts': {
                'Epithelial cells': int(n_obs * 0.45),
                'Endothelial cells': int(n_obs * 0.20),
                'Macrophages': int(n_obs * 0.15),
                'T cells': int(n_obs * 0.10),
                'Fibroblasts': int(n_obs * 0.10)
            }
        },
        'spatial_heterogeneity': {
            'morans_i_mean': 0.65
        },
        'uncertainty': {
            'mean_prediction_entropy': 1.2
        }
    }

    return features


def main():
    parser = argparse.ArgumentParser(description="Test multimodal MedGemma")
    parser.add_argument('--h5ad', required=True, help="Path to H5AD file with spatial images")
    parser.add_argument('--use-images', action='store_true', help="Use multimodal (default: text-only)")
    parser.add_argument('--output', default='/tmp/multimodal_test_report.txt', help="Output path")
    parser.add_argument('--device', default='mps', choices=['mps', 'cuda', 'cpu'], help="Device")

    args = parser.parse_args()

    print("="*80)
    print("MEDGEMMA MULTIMODAL TEST")
    print("="*80)

    print(f"\n[1/4] Loading H5AD: {args.h5ad}")
    try:
        adata = sc.read_h5ad(args.h5ad)
        print(f"  ✓ Loaded: {adata.n_obs} spots, {adata.n_vars} genes")
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return 1

    print("\n[2/4] Extracting features...")
    features = extract_features_from_adata(adata)
    print(f"  ✓ Features ready")

    print(f"\n[3/4] Generating report (mode: {'multimodal' if args.use_images else 'text-only'})...")

    try:
        if args.use_images:
            report, metadata = generate_multimodal_report(
                adata,
                features,
                device=args.device
            )
        else:
            report, metadata = generate_textonly_fallback(
                features,
                device=args.device
            )

        print(f"  ✓ Report generated")

    except Exception as e:
        print(f"  ✗ Generation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    print("\n[4/4] Saving report...")
    with open(args.output, 'w') as f:
        f.write("MULTIMODAL MEDGEMMA REPORT (TEST)\n")
        f.write("="*80 + "\n\n")
        f.write(f"Mode: {metadata['mode']}\n")
        f.write(f"Model: {metadata['model_id']}\n")
        f.write(f"Device: {metadata['device']}\n\n")
        f.write("REPORT:\n")
        f.write("-"*80 + "\n")
        f.write(report)
        f.write("\n" + "="*80 + "\n")

    print(f"  ✓ Saved to: {args.output}")

    print("\n" + "="*80)
    print("REPORT PREVIEW:")
    print("="*80)
    print(report[:500] + "..." if len(report) > 500 else report)
    print("="*80)

    print(f"\n✓ TEST COMPLETE")
    print(f"  Full report: {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

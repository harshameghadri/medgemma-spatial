#!/usr/bin/env python3
"""
Test CellTypist annotation on spatial data subset.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'notebooks'))

import scanpy as sc
from uncertainty_spatial_analysis import annotate_cell_types_celltypist

def test_celltypist_on_subset():
    """Test CellTypist on a small subset of spatial data."""
    print("="*80)
    print("CELLTYPIST ANNOTATION TEST")
    print("="*80)

    # Load the h5ad file from our successful test
    h5ad_path = "outputs/test_quick/square_008um.h5ad"

    # If that doesn't exist, use the extracted one
    if not os.path.exists(h5ad_path):
        h5ad_path = "outputs/test_final_success_20260130_151551/e2735493904baddea063a55d5e676d24/extracted/square_008um.h5ad"

    if not os.path.exists(h5ad_path):
        print(f"ERROR: Could not find h5ad file")
        print(f"Looked for: {h5ad_path}")
        return

    print(f"\n[1/4] Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

    # Subsample to 5000 spots for quick testing
    print(f"\n[2/4] Subsampling to 5,000 spots for quick test...")
    if adata.n_obs > 5000:
        sc.pp.subsample(adata, n_obs=5000, random_state=42)
    print(f"  Subsampled: {adata.n_obs:,} spots")

    # Test CellTypist annotation
    print(f"\n[3/4] Running CellTypist annotation...")
    try:
        adata, metrics = annotate_cell_types_celltypist(
            adata,
            model='Immune_All_Low.pkl',  # Good for cancer samples
            majority_voting=True
        )

        print(f"\n[4/4] Results:")
        print(f"  ✓ CellTypist available: {metrics.get('celltypist_available', False)}")

        if metrics.get('celltypist_available'):
            print(f"  ✓ Number of cell types: {metrics['n_cell_types']}")
            print(f"  ✓ Mean confidence: {metrics['mean_confidence']:.3f}")

            print(f"\n  Cell type distribution:")
            if 'cell_type_counts' in metrics:
                for ct, count in list(metrics['cell_type_counts'].items())[:10]:
                    pct = 100 * count / adata.n_obs
                    print(f"    - {ct}: {count:,} ({pct:.1f}%)")

            # Save annotated data
            output_path = "outputs/test_celltypist_annotated.h5ad"
            adata.write_h5ad(output_path)
            print(f"\n  Saved annotated data to: {output_path}")

            print(f"\n{'='*80}")
            print(f"✅ SUCCESS! CellTypist annotation completed")
            print(f"{'='*80}")

        else:
            print(f"\n{'='*80}")
            print(f"❌ CellTypist annotation failed")
            if 'error' in metrics:
                print(f"Error: {metrics['error']}")
            print(f"{'='*80}")

    except Exception as e:
        print(f"\n{'='*80}")
        print(f"❌ ERROR: {e}")
        print(f"{'='*80}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_celltypist_on_subset()

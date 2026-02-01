#!/usr/bin/env python3
"""
Test spatial region annotation on a subsample of Visium data.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'notebooks'))

import scanpy as sc
from uncertainty_spatial_analysis import annotate_spatial_regions

def test_spatial_annotation():
    """Test Leiden clustering-based spatial annotation."""
    print("="*80)
    print("SPATIAL ANNOTATION TEST (Leiden Clustering)")
    print("="*80)

    # Load the h5ad file from successful test
    h5ad_path = "outputs/test_final_success_20260130_151551/e2735493904baddea063a55d5e676d24/extracted/square_008um.h5ad"

    if not os.path.exists(h5ad_path):
        print(f"ERROR: Could not find h5ad file at {h5ad_path}")
        return

    print(f"\n[1/4] Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

    # Subsample to 5000 spots for quick testing
    print(f"\n[2/4] Subsampling to 5,000 spots for quick test...")
    if adata.n_obs > 5000:
        sc.pp.subsample(adata, n_obs=5000, random_state=42)
    print(f"  Subsampled: {adata.n_obs:,} spots")

    # Test spatial annotation
    print(f"\n[3/4] Running spatial region annotation (Leiden clustering)...")
    try:
        adata, metrics = annotate_spatial_regions(adata, resolution=0.5)

        print(f"\n[4/4] Results:")
        print(f"  ✓ Method: {metrics.get('annotation_method', 'unknown')}")
        print(f"  ✓ Number of regions: {metrics['n_regions']}")
        print(f"  ✓ Resolution: {metrics['resolution']}")
        print(f"  ✓ Region sizes - min: {metrics['region_sizes_min']}, max: {metrics['region_sizes_max']}, mean: {metrics['region_sizes_mean']:.1f}")

        print(f"\n  Region distribution:")
        region_counts = adata.obs['spatial_region'].value_counts().sort_index()
        for region, count in region_counts.items():
            pct = 100 * count / adata.n_obs
            print(f"    - Region {region}: {count:,} spots ({pct:.1f}%)")

        # Check cell_type assignment
        print(f"\n  Cell type labels (first 10):")
        print(f"    {adata.obs['cell_type'].value_counts().head(10).to_dict()}")

        # Save annotated data
        output_path = "outputs/test_spatial_annotation.h5ad"
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        adata.write_h5ad(output_path)
        print(f"\n  Saved annotated data to: {output_path}")

        print(f"\n{'='*80}")
        print(f"✅ SUCCESS! Spatial annotation completed")
        print(f"{'='*80}")

    except Exception as e:
        print(f"\n{'='*80}")
        print(f"❌ ERROR: {e}")
        print(f"{'='*80}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_spatial_annotation()

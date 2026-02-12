#!/usr/bin/env python3
"""
Test marker-based cell type annotation on spatial data.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'notebooks'))

import scanpy as sc
from uncertainty_spatial_analysis import annotate_spatial_regions

def test_marker_annotation():
    """Test marker-based annotation on subsample."""
    print("="*80)
    print("MARKER-BASED CELL TYPE ANNOTATION TEST")
    print("="*80)

    # Load h5ad file
    h5ad_path = "outputs/test_final_success_20260130_151551/e2735493904baddea063a55d5e676d24/extracted/square_008um.h5ad"

    if not os.path.exists(h5ad_path):
        print(f"ERROR: Could not find h5ad file at {h5ad_path}")
        return

    print(f"\n[1/5] Loading data from: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

    # Subsample for quick test
    print(f"\n[2/5] Subsampling to 5,000 spots for quick test...")
    if adata.n_obs > 5000:
        sc.pp.subsample(adata, n_obs=5000, random_state=42)
    print(f"  Subsampled: {adata.n_obs:,} spots")

    # Test WITHOUT markers (baseline)
    print(f"\n[3/5] Testing baseline annotation (no markers)...")
    adata_baseline = adata.copy()
    adata_baseline, metrics_baseline = annotate_spatial_regions(
        adata_baseline,
        resolution=0.5,
        use_markers=False
    )

    print(f"\n  Baseline results:")
    print(f"    Method: {metrics_baseline['annotation_method']}")
    print(f"    Regions: {metrics_baseline['n_regions']}")
    baseline_types = adata_baseline.obs['cell_type'].value_counts().head(5)
    print(f"    Top 5 labels: {baseline_types.to_dict()}")

    # Test WITH markers (enhanced)
    print(f"\n[4/5] Testing marker-based annotation...")
    adata_markers, metrics_markers = annotate_spatial_regions(
        adata,
        resolution=0.5,
        use_markers=True,
        marker_file="data/PanglaoDB_markers_27_Mar_2020.tsv",
        tissue="Colon"
    )

    print(f"\n[5/5] Comparison Results:")
    print(f"\n  BASELINE (no markers):")
    print(f"    Method: {metrics_baseline['annotation_method']}")
    for ct, count in baseline_types.items():
        pct = 100 * count / adata_baseline.n_obs
        print(f"      {ct}: {count:,} ({pct:.1f}%)")

    print(f"\n  ENHANCED (with markers):")
    print(f"    Method: {metrics_markers.get('annotation_method', 'unknown')}")
    if 'pct_annotated' in metrics_markers:
        print(f"    Annotated: {metrics_markers['n_annotated']}/{metrics_markers['n_clusters']} "
              f"({metrics_markers['pct_annotated']:.1f}%)")

    marker_types = adata_markers.obs['cell_type'].value_counts().head(10)
    for ct, count in marker_types.items():
        pct = 100 * count / adata_markers.n_obs
        print(f"      {ct}: {count:,} ({pct:.1f}%)")

    # Save results
    output_dir = "outputs/test_marker_annotation"
    os.makedirs(output_dir, exist_ok=True)

    baseline_path = f"{output_dir}/baseline_no_markers.h5ad"
    markers_path = f"{output_dir}/enhanced_with_markers.h5ad"

    adata_baseline.write_h5ad(baseline_path)
    adata_markers.write_h5ad(markers_path)

    print(f"\n  Saved results:")
    print(f"    Baseline: {baseline_path}")
    print(f"    Enhanced: {markers_path}")

    print(f"\n{'='*80}")
    print(f"✅ SUCCESS! Marker annotation test complete")
    print(f"{'='*80}")

    # Summary
    print(f"\n📊 Summary:")
    n_unknown = (adata_markers.obs['cell_type'] == 'Unknown').sum()
    n_annotated = adata_markers.n_obs - n_unknown
    pct_biological = 100 * n_annotated / adata_markers.n_obs

    print(f"  Baseline: All spots labeled as generic regions (Region_0, Region_1, etc.)")
    print(f"  Enhanced: {n_annotated:,}/{adata_markers.n_obs:,} spots "
          f"({pct_biological:.1f}%) have biological labels")
    print(f"  Improvement: Transformed generic clusters into interpretable cell types!")

if __name__ == "__main__":
    test_marker_annotation()

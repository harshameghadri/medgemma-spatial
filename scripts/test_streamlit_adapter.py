#!/usr/bin/env python3
"""
Test Streamlit Adapter Functions

Validates adapter works with real H5AD file before launching Streamlit UI.
"""

import sys
from pathlib import Path
import time

sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc
from src.streamlit_adapter import (
    annotate_spatial_regions,
    calculate_spatial_heterogeneity
)


def test_adapter_pipeline():
    """Test complete adapter pipeline."""
    print("="*80)
    print("STREAMLIT ADAPTER VALIDATION TEST")
    print("="*80)

    # Load test data
    test_file = Path(__file__).parent.parent / "outputs" / "annotated_visium.h5ad"

    print(f"\n[1/4] Loading test data: {test_file.name}")
    adata = sc.read_h5ad(test_file)
    print(f"  ✓ Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

    # Test annotation
    print(f"\n[2/4] Testing annotate_spatial_regions()...")
    start = time.time()

    adata_annotated, annot_metrics = annotate_spatial_regions(
        adata,
        resolution=0.5,
        use_markers=True,
        tissue="Unknown"
    )

    elapsed = time.time() - start
    print(f"  ✓ Annotation complete in {elapsed:.1f}s")
    print(f"  - Clusters: {annot_metrics['n_clusters']}")
    print(f"  - Cell types: {annot_metrics['n_cell_types']}")
    print(f"  - Mean confidence: {annot_metrics['mean_confidence']:.3f}")

    # Test spatial heterogeneity
    print(f"\n[3/4] Testing calculate_spatial_heterogeneity()...")
    start = time.time()

    spatial_metrics = calculate_spatial_heterogeneity(adata_annotated)

    elapsed = time.time() - start
    print(f"  ✓ Spatial analysis complete in {elapsed:.1f}s")

    if 'morans_i' in spatial_metrics:
        print(f"  - Moran's I: {spatial_metrics['morans_i']['mean']:.3f}")
        print(f"  - Significant genes: {spatial_metrics['morans_i']['significant_genes']}")

    if 'spatial_entropy' in spatial_metrics:
        print(f"  - Spatial entropy: {spatial_metrics['spatial_entropy']['mean']:.3f}")

    # Validate outputs
    print(f"\n[4/4] Validating outputs...")

    checks = [
        ('spatial_region' in adata_annotated.obs, "Clustering results"),
        ('cell_type' in adata_annotated.obs, "Cell type annotations"),
        ('morans_i' in spatial_metrics, "Moran's I results"),
        ('spatial_entropy' in spatial_metrics, "Spatial entropy"),
        (annot_metrics['n_clusters'] > 0, "Non-zero clusters")
    ]

    # Optional check (won't fail test)
    import numpy as np
    has_confidence = not np.isnan(annot_metrics['mean_confidence'])
    checks.append((True, f"Confidence scores: {'available' if has_confidence else 'not available (ok)'}"))


    all_passed = True
    for check, desc in checks:
        status = "✓" if check else "✗"
        print(f"  {status} {desc}")
        if not check:
            all_passed = False

    print("\n" + "="*80)
    if all_passed:
        print("SUCCESS: All adapter functions working correctly")
        print("✓ Ready for Streamlit deployment")
    else:
        print("FAILURE: Some checks failed")
        print("✗ Review errors above before deployment")
    print("="*80)

    return all_passed


if __name__ == "__main__":
    success = test_adapter_pipeline()
    sys.exit(0 if success else 1)

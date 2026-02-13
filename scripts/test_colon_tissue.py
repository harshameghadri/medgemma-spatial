"""
Multi-tissue validation: colon cancer Visium HD
Proves the tissue-agnostic annotation pipeline works on an unknown tissue type.
Run as: python scripts/test_colon_tissue.py
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc
import numpy as np
from src.streamlit_adapter import annotate_spatial_regions, calculate_spatial_heterogeneity

H5AD = Path("outputs/test_final_success_20260130_151551/e2735493904baddea063a55d5e676d24/extracted/square_008um.h5ad")

if not H5AD.exists():
    print(f"ERROR: H5AD not found at {H5AD}")
    print("Expected colon cancer Visium HD at this path.")
    sys.exit(1)

print("=" * 70)
print("MULTI-TISSUE VALIDATION: Colon Cancer Visium HD (tissue='Unknown')")
print("=" * 70)

print(f"\n[1/4] Loading {H5AD.name}...")
t0 = time.time()
adata = sc.read_h5ad(H5AD)
print(f"  Loaded: {adata.n_obs:,} spots x {adata.n_vars:,} genes ({time.time()-t0:.1f}s)")

# Filter for well-covered spots first (HD bins are 8µm - many sparse or empty)
counts_per_spot = np.asarray(adata.X.sum(axis=1)).flatten()
MIN_COUNTS = 200
n_before = adata.n_obs
adata = adata[counts_per_spot >= MIN_COUNTS].copy()
print(f"  Filtered to spots >= {MIN_COUNTS} counts: {adata.n_obs:,} / {n_before:,} ({adata.n_obs/n_before*100:.1f}%)")

# Subsample to 10K for speed
N_SUBSAMPLE = 10000
if adata.n_obs > N_SUBSAMPLE:
    print(f"  Subsampling to {N_SUBSAMPLE:,} spots...")
    sc.pp.subsample(adata, n_obs=N_SUBSAMPLE, random_state=42)
    print(f"  Subsampled: {adata.n_obs:,} spots")

print(f"\n[2/4] Running annotation (tissue='Unknown')...")
t0 = time.time()
try:
    adata_ann, metrics = annotate_spatial_regions(
        adata,
        resolution=0.5,
        use_markers=True,
        tissue="Unknown"
    )
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")
    print(f"  Clusters: {metrics['n_clusters']}")
    print(f"  Cell types: {metrics['n_cell_types']}")
    print(f"  Mean confidence: {metrics.get('mean_confidence', 'N/A')}")
    print(f"\n  Cell type distribution (top 10):")
    cell_type_counts = metrics.get('cell_type_counts', {})
    for ct, n in sorted(cell_type_counts.items(), key=lambda x: -x[1])[:10]:
        bar = "#" * int(n / max(cell_type_counts.values()) * 20)
        print(f"    {ct:<25} {n:>5}  {bar}")
    annotation_ok = True
except Exception as e:
    print(f"  ERROR: {e}")
    import traceback
    traceback.print_exc()
    annotation_ok = False
    adata_ann = adata
    metrics = {}

print(f"\n[3/4] Running spatial heterogeneity analysis...")
t0 = time.time()
try:
    spatial_metrics = calculate_spatial_heterogeneity(adata_ann)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")
    if 'morans_i' in spatial_metrics:
        mi = spatial_metrics['morans_i']
        print(f"  Moran's I: mean={mi['mean']:.3f}, sig_genes={mi['significant_genes']}")
    if 'spatial_entropy' in spatial_metrics:
        se = spatial_metrics['spatial_entropy']
        print(f"  Spatial entropy: {se['mean']:.3f} ± {se['std']:.3f}")
    if 'neighborhood_enrichment' in spatial_metrics:
        ne = spatial_metrics['neighborhood_enrichment']
        print(f"  Enriched cell type pairs: {ne['n_enriched_pairs']}")
    spatial_ok = True
except Exception as e:
    print(f"  ERROR: {e}")
    spatial_ok = False
    spatial_metrics = {}

print(f"\n[4/4] Validating results...")

# Check that colon-relevant cell types were found
colon_expected = ['Epithelial', 'Stromal_CAF', 'T_cells', 'Macrophages', 'Plasma_cells',
                  'B_cells', 'Endothelial', 'NK_cells', 'Luminal']
cell_type_keys = list(metrics.get('cell_type_counts', {}).keys())
found = [ct for ct in colon_expected if any(ct.lower() in k.lower() for k in cell_type_keys)]
print(f"  Colon-relevant cell types found: {found}")

checks = [
    (annotation_ok, "Annotation completed without error"),
    (metrics.get('n_clusters', 0) > 2, f"Clusters found: {metrics.get('n_clusters', 0)} > 2"),
    (metrics.get('n_cell_types', 0) > 3, f"Cell types found: {metrics.get('n_cell_types', 0)} > 3"),
    (len(found) >= 2, f"Colon-relevant types: {len(found)} >= 2"),
    (spatial_ok, "Spatial heterogeneity computed"),
    (spatial_metrics.get('morans_i', {}).get('mean', 0) >= 0, "Moran's I >= 0 (spatial structure)"),
]

print()
all_pass = True
for ok, desc in checks:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {desc}")
    if not ok:
        all_pass = False

print()
print("=" * 70)
if all_pass:
    print("RESULT: PASS - Tissue-agnostic pipeline works on colon cancer data")
    print(f"Cell types identified without tissue label: {cell_type_keys[:5]}")
else:
    print("RESULT: FAIL - See failures above")
print("=" * 70)

sys.exit(0 if all_pass else 1)

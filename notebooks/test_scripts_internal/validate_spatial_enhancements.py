"""
Quick validation script for spatial enhancements
Run after executing 01_scanpy_baseline.ipynb
"""

import anndata as ad
from pathlib import Path
import json

PROJECT_ROOT = Path.cwd().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"

print("="*60)
print("VALIDATION: Spatial Transcriptomics Enhancements")
print("="*60)

print("\n1. Checking processed h5ad file...")
h5ad_file = OUTPUT_DIR / "processed_visium.h5ad"
if h5ad_file.exists():
    adata = ad.read_h5ad(h5ad_file)

    checks = {
        "adata.obsm['spatial'] exists": 'spatial' in adata.obsm,
        "adata.obsp['spatial_connectivities'] exists": 'spatial_connectivities' in adata.obsp,
        "adata.uns['moranI'] exists": 'moranI' in adata.uns,
        "adata.uns['leiden_co_occurrence'] exists": 'leiden_co_occurrence' in adata.uns,
        "adata.obs['in_tissue'] exists": 'in_tissue' in adata.obs,
    }

    for check, result in checks.items():
        status = "✅" if result else "❌"
        print(f"  {status} {check}")

    if 'spatial' in adata.obsm:
        print(f"\n  Spatial coordinates shape: {adata.obsm['spatial'].shape}")

    if 'moranI' in adata.uns:
        morans_i = adata.uns['moranI']
        print(f"  Moran's I values range: [{morans_i['I'].min():.3f}, {morans_i['I'].max():.3f}]")
        print(f"  Mean Moran's I: {morans_i['I'].mean():.3f}")
else:
    print("  ❌ processed_visium.h5ad not found - run notebook first")

print("\n2. Checking JSON output...")
json_file = OUTPUT_DIR / "scanpy_features_spatial.json"
if json_file.exists():
    with open(json_file, 'r') as f:
        features = json.load(f)

    checks = {
        "metadata.analysis_type == 'spatial_transcriptomics'":
            features.get('metadata', {}).get('analysis_type') == 'spatial_transcriptomics',
        "spatial_statistics section exists": 'spatial_statistics' in features,
        "morans_i section exists":
            'morans_i' in features.get('spatial_statistics', {}),
        "spatial_cooccurrence section exists":
            'spatial_cooccurrence' in features.get('spatial_statistics', {}),
    }

    for check, result in checks.items():
        status = "✅" if result else "❌"
        print(f"  {status} {check}")

    if 'spatial_statistics' in features:
        spatial_stats = features['spatial_statistics']
        if 'morans_i' in spatial_stats:
            print(f"\n  Mean Moran's I (JSON): {spatial_stats['morans_i']['mean']:.3f}")
            print(f"  Significant genes: {spatial_stats['morans_i']['significant_genes_count']}")
        if 'spatial_cooccurrence' in spatial_stats:
            print(f"  Co-occurrence entries: {len(spatial_stats['spatial_cooccurrence'])}")
else:
    print("  ❌ scanpy_features_spatial.json not found - run notebook first")

print("\n3. Checking output images...")
expected_images = [
    "spatial_tissue_overview.png",
    "spatial_cooccurrence.png",
    "qc_violin.png",
    "highly_variable_genes.png",
    "pca_variance.png",
    "clustering_overview.png"
]

for img in expected_images:
    img_path = OUTPUT_DIR / img
    status = "✅" if img_path.exists() else "❌"
    print(f"  {status} {img}")

print("\n" + "="*60)
print("VALIDATION COMPLETE")
print("="*60)

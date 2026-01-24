#!/usr/bin/env python
"""
Execute Scanpy baseline spatial transcriptomics analysis
Generated from 01_scanpy_baseline.ipynb
"""

import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import squidpy as sq
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
from datetime import datetime
import psutil
import os

def main():
    print(f"scanpy version: {sc.__version__}")
    print(f"squidpy version: {sq.__version__}")
    print(f"anndata version: {ad.__version__}")

    # Configuration
    SEED = 42
    np.random.seed(SEED)

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "sample"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

print(f"Project root: {PROJECT_ROOT}")
print(f"Data directory: {DATA_DIR}")
print(f"Output directory: {OUTPUT_DIR}")

# Load Data
print("\n" + "="*60)
print("LOADING DATA")
print("="*60)
h5_file = list(DATA_DIR.glob("*.h5"))[0]
print(f"Loading: {h5_file.name}")

adata = sc.read_10x_h5(h5_file)
adata.var_names_make_unique()

print(f"\nDataset shape: {adata.shape}")
print(f"Spots (observations): {adata.n_obs}")
print(f"Genes (variables): {adata.n_vars}")

# Load Spatial Coordinates
print("\n" + "="*60)
print("LOADING SPATIAL COORDINATES")
print("="*60)
spatial_dir = DATA_DIR / "spatial"
tissue_positions = spatial_dir / "tissue_positions_list.csv"

spatial_coords = pd.read_csv(
    tissue_positions,
    header=None,
    names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
)

spatial_coords = spatial_coords.set_index('barcode')
spatial_coords = spatial_coords.loc[adata.obs_names]

adata.obsm['spatial'] = spatial_coords[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
adata.obs['in_tissue'] = spatial_coords['in_tissue'].values

import PIL
tissue_img = PIL.Image.open(spatial_dir / "tissue_hires_image.png")
adata.uns['spatial'] = {
    'Visium_Human_Breast_Cancer': {
        'images': {'hires': np.array(tissue_img)},
        'scalefactors': {
            'tissue_hires_scalef': 1.0,
            'spot_diameter_fullres': 89.43
        }
    }
}

print(f"Spatial coordinates loaded: {adata.obsm['spatial'].shape}")
print(f"Spots in tissue: {adata.obs['in_tissue'].sum()} / {len(adata)}")

# QC
print("\n" + "="*60)
print("QUALITY CONTROL")
print("="*60)
adata.var['mt'] = adata.var_names.str.startswith('MT-')

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True
)

print(f"Mean genes per spot: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"Mean counts per spot: {adata.obs['total_counts'].mean():.0f}")
print(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, ax=axes, show=False)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "qc_violin.png", dpi=150, bbox_inches='tight')
plt.close()

# Filtering
print("\n" + "="*60)
print("FILTERING")
print("="*60)
print(f"Before filtering: {adata.n_obs} spots, {adata.n_vars} genes")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} spots, {adata.n_vars} genes")

# Normalization
print("\n" + "="*60)
print("NORMALIZATION & PREPROCESSING")
print("="*60)
adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

n_hvg = adata.var['highly_variable'].sum()
print(f"Highly variable genes: {n_hvg}")

sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(OUTPUT_DIR / "highly_variable_genes.png", dpi=150, bbox_inches='tight')
plt.close()

# PCA
print("\n" + "="*60)
print("DIMENSIONALITY REDUCTION")
print("="*60)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack', random_state=SEED)

sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
plt.savefig(OUTPUT_DIR / "pca_variance.png", dpi=150, bbox_inches='tight')
plt.close()

# Clustering
print("\n" + "="*60)
print("CLUSTERING")
print("="*60)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=SEED)

sc.tl.umap(adata, random_state=SEED)

sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

cluster_counts = adata.obs['leiden'].value_counts().sort_index()
print(f"\nClusters found: {len(cluster_counts)}")
print("\nCluster sizes:")
print(cluster_counts)

# Build Spatial Neighbor Graph
print("\n" + "="*60)
print("SPATIAL NEIGHBOR GRAPH & MORAN'S I")
print("="*60)
print("Building spatial neighbor graph...")

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

print(f"Spatial connectivities: {adata.obsp['spatial_connectivities'].shape}")
print(f"Spatial distances: {adata.obsp['spatial_distances'].shape}")

print("\nComputing Moran's I for spatial autocorrelation...")

top_hvgs = adata.var_names[adata.var['highly_variable']][:100]

sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=top_hvgs,
    n_perms=100,
    n_jobs=1
)

morans_i = adata.uns['moranI'].sort_values('I', ascending=False)

print(f"\nTop 5 spatially autocorrelated genes:")
print(morans_i.head())
print(f"\nMean Moran's I: {morans_i['I'].mean():.4f}")
print(f"Significant genes (p < 0.05): {(morans_i['pval_norm'] < 0.05).sum()}")

# Spatial Visualizations
print("\n" + "="*60)
print("SPATIAL VISUALIZATIONS")
print("="*60)
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

sq.pl.spatial_scatter(
    adata,
    color='leiden',
    ax=axes[0, 0],
    title='Spatial Clusters on Tissue',
    size=1.5,
    legend_loc='right margin'
)

sq.pl.spatial_scatter(
    adata,
    color='total_counts',
    ax=axes[0, 1],
    title='Total Counts per Spot',
    size=1.5,
    cmap='viridis'
)

top_gene = morans_i.index[0]
sq.pl.spatial_scatter(
    adata,
    color=top_gene,
    ax=axes[1, 0],
    title=f'{top_gene} Expression (Moran I={morans_i.loc[top_gene, "I"]:.3f})',
    size=1.5,
    cmap='Reds'
)

sq.pl.spatial_scatter(
    adata,
    color='n_genes_by_counts',
    ax=axes[1, 1],
    title='Genes Detected per Spot',
    size=1.5,
    cmap='Blues'
)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "spatial_tissue_overview.png", dpi=150, bbox_inches='tight')
plt.close()

print(f"✅ Spatial tissue visualizations saved")

# Co-occurrence Analysis
print("\n" + "="*60)
print("CO-OCCURRENCE ANALYSIS")
print("="*60)
print("Analyzing spatial co-occurrence between clusters...")

sq.gr.co_occurrence(
    adata,
    cluster_key='leiden',
    spatial_key='spatial'
)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sq.pl.co_occurrence(
    adata,
    cluster_key='leiden',
    clusters=None,
    ax=axes[0]
)
axes[0].set_title('Spatial Co-occurrence Matrix')

co_occur = adata.uns['leiden_co_occurrence']
occurrence = co_occur['occurrence']
expected = np.outer(occurrence.sum(axis=1), occurrence.sum(axis=0)) / occurrence.sum()
enrichment = np.log2((occurrence + 1) / (expected + 1))

im = axes[1].imshow(enrichment, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
axes[1].set_xticks(range(len(occurrence)))
axes[1].set_yticks(range(len(occurrence)))
axes[1].set_xticklabels(occurrence.columns, rotation=90)
axes[1].set_yticklabels(occurrence.index)
axes[1].set_title('Spatial Enrichment (log2 obs/exp)')
axes[1].set_xlabel('Cluster')
axes[1].set_ylabel('Cluster')
plt.colorbar(im, ax=axes[1], label='log2(enrichment)')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "spatial_cooccurrence.png", dpi=150, bbox_inches='tight')
plt.close()

print(f"✅ Co-occurrence analysis complete")

# Export Features
print("\n" + "="*60)
print("EXPORTING FEATURES")
print("="*60)
cluster_info = {}
for cluster in adata.obs['leiden'].unique():
    cluster_mask = adata.obs['leiden'] == cluster
    cluster_info[str(cluster)] = {
        "count": int(cluster_mask.sum()),
        "mean_genes": float(adata.obs.loc[cluster_mask, 'n_genes_by_counts'].mean()),
        "mean_counts": float(adata.obs.loc[cluster_mask, 'total_counts'].mean())
    }

morans_top = morans_i.head(10)
morans_dict = {
    gene: {
        "morans_i": float(morans_top.loc[gene, 'I']),
        "p_value": float(morans_top.loc[gene, 'pval_norm'])
    }
    for gene in morans_top.index
}

co_occur_matrix = adata.uns['leiden_co_occurrence']['occurrence']
max_cooccur = []
for cluster in co_occur_matrix.index:
    others = co_occur_matrix.loc[cluster].drop(cluster)
    max_neighbor = others.idxmax()
    max_cooccur.append({
        "cluster": str(cluster),
        "most_colocalized_with": str(max_neighbor),
        "cooccurrence_count": int(others.max())
    })

features = {
    "metadata": {
        "analysis_date": datetime.now().isoformat(),
        "scanpy_version": sc.__version__,
        "squidpy_version": sq.__version__,
        "dataset": h5_file.name,
        "random_seed": SEED,
        "analysis_type": "spatial_transcriptomics"
    },
    "dataset_summary": {
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_highly_variable_genes": int(adata.var['highly_variable'].sum()),
        "spots_in_tissue": int(adata.obs['in_tissue'].sum())
    },
    "qc_metrics": {
        "mean_genes_per_spot": float(adata.obs['n_genes_by_counts'].mean()),
        "mean_counts_per_spot": float(adata.obs['total_counts'].mean()),
        "mean_pct_mt": float(adata.obs['pct_counts_mt'].mean())
    },
    "clustering": {
        "n_clusters": int(len(cluster_counts)),
        "resolution": 0.5,
        "clusters": cluster_info
    },
    "spatial_statistics": {
        "spatial_neighbors": {
            "n_neighbors": 6,
            "coord_type": "generic"
        },
        "morans_i": {
            "mean": float(morans_i['I'].mean()),
            "std": float(morans_i['I'].std()),
            "significant_genes_count": int((morans_i['pval_norm'] < 0.05).sum()),
            "top_genes": morans_dict
        },
        "spatial_cooccurrence": max_cooccur
    }
}

output_file = OUTPUT_DIR / "scanpy_features_spatial.json"
with open(output_file, 'w') as f:
    json.dump(features, f, indent=2)

print(f"\n✅ Enhanced spatial features exported to: {output_file}")
print(f"\nSpatial feature summary:")
print(f"  - Spatial neighbors computed: ✅")
print(f"  - Moran's I genes analyzed: {len(morans_i)}")
print(f"  - Significant spatial genes: {(morans_i['pval_norm'] < 0.05).sum()}")
print(f"  - Co-occurrence matrix: {co_occur_matrix.shape}")

# Save h5ad
output_h5ad = OUTPUT_DIR / "processed_visium.h5ad"
adata.write(output_h5ad)
print(f"✅ Processed AnnData saved to: {output_h5ad}")

# Final Summary
process = psutil.Process(os.getpid())
memory_mb = process.memory_info().rss / 1024 / 1024

print("\n" + "="*60)
print("SPATIAL TRANSCRIPTOMICS ANALYSIS SUMMARY")
print("="*60)
print(f"\nDataset: {h5_file.name}")
print(f"Spots analyzed: {adata.n_obs}")
print(f"Spots in tissue: {adata.obs['in_tissue'].sum()}")
print(f"Genes analyzed: {adata.n_vars}")
print(f"Clusters identified: {len(cluster_counts)}")
print(f"\nSpatial Analysis:")
print(f"  - Spatial neighbor graph: ✅")
print(f"  - Mean Moran's I: {morans_i['I'].mean():.4f}")
print(f"  - Significant spatial genes: {(morans_i['pval_norm'] < 0.05).sum()}")
print(f"  - Co-occurrence clusters: {len(co_occur_matrix)}")
print(f"\nMemory usage: {memory_mb:.0f} MB")
print(f"\nOutputs saved:")
print(f"  - Spatial features JSON: {output_file}")
print(f"  - Processed h5ad: {output_h5ad}")
print(f"  - Spatial tissue overview: {OUTPUT_DIR}/spatial_tissue_overview.png")
print(f"  - Co-occurrence analysis: {OUTPUT_DIR}/spatial_cooccurrence.png")
print(f"  - Other figures: {OUTPUT_DIR}/*.png")
print("\n✅ Spatial transcriptomics analysis complete!")
print("="*60)

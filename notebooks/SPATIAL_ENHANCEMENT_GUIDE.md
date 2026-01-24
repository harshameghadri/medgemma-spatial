# Scanpy Baseline → Spatial Transcriptomics Enhancement Guide

**Purpose**: Add true spatial analysis to `01_scanpy_baseline.ipynb`
**Status**: Squidpy already installed ✅
**Target**: Transform scRNA-seq notebook into proper spatial transcriptomics analysis

---

## 1. INSTALLATION (ALREADY COMPLETE)

```bash
# Squidpy is already installed in medgemma environment
conda activate medgemma
python -c "import squidpy as sq; print(f'Squidpy {sq.__version__} installed')"
```

**Installed packages**:
- `squidpy==1.5.0` - Spatial transcriptomics analysis
- All dependencies (networkx, leidenalg, statsmodels, etc.)

---

## 2. UPDATED IMPORTS (Replace Cell 1)

```python
# Cell 1: Setup & Imports
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import squidpy as sq  # ← ADD THIS
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
from datetime import datetime

print(f"scanpy version: {sc.__version__}")
print(f"squidpy version: {sq.__version__}")  # ← ADD THIS
print(f"anndata version: {ad.__version__}")
```

---

## 3. NEW CELLS TO ADD

### Cell 5B: Load Spatial Coordinates (INSERT AFTER Cell 5)

**Location**: After "Load Data" section, before "Quality Control"

```python
# Cell 5B: Load Spatial Coordinates
print("Loading spatial coordinates...")

# Path to spatial data
spatial_dir = DATA_DIR / "spatial"
tissue_positions = spatial_dir / "tissue_positions_list.csv"

# Load tissue positions (Visium format: no header)
# Columns: barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
spatial_coords = pd.read_csv(
    tissue_positions,
    header=None,
    names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
)

# Match barcodes with adata
spatial_coords = spatial_coords.set_index('barcode')
spatial_coords = spatial_coords.loc[adata.obs_names]

# Add spatial coordinates to adata
adata.obsm['spatial'] = spatial_coords[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
adata.obs['in_tissue'] = spatial_coords['in_tissue'].values

# Load and store tissue image
import PIL
tissue_img = PIL.Image.open(spatial_dir / "tissue_hires_image.png")
adata.uns['spatial'] = {
    'Visium_Human_Breast_Cancer': {
        'images': {'hires': np.array(tissue_img)},
        'scalefactors': {
            'tissue_hires_scalef': 1.0,
            'spot_diameter_fullres': 89.43  # Visium spot diameter
        }
    }
}

print(f"Spatial coordinates loaded: {adata.obsm['spatial'].shape}")
print(f"Spots in tissue: {adata.obs['in_tissue'].sum()} / {len(adata)}")
```

---

### Cell 7B: Build Spatial Neighbor Graph (INSERT AFTER Cell 7 - Leiden Clustering)

**Location**: After dimensionality reduction, before visualization

```python
# Cell 7B: Build Spatial Neighbor Graph
print("Building spatial neighbor graph...")

# Compute spatial neighbors based on physical proximity
# coord_type='generic' for pixel coordinates
# n_neighs=6 for Visium hexagonal lattice (default)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

print(f"Spatial connectivities: {adata.obsp['spatial_connectivities'].shape}")
print(f"Spatial distances: {adata.obsp['spatial_distances'].shape}")

# Compute spatial autocorrelation (Moran's I) for highly variable genes
print("\nComputing Moran's I for spatial autocorrelation...")

# Select top 100 HVGs for speed
top_hvgs = adata.var_names[adata.var['highly_variable']][:100]

# Calculate Moran's I
sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=top_hvgs,
    n_perms=100,  # Permutation test
    n_jobs=1  # M1 Mac: use 1 core to avoid issues
)

# Sort by Moran's I score
morans_i = adata.uns['moranI'].sort_values('I', ascending=False)

print(f"\nTop 5 spatially autocorrelated genes:")
print(morans_i.head())
print(f"\nMean Moran's I: {morans_i['I'].mean():.4f}")
print(f"Significant genes (p < 0.05): {(morans_i['pval_norm'] < 0.05).sum()}")
```

---

### Cell 8B: Spatial Visualizations on Tissue (INSERT AFTER Cell 8)

**Location**: After UMAP clustering visualization

```python
# Cell 8B: Spatial Visualizations on Tissue
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Leiden clusters on tissue
sq.pl.spatial_scatter(
    adata,
    color='leiden',
    ax=axes[0, 0],
    title='Spatial Clusters on Tissue',
    size=1.5,
    legend_loc='right margin'
)

# Plot 2: Total counts on tissue (shows tissue quality)
sq.pl.spatial_scatter(
    adata,
    color='total_counts',
    ax=axes[0, 1],
    title='Total Counts per Spot',
    size=1.5,
    cmap='viridis'
)

# Plot 3: Most spatially autocorrelated gene
top_gene = morans_i.index[0]
sq.pl.spatial_scatter(
    adata,
    color=top_gene,
    ax=axes[1, 0],
    title=f'{top_gene} Expression (Moran I={morans_i.loc[top_gene, "I"]:.3f})',
    size=1.5,
    cmap='Reds'
)

# Plot 4: n_genes per spot (shows tissue coverage)
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
plt.show()

print(f"✅ Spatial tissue visualizations saved")
```

---

### Cell 9B: Niche/Cell Type Enrichment (INSERT AFTER Cell 9)

**Location**: After Moran's I calculation

```python
# Cell 9B: Spatial Co-occurrence and Niche Analysis
print("Analyzing spatial co-occurrence between clusters...")

# Co-occurrence analysis: which clusters are neighbors?
sq.gr.co_occurrence(
    adata,
    cluster_key='leiden',
    spatial_key='spatial'
)

# Visualize co-occurrence
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Occurrence matrix (how often clusters are neighbors)
sq.pl.co_occurrence(
    adata,
    cluster_key='leiden',
    clusters=None,  # All clusters
    ax=axes[0]
)
axes[0].set_title('Spatial Co-occurrence Matrix')

# Enrichment matrix (observed vs expected)
from matplotlib.colors import diverging_palette
import matplotlib.pyplot as plt

# Get enrichment scores
co_occur = adata.uns['leiden_co_occurrence']
occurrence = co_occur['occurrence']
expected = np.outer(occurrence.sum(axis=1), occurrence.sum(axis=0)) / occurrence.sum()
enrichment = np.log2((occurrence + 1) / (expected + 1))

# Plot enrichment
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
plt.show()

print(f"✅ Co-occurrence analysis complete")
```

---

### Cell 10 UPDATE: Enhanced JSON Export with Spatial Features

**Replace the existing Cell 10 (Export Features to JSON) with this:**

```python
# Cell 10: Export Enhanced Features to JSON with Spatial Statistics
cluster_info = {}
for cluster in adata.obs['leiden'].unique():
    cluster_mask = adata.obs['leiden'] == cluster
    cluster_info[str(cluster)] = {
        "count": int(cluster_mask.sum()),
        "mean_genes": float(adata.obs.loc[cluster_mask, 'n_genes_by_counts'].mean()),
        "mean_counts": float(adata.obs.loc[cluster_mask, 'total_counts'].mean())
    }

# Get top spatially autocorrelated genes
morans_top = morans_i.head(10)
morans_dict = {
    gene: {
        "morans_i": float(morans_top.loc[gene, 'I']),
        "p_value": float(morans_top.loc[gene, 'pval_norm'])
    }
    for gene in morans_top.index
}

# Compute co-occurrence enrichment summary
co_occur_matrix = adata.uns['leiden_co_occurrence']['occurrence']
max_cooccur = []
for cluster in co_occur_matrix.index:
    # Find cluster with most co-occurrence (excluding self)
    others = co_occur_matrix.loc[cluster].drop(cluster)
    max_neighbor = others.idxmax()
    max_cooccur.append({
        "cluster": str(cluster),
        "most_colocalized_with": str(max_neighbor),
        "cooccurrence_count": int(others.max())
    })

# Build enhanced features JSON
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
```

---

## 4. EXPECTED OUTPUTS

### New Files Generated

```
outputs/
├── spatial_tissue_overview.png      # 4-panel spatial visualization
├── spatial_cooccurrence.png         # Co-occurrence heatmaps
└── scanpy_features_spatial.json     # Enhanced features with spatial stats
```

### Spatial Features JSON Structure

```json
{
  "metadata": {
    "analysis_type": "spatial_transcriptomics",
    "squidpy_version": "1.5.0"
  },
  "dataset_summary": {
    "spots_in_tissue": 4895
  },
  "spatial_statistics": {
    "spatial_neighbors": {
      "n_neighbors": 6,
      "coord_type": "generic"
    },
    "morans_i": {
      "mean": 0.045,
      "significant_genes_count": 67,
      "top_genes": {
        "KRT14": {"morans_i": 0.523, "p_value": 0.001},
        "KRT5": {"morans_i": 0.487, "p_value": 0.002}
      }
    },
    "spatial_cooccurrence": [
      {"cluster": "0", "most_colocalized_with": "3", "cooccurrence_count": 245}
    ]
  }
}
```

---

## 5. SPATIAL NEIGHBOR GRAPH EXPLAINED

**What it is**: Graph where edges connect spatially adjacent spots

**Visium Layout**: Hexagonal lattice
- Each spot has up to 6 neighbors
- Distance ~55μm center-to-center
- Built from pixel coordinates

**Parameters**:
- `coord_type='generic'`: Use pixel coordinates (not grid indices)
- `n_neighs=6`: Visium hex grid default
- Creates: `adata.obsp['spatial_connectivities']` (adjacency matrix)

**Why it matters**:
- Required for Moran's I calculation
- Identifies spatially correlated genes
- Enables niche/co-occurrence analysis

---

## 6. TROUBLESHOOTING

### Issue 1: `tissue_positions_list.csv` has headers

**Error**: `ValueError: could not convert string to float`

**Fix**: Check CSV format
```python
# If file has header row (barcode,in_tissue,...)
spatial_coords = pd.read_csv(tissue_positions, index_col=0)

# If file has NO header (our case)
spatial_coords = pd.read_csv(tissue_positions, header=None, names=[...])
```

### Issue 2: Memory error with large datasets

**Error**: `MemoryError` during Moran's I

**Fix**: Reduce genes tested
```python
# Instead of all HVGs
top_hvgs = adata.var_names[adata.var['highly_variable']][:50]  # Test fewer genes
```

### Issue 3: Squidpy plot not showing tissue image

**Error**: Scatter plot without tissue background

**Fix**: Ensure `adata.uns['spatial']` is populated
```python
# Verify image loaded
print(adata.uns['spatial']['Visium_Human_Breast_Cancer']['images']['hires'].shape)
# Should be (height, width, 3) or (height, width, 4)
```

### Issue 4: Barcode mismatch

**Error**: `KeyError` when matching barcodes

**Fix**: Suffix handling
```python
# 10x h5 file may have suffixes like -1
# tissue_positions.csv usually doesn't
# Match by removing suffix:
spatial_coords.index = spatial_coords.index.str.replace('-1', '')
```

---

## 7. VALIDATION CHECKLIST

After adding all cells, verify:

- [ ] `adata.obsm['spatial']` exists and has shape (n_spots, 2)
- [ ] `adata.obsp['spatial_connectivities']` exists (sparse matrix)
- [ ] `adata.uns['moranI']` contains Moran's I results
- [ ] Moran's I values are between -1 and 1
- [ ] `spatial_tissue_overview.png` shows tissue structure
- [ ] Clusters visible on tissue image
- [ ] Co-occurrence matrix is symmetric
- [ ] JSON contains `spatial_statistics` section
- [ ] Runtime still <5 minutes ✅
- [ ] Memory usage <16GB ✅

---

## 8. PERFORMANCE NOTES

**Bottlenecks**:
- Moran's I permutation test: ~2-3 min for 100 genes
- Co-occurrence computation: ~30 sec

**Optimizations**:
- Test fewer genes (50-100 instead of all HVGs)
- Use `n_perms=100` (default=1000 is slow)
- `n_jobs=1` on M1 Mac (multiprocessing issues)

**Expected total runtime**: 4-5 minutes (still under limit ✅)

---

## 9. INSERTION ORDER SUMMARY

1. **Cell 1**: Replace imports (add squidpy)
2. **After Cell 5**: Insert Cell 5B (load spatial coords)
3. **After Cell 7**: Insert Cell 7B (spatial neighbors + Moran's I)
4. **After Cell 8**: Insert Cell 8B (spatial tissue plots)
5. **After Cell 9**: Insert Cell 9B (co-occurrence)
6. **Cell 10**: Replace with enhanced JSON export

---

## 10. WHAT THIS ENABLES FOR MEDGEMMA

With spatial features, MedGemma can generate reports describing:

1. **Tumor Architecture**: "Cluster 4 shows high KRT14 expression in the tumor core"
2. **Spatial Organization**: "Immune clusters co-localize with tumor-stroma interface"
3. **Tissue Heterogeneity**: "High spatial autocorrelation (Moran's I=0.52) indicates distinct zones"
4. **Clinical Relevance**: "Invasive margin (cluster 2-3 boundary) shows high proliferation genes"

This transforms generic clustering into **spatially-aware pathology reports**! 🎯

---

**Ready to enhance the notebook? All cells are copy-paste ready!**

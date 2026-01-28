# Enhanced Spatial Statistics Analysis

Comprehensive spatial metrics generation for MedGemma clinical reports.

## Files

### Jupyter Notebook
- **05_spatial_statistics.ipynb**: Interactive notebook with full analysis and visualizations
  - Requires: `squidpy` package
  - Best for: Exploration and development

### Python Scripts

#### 1. run_spatial_stats.py (Preferred)
Full-featured version using squidpy for advanced spatial statistics.

**Requirements:**
```bash
pip install squidpy
```

**Usage:**
```bash
# Default (uses annotated_visium_enhanced.h5ad)
python3 notebooks/run_spatial_stats.py

# Custom input
python3 notebooks/run_spatial_stats.py /path/to/input.h5ad /path/to/output_dir
```

**Features:**
- Ripley's K function (multi-scale clustering)
- Neighborhood enrichment (squidpy implementation)
- Spatial entropy
- Nearest neighbor distances
- Extended spatial autocorrelation
- Cluster compactness

#### 2. run_spatial_stats_scanpy.py (Fallback)
Scanpy-only version without squidpy dependency.

**Requirements:**
```bash
pip install scanpy scipy scikit-learn
```

**Usage:**
```bash
# Default
python3 notebooks/run_spatial_stats_scanpy.py

# Custom input
python3 notebooks/run_spatial_stats_scanpy.py /path/to/input.h5ad /path/to/output_dir
```

**Features:**
- Neighborhood enrichment (manual implementation)
- Spatial entropy
- Nearest neighbor distances
- Moran's I (manual calculation)
- Cluster compactness

**Note:** Use this version if you encounter import errors with squidpy (e.g., zarr compatibility issues).

## Outputs

All scripts generate the following files in `outputs/`:

### JSON Files
1. **spatial_statistics_enhanced.json**
   - Comprehensive spatial metrics
   - Cell type distributions
   - Cluster statistics
   - Ready for MedGemma prompt inclusion

### Text Files
2. **clinical_spatial_summary.txt**
   - Human-readable clinical interpretation
   - Spatial patterns summary
   - Neighborhood interactions
   - Top spatially variable genes

### Visualizations
3. **neighborhood_enrichment.png** - Cell type co-location heatmap
4. **spatial_entropy.png** - Local heterogeneity map
5. **nearest_neighbor_distances.png** - Dispersion patterns
6. **spatial_autocorrelation_extended.png** - Gene spatial patterns
7. **cluster_compactness.png** - Cluster organization metrics

### Data Files
8. **annotated_visium_spatial_stats.h5ad** - Updated AnnData with new columns:
   - `spatial_entropy`
   - `nn_distance`

## Spatial Metrics Explained

### 1. Spatial Entropy
**What it measures:** Local diversity of cell types
**Interpretation:**
- High (>1.0): Highly mixed neighborhoods
- Moderate (0.5-1.0): Some mixing
- Low (<0.5): Homogeneous regions

### 2. Nearest Neighbor Distance
**What it measures:** How clustered each cell type is
**Interpretation:**
- Tightly clustered (<50 pixels): Strong spatial organization
- Moderate (50-150 pixels): Mixed distribution
- Dispersed (>150 pixels): Scattered throughout tissue

### 3. Neighborhood Enrichment
**What it measures:** Which cell types prefer to be neighbors
**Interpretation:**
- Z-score > 2: Significant enrichment
- Z-score < -2: Significant depletion
- Identifies spatial niches and interfaces

### 4. Cluster Compactness
**What it measures:** Spatial organization of Leiden clusters
**Interpretation:**
- Compactness score > 1.5: Tightly organized
- Score 0.8-1.5: Moderate organization
- Score < 0.8: Dispersed/fragmented

### 5. Spatial Autocorrelation (Moran's I)
**What it measures:** Genes with spatially correlated expression
**Interpretation:**
- Moran's I > 0.3: Strong spatial pattern
- Moran's I near 0: Random distribution
- Moran's I < 0: Checkerboard pattern (rare)

## Integration with MedGemma

The generated `spatial_statistics_enhanced.json` can be directly included in MedGemma prompts:

```python
import json

# Load spatial metrics
with open('outputs/spatial_statistics_enhanced.json') as f:
    spatial_metrics = json.load(f)

# Create MedGemma prompt
prompt = f"""
Based on the following spatial transcriptomics analysis:

Sample: {spatial_metrics['sample_info']['n_spots']} spots
Cell types: {', '.join(spatial_metrics['sample_info']['cell_types'])}

Spatial heterogeneity: {spatial_metrics['spatial_entropy']['overall']['interpretation']}
Top spatially variable genes: {', '.join(spatial_metrics['spatial_autocorrelation']['top_genes'][:5])}

Generate a clinical pathology report...
"""
```

## Troubleshooting

### Import Error: squidpy
```
ImportError: cannot import name 'ArrayNotFoundError' from 'zarr.errors'
```
**Solution:** Use `run_spatial_stats_scanpy.py` instead

### Memory Error
**Solution:** Reduce number of genes tested in Moran's I calculation
- Edit script: Change `gene_list = adata.var_names[:100]` to use fewer genes

### Missing spatial coordinates
```
KeyError: 'spatial'
```
**Solution:** Ensure input h5ad has `adata.obsm['spatial']`
- Use output from `01_spatial_analysis.ipynb` or similar

## Performance

**Typical runtime on M1 Mac (64GB RAM):**
- 4895 spots, 2000 genes
- ~2-3 minutes total
- Memory usage: <8GB

## Next Steps

After running spatial statistics:
1. Review `clinical_spatial_summary.txt` for interpretation
2. Integrate `spatial_statistics_enhanced.json` into MedGemma prompts
3. Use generated visualizations in reports
4. Continue to MedGemma report generation (Week 2)

## Citation

If using for publication, cite:
- **Scanpy**: Wolf et al., Genome Biology 2018
- **Squidpy**: Palla et al., Nature Methods 2022
- **Spatial entropy method**: Custom implementation based on Shannon entropy

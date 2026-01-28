# Spatial Metrics Quick Reference Card

**For:** MedGemma Clinical Report Integration
**Updated:** 2026-01-28

---

## 1. Spatial Entropy

**What:** Shannon entropy of cell type diversity in local neighborhoods

**Value Range:** 0.0 - 2.5
- **Low (<0.5):** Homogeneous regions, single cell type dominates
- **Moderate (0.5-1.0):** Mixed neighborhoods, 2-3 cell types
- **High (>1.0):** Highly heterogeneous, diverse cell types

**Clinical Interpretation:**
- High entropy in tumor → heterogeneous tumor microenvironment
- Low entropy → well-defined tissue compartments
- Entropy gradients → interface regions (e.g., tumor-stroma boundary)

**Example JSON:**
```json
{
  "spatial_entropy": {
    "overall": {
      "median": 0.85,
      "interpretation": "moderate"
    },
    "by_cell_type": {
      "Tumor": {"mean_entropy": 0.72},
      "Immune": {"mean_entropy": 1.23}
    }
  }
}
```

**Prompt Usage:**
> "The tissue shows MODERATE spatial heterogeneity (entropy: 0.85), with immune cells in highly mixed regions (1.23) and tumor cells in more homogeneous areas (0.72)."

---

## 2. Nearest Neighbor Distance

**What:** Distance (in pixels) to nearest spot of same cell type

**Value Range:** 0 - 500+ pixels
- **Tightly clustered (<50):** Strong spatial organization
- **Moderate (50-150):** Mixed distribution
- **Dispersed (>150):** Scattered throughout tissue

**Clinical Interpretation:**
- Tightly clustered tumor → compact tumor mass
- Dispersed immune cells → infiltrating T cells
- Moderate stromal cells → fibrous tissue bands

**Example JSON:**
```json
{
  "nearest_neighbor_distances": {
    "Tumor": {
      "median_distance": 45.3,
      "spatial_pattern": "tightly_clustered"
    },
    "Immune": {
      "median_distance": 178.9,
      "spatial_pattern": "dispersed"
    }
  }
}
```

**Prompt Usage:**
> "Tumor cells form a TIGHTLY CLUSTERED mass (median NN distance: 45px), while immune cells are DISPERSED throughout (178px), suggesting infiltration."

---

## 3. Neighborhood Enrichment

**What:** Z-scores indicating which cell types prefer to be neighbors

**Value Range:** -5 to +5 (Z-score)
- **Enriched (>2):** Significantly more neighbors than expected
- **Random (-2 to 2):** Expected co-location
- **Depleted (<-2):** Significantly fewer neighbors than expected

**Clinical Interpretation:**
- Tumor-Immune enrichment → immune infiltration
- Tumor-Stroma enrichment → desmoplastic reaction
- Immune-Immune enrichment → tertiary lymphoid structures
- Tumor-Tumor depletion → invasive phenotype

**Example JSON:**
```json
{
  "neighborhood_enrichment": {
    "Tumor": {
      "enriched_neighbors": ["Immune", "Stroma"],
      "depleted_neighbors": [],
      "n_enriched": 2
    },
    "Immune": {
      "enriched_neighbors": ["Tumor"],
      "depleted_neighbors": ["Epithelial"],
      "n_enriched": 1
    }
  }
}
```

**Prompt Usage:**
> "Significant spatial niches identified: Tumor cells enriched near Immune and Stroma (p<0.05), suggesting active tumor-immune interface with desmoplastic response."

---

## 4. Cluster Compactness

**What:** Spatial organization metrics for Leiden clusters

**Components:**
- **Compactness Score:** mean_distance / std_distance
- **Density:** spots / convex_hull_area
- **Convex Hull Area:** Total spatial extent

**Compactness Score Range:** 0.5 - 3.0
- **Tightly organized (>1.5):** Compact, well-defined cluster
- **Moderate (0.8-1.5):** Some spatial coherence
- **Dispersed (<0.8):** Fragmented, scattered spots

**Clinical Interpretation:**
- High compactness → distinct tissue compartment
- Low compactness → invasive or fragmented region
- High density + high compactness → dense cellular infiltrate

**Example JSON:**
```json
{
  "cluster_compactness": {
    "0": {
      "n_spots": 1247,
      "compactness_score": 1.82,
      "density": 45.3,
      "convex_hull_area": 27543.2
    },
    "1": {
      "n_spots": 892,
      "compactness_score": 0.65,
      "density": 18.7,
      "convex_hull_area": 47721.9
    }
  }
}
```

**Prompt Usage:**
> "Cluster 0 (n=1247) shows TIGHTLY ORGANIZED structure (compactness: 1.82), while Cluster 1 (n=892) is DISPERSED (0.65), suggesting invasive phenotype."

---

## 5. Spatial Autocorrelation (Moran's I)

**What:** Measure of gene expression spatial clustering

**Value Range:** -1.0 to +1.0
- **Strong positive (>0.3):** Gene expression spatially clustered
- **Random (-0.1 to 0.1):** No spatial pattern
- **Negative (<-0.1):** Checkerboard pattern (rare in biology)

**Clinical Interpretation:**
- High Moran's I → zonation, compartmentalization
- Spatially variable genes → functional domains
- CD8A high I → organized immune infiltrate
- EPCAM high I → epithelial compartment boundaries

**Example JSON:**
```json
{
  "spatial_autocorrelation": {
    "n_genes_tested": 100,
    "n_significant": 47,
    "mean_morans_i": 0.34,
    "top_genes": ["CD8A", "EPCAM", "KRT19", "CD3D", "VIM"],
    "top_morans_i": [0.82, 0.78, 0.76, 0.74, 0.71]
  }
}
```

**Prompt Usage:**
> "47% of genes show significant spatial patterns. Top spatially variable genes (CD8A: 0.82, EPCAM: 0.78) indicate organized immune infiltrate and epithelial compartments."

---

## Combining Metrics for Clinical Interpretation

### Pattern 1: Hot Tumor (Inflamed)
```
✓ Tumor-Immune neighborhood enrichment (Z>2)
✓ High spatial entropy in tumor regions (>1.0)
✓ Dispersed immune cells (NN distance >150px)
✓ CD8A high Moran's I (>0.5)
→ "Active immune infiltration with T cell presence"
```

### Pattern 2: Cold Tumor (Desert)
```
✓ Tumor-Immune depletion (Z<-2)
✓ Low entropy in tumor (entropy <0.5)
✓ Tightly clustered tumor (NN distance <50px)
✓ Immune compactness low (<0.8)
→ "Immune-excluded phenotype with compact tumor mass"
```

### Pattern 3: Desmoplastic Response
```
✓ Tumor-Stroma enrichment (Z>2)
✓ Stroma tightly clustered (NN <50px)
✓ High VIM/COL1A1 Moran's I (>0.6)
✓ Stroma compactness high (>1.5)
→ "Prominent desmoplastic reaction with dense fibrosis"
```

### Pattern 4: Invasive Front
```
✓ Low tumor compactness (<0.8)
✓ High tumor entropy (>1.0)
✓ Multiple cell type enrichments
✓ High n_enriched for tumor (≥3)
→ "Invasive tumor front with complex microenvironment"
```

---

## MedGemma Prompt Template

```python
prompt = f"""
Generate a clinical pathology report for spatial transcriptomics analysis:

SPATIAL ORGANIZATION:
- Heterogeneity: {metrics['spatial_entropy']['overall']['interpretation'].upper()}
  (entropy: {metrics['spatial_entropy']['overall']['median']:.2f})

CELL TYPE PATTERNS:
{format_cell_patterns(metrics['nearest_neighbor_distances'])}

SPATIAL NICHES:
{format_enrichment(metrics['neighborhood_enrichment'])}

GENE EXPRESSION PATTERNS:
- {metrics['spatial_autocorrelation']['n_significant']} spatially variable genes
- Key markers: {', '.join(metrics['spatial_autocorrelation']['top_genes'][:5])}

CLUSTER ORGANIZATION:
{format_compactness(metrics['cluster_compactness'])}

Synthesize these findings into a 200-word clinical report focusing on:
1. Tumor microenvironment architecture
2. Immune infiltration patterns
3. Stromal organization
4. Clinical/prognostic implications
"""
```

---

## Quick Troubleshooting

| Issue | Check | Fix |
|-------|-------|-----|
| All entropy = 0 | Single cell type | Expected if homogeneous |
| NN distance = NaN | <2 spots of type | Filter rare types |
| No enrichment | Small sample | Expected for <500 spots |
| Moran's I all ~0 | No spatial structure | Check if pre-clustered |
| Compactness = NaN | Cluster too small | Filter clusters <10 spots |

---

## File References

**Input:** `outputs/annotated_visium_enhanced.h5ad`
**Output:** `outputs/spatial_statistics_enhanced.json`
**Script:** `notebooks/run_spatial_stats_scanpy.py`
**Docs:** `notebooks/README_spatial_stats.md`

---

**Next Action:** Integrate these metrics into Week 2 MedGemma prompt engineering

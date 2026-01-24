# Spatial Transcriptomics Analysis: Algorithm Explanations

**Purpose**: Explain key algorithm choices in spatial transcriptomics workflow
**Audience**: Beginners in spatial biology and computational methods

---

## 1. Highly Variable Genes (HVGs)

### What are they?

Genes that show high variability in expression across spots, beyond technical noise.

### Why identify them?

**Problem**: 20,000+ genes, but most are:
- Housekeeping genes (constant expression)
- Low-expressed noise (technical variation)
- Biologically uninformative

**Solution**: Focus on **2,000 highly variable genes** that capture biological differences

### Algorithm: Seurat Method (`flavor='seurat'`)

**Steps**:
1. **Bin genes** by mean expression level (20 bins)
2. **Calculate variance** within each bin
3. **Standardize variance** to Z-scores (variance/expected)
4. **Select top genes** with highest standardized variance

**Parameters in our code**:
```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
```

- `n_top_genes=2000`: Standard choice balancing information vs computation
- `flavor='seurat'`: Robust to mean-variance relationship

### Why 2,000 genes?

| # Genes | Pros | Cons |
|---------|------|------|
| 500 | Fast, less noise | May miss biology |
| **2,000** | **Standard, well-tested** | **Balanced** ✅ |
| 5,000 | More complete | Slower, more noise |

**Benchmark studies** (Luecken & Theis 2019) show 2,000 captures 90%+ of variation.

### What it looks like

**High variance gene**: KRT14 (keratin)
- Basal cells: 1000 counts
- Luminal cells: 10 counts
- **Variance**: 📈 HIGH → Selected as HVG

**Low variance gene**: ACTB (actin, housekeeping)
- All cells: ~500 counts
- **Variance**: 📉 LOW → Not selected

---

## 2. Moran's I Spatial Autocorrelation

### What is it?

**Statistical measure** of spatial autocorrelation: "Do nearby spots have similar gene expression?"

**Range**: -1 to +1
- **+1**: Perfect positive correlation (neighbors are similar)
- **0**: Random spatial pattern
- **-1**: Perfect negative correlation (neighbors are opposite)

### Why use it?

**Goal**: Identify genes with **spatial structure** (not just high/low expression)

**Example**:
- **Gene A**: High in some spots, low in others, but **randomly distributed** → Moran's I ≈ 0
- **Gene B**: High in tumor core, low at edges, **spatially organized** → Moran's I = 0.6 ✅

Gene B is more biologically interesting for spatial analysis!

### Algorithm

**Formula**:
```
I = (N / W) * Σᵢ Σⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄) / Σᵢ (xᵢ - x̄)²
```

Where:
- `N` = number of spots
- `wᵢⱼ` = spatial weight (1 if neighbors, 0 otherwise)
- `xᵢ` = expression in spot i
- `x̄` = mean expression

**In plain English**:
1. For each spot, look at its neighbors (typically 6 in Visium)
2. If spot is high-expressing and neighbors are also high → positive contribution
3. If spot is high but neighbors are low → negative contribution
4. Average across all spots

### Permutation Test (Statistical Significance)

**Problem**: Is Moran's I = 0.25 significant or just random?

**Solution**: Permutation test
1. Randomly shuffle gene expression across spots (break spatial structure)
2. Recalculate Moran's I
3. Repeat 100-1000 times
4. Compare observed I to null distribution
5. If observed I > 95% of permutations → p < 0.05 (significant)

**Parameters in our code**:
```python
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)
```

- `n_perms=100`: Balance between speed and accuracy (default=1000 is slower)
- `n_jobs=1`: M1 Mac compatibility (avoid multiprocessing issues)

### Interpretation

| Moran's I | p-value | Interpretation |
|-----------|---------|----------------|
| 0.6 | 0.001 | **Strong spatial pattern** (e.g., tumor marker) |
| 0.3 | 0.02 | **Moderate pattern** (e.g., ECM gene) |
| 0.1 | 0.15 | **Weak/random** (housekeeping gene) |
| -0.2 | 0.05 | **Negative correlation** (rare, e.g., mutually exclusive cell types) |

### Real Example: Breast Cancer

**High Moran's I genes**:
- **KRT14** (I=0.52): Basal marker, forms contiguous regions
- **ERBB2** (I=0.45): HER2, clusters in tumor subregions
- **DCN** (I=0.38): Decorin, marks stromal areas

**Low Moran's I genes**:
- **GAPDH** (I=0.05): Housekeeping, ubiquitous
- **Ribosomal proteins** (I≈0): Random distribution

---

## 3. Spatial Neighbor Graph

### What is it?

**Graph structure** where:
- **Nodes** = Spots
- **Edges** = Connect spatially adjacent spots

### Why build it?

**Foundation for all spatial analyses**:
- Moran's I calculation
- Spatial smoothing
- Co-occurrence analysis
- Niche identification

### Visium-Specific Geometry

**Hexagonal Lattice**:
```
    ●   ●   ●
  ●   ●   ●   ●
    ●   ●   ●
```

Each spot has **up to 6 neighbors** (hexagonal packing)

**Spot specifications**:
- Diameter: 55 μm
- Center-to-center: 100 μm
- Each spot: ~10 cells

### Algorithm: k-Nearest Neighbors (kNN)

**Steps**:
1. Extract spatial coordinates (x, y pixels) from `tissue_positions_list.csv`
2. For each spot, find k nearest neighbors by Euclidean distance
3. Create adjacency matrix: 1 if neighbors, 0 otherwise

**Parameters in our code**:
```python
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
```

- `coord_type='generic'`: Use pixel coordinates (not grid indices)
- `n_neighs=6`: Match Visium hexagonal geometry

### Distance Types

| coord_type | Use Case | Our Choice |
|------------|----------|------------|
| `'grid'` | Regular square lattice | ❌ (Visium is hexagonal) |
| `'generic'` | **Pixel coordinates** | **✅** (Correct for Visium) |

### What it creates

**Outputs stored in AnnData**:

1. **`adata.obsp['spatial_connectivities']`**
   - Sparse matrix (4895 × 4895)
   - Binary: 1 if neighbors, 0 otherwise
   - Symmetric (if A→B then B→A)

2. **`adata.obsp['spatial_distances']`**
   - Euclidean distances (in pixels)
   - Only for connected neighbors
   - Used for distance-weighted analyses

**Example**:
```
Spot AAACAAGTATCTCCCA-1 neighbors:
  - AAACACCAATAACTGC-1 (dist: 123.4 pixels)
  - AAACAGAGCGACTCCT-1 (dist: 118.9 pixels)
  - AAACCGTGCTTCCGAT-1 (dist: 125.1 pixels)
  ... (up to 6 total)
```

---

## 4. Spatial Co-occurrence Analysis

### What is it?

**Quantify which clusters are spatially close** to each other.

**Question**: Do tumor cells neighbor immune cells more than expected by chance?

### Algorithm

**Steps**:
1. For each spot in cluster A, count neighbors in cluster B
2. Build **occurrence matrix**: O[A,B] = # of A spots with B neighbors
3. Calculate **expected**: E[A,B] = (size of A) × (size of B) / total
4. Compute **enrichment**: log₂(O[A,B] / E[A,B])

**Parameters in our code**:
```python
sq.gr.co_occurrence(adata, cluster_key='leiden', spatial_key='spatial')
```

### Interpretation

**Occurrence Matrix**:
```
      Cluster 0  Cluster 1  Cluster 2
  0      245         89         12
  1       87        156         34
  2       15         32         78
```
Row: focal cluster, Column: neighbor cluster

**Enrichment** (log₂ obs/exp):
- **+2**: 4× more co-occurrence than expected (strong attraction)
- **0**: Random (observed = expected)
- **-2**: 4× less co-occurrence (avoidance/segregation)

### Biological Examples

**Tumor-Immune Interface**:
```
Cluster 3 (Tumor) ↔ Cluster 7 (Immune): enrichment = +1.8
→ Immune cells infiltrate tumor regions
```

**Mutually Exclusive**:
```
Cluster 2 (Basal) ↔ Cluster 4 (Luminal): enrichment = -2.1
→ Cell types segregate spatially
```

---

## 5. Why Not Just Use Correlation?

### Problem with Pearson Correlation

**Our old code** (from original notebook):
```python
correlation = np.corrcoef(expression[:-1], expression[1:])[0, 1]
```

**Issues**:
1. **Assumes linear order**: Treats spots as 1D sequence
2. **Ignores spatial structure**: Spot #100 may be far from #101
3. **Not a real autocorrelation**: Just correlation of shifted values
4. **Wrong null distribution**: Can't test significance

**Result**: Meaningless values, not interpretable as spatial statistics

### Why Moran's I is Correct

✅ **Uses true spatial neighbors** (from graph)
✅ **Accounts for 2D geometry** (x, y coordinates)
✅ **Permutation test** for significance
✅ **Standard metric** in spatial statistics (published literature)

---

## 6. Leiden Clustering

### What is it?

**Graph-based clustering** that groups spots with similar gene expression.

### Why Leiden (not k-means)?

| Method | Pros | Cons |
|--------|------|------|
| k-means | Fast, simple | Requires # clusters, assumes spherical |
| **Leiden** | **No # clusters needed**, **detects modules** | **More complex** |
| Louvain | Popular, fast | Can miss small communities |

**Leiden** = Improved Louvain algorithm (guarantees well-connected communities)

### Algorithm

**Steps**:
1. Build **k-NN graph** from PCA space (not spatial!)
   - Each spot connected to k similar spots by gene expression
2. **Optimize modularity**: Group spots to maximize within-group connections
3. **Iterate** until modularity stops improving

**Parameters in our code**:
```python
sc.tl.leiden(adata, resolution=0.5, random_state=SEED)
```

- `resolution=0.5`: Controls granularity (higher = more clusters)
- `random_state=42`: Reproducibility

### Resolution Parameter

| Resolution | # Clusters | Use Case |
|------------|------------|----------|
| 0.1 | 3-5 | Broad cell types |
| **0.5** | **8-12** | **Balanced** ✅ |
| 1.0 | 15-25 | Fine subtypes |
| 2.0 | 30+ | Over-clustering |

**Rule of thumb**: Start with 0.5, adjust based on biology

---

## 7. PCA (Principal Component Analysis)

### Why reduce dimensions?

**Problem**: 2,000 HVGs → 2,000-dimensional space
**Issue**: Curse of dimensionality, noise, computation

**Solution**: PCA finds **40 principal components** capturing most variation

### How it works

**Intuition**: Find new axes that explain maximum variance

**Steps**:
1. Center data (subtract mean)
2. Compute covariance matrix (2000×2000)
3. Find eigenvectors (directions of maximum variance)
4. Project data onto top 40 eigenvectors

**Result**: 2,000 genes → 40 PCs, retaining ~80% variance

### Why 40 components?

**Elbow method**: Plot variance explained vs # PCs

```
Variance
  │  ●
  │   ●
  │    ●
  │     ●___
  │         ●___●___●___
  └────────────────────── PCs
         40
```

After ~40 PCs, adding more gives diminishing returns.

**Default**: 40-50 PCs is standard for scRNA-seq/Visium

---

## 8. Normalization: Why log1p?

### Raw Counts Problem

**Spots have different sequencing depth**:
- Spot A: 5,000 total counts
- Spot B: 20,000 total counts

**Gene expression confounded by depth!**

### Solution: Normalize + Log Transform

**Step 1**: Library size normalization
```python
sc.pp.normalize_total(adata, target_sum=1e4)
```
Scale each spot to 10,000 total counts

**Step 2**: Log transform
```python
sc.pp.log1p(adata)
```

**Why log(x+1)?**
- **Variance stabilization**: High counts dominate less
- **Approximate normality**: Better for downstream stats
- **+1**: Avoid log(0) = -∞

**Effect**:
```
Raw:  [0, 1, 10, 100, 1000]
Log:  [0, 0.7, 2.4, 4.6, 6.9]
```
Compresses range, makes differences more comparable

---

## 9. QC Filtering Thresholds

### Why filter?

**Low-quality spots**:
- Tissue edge/damage
- Low cell density
- Failed library prep

**Can distort analysis** (artificial clusters of bad spots)

### Our Thresholds

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| `min_genes` | 200 | Typical Visium: 1000-5000 genes/spot |
| `min_cells` | 3 | Gene in ≥3 spots (avoid noise) |
| `pct_mt` | (not used) | <20% for scRNA-seq, less critical for Visium |

**Visium-specific**: Higher min_genes than scRNA-seq (spots = ~10 cells)

---

## 10. Random Seed: Why 42?

### What is it?

**Seed for random number generator** → reproducible results

### Where randomness occurs

- PCA initialization (ARPACK algorithm)
- UMAP embedding
- Leiden clustering (stochastic optimization)
- Moran's I permutation test

### Why 42?

**Hitchhiker's Guide to the Galaxy reference** 😄

**Real reason**: Any fixed seed works, 42 is community convention

**Important**: Same seed + same data + same code = **identical results**

---

## Summary Table

| Algorithm | Purpose | Key Parameter | Why This Choice? |
|-----------|---------|---------------|-------------------|
| **HVG Selection** | Find informative genes | `n_top_genes=2000` | Standard, balances info/speed |
| **PCA** | Dimension reduction | `n_pcs=40` | Captures ~80% variance |
| **k-NN Graph** | Build spot similarity | `n_neighbors=10` | Robust community detection |
| **Leiden** | Clustering | `resolution=0.5` | Balanced granularity |
| **Spatial Neighbors** | Spatial graph | `n_neighs=6` | Visium hexagonal lattice |
| **Moran's I** | Spatial autocorr | `n_perms=100` | Speed vs accuracy tradeoff |
| **Normalization** | Scale counts | `target_sum=1e4` | Standard CPM (counts per 10k) |

---

## Further Reading

**Scanpy tutorials**: https://scanpy-tutorials.readthedocs.io/
**Squidpy tutorials**: https://squidpy.readthedocs.io/
**Best practices**: https://www.sc-best-practices.org/
**Moran's I paper**: Moran (1950) Biometrika
**Leiden paper**: Traag et al (2019) Scientific Reports

---

**Last Updated**: 2026-01-24

# Notebooks Directory

Development notebooks for spatial transcriptomics analysis.

---

## 📁 Notebooks

### Week 1: Exploration

**01_scanpy_baseline.ipynb** ✅
- **Status**: Ready
- **Purpose**: Baseline spatial analysis using Scanpy
- **Runtime**: <5 minutes
- **Memory**: <16GB
- **Outputs**:
  - `outputs/scanpy_features.json`
  - `outputs/processed_visium.h5ad`
  - QC and clustering plots

**02_loki_test.ipynb** ⏳
- **Status**: Not started
- **Purpose**: Test Loki spatial foundation model
- **Environment**: `loki_env` (separate conda environment)
- **Timeline**: Week 1 Day 3-4

**03_medgemma_integration.ipynb** ⏳
- **Status**: Not started
- **Purpose**: Integrate MedGemma for report generation
- **Input**: Features JSON from 01 or 02
- **Output**: Clinical pathology report
- **Timeline**: Week 1 Day 5-7

**kaggle_submission.ipynb** ⏳
- **Status**: Not started
- **Purpose**: Final Kaggle competition submission
- **Timeline**: Week 4

---

## 🚀 Quick Start

### 1. Download Sample Data

```bash
# From project root
bash scripts/download_data.sh

# Or manually:
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5 -P data/sample/
```

### 2. Activate Environment

```bash
conda activate medgemma
```

### 3. Launch Jupyter

```bash
# From project root
jupyter notebook

# Or from notebooks directory
cd notebooks
jupyter notebook
```

### 4. Run Notebook

1. Open `01_scanpy_baseline.ipynb`
2. Select kernel: `medgemma`
3. Run all cells: Cell → Run All

---

## 📊 Expected Outputs

After running `01_scanpy_baseline.ipynb`:

```
outputs/
├── scanpy_features.json       # Feature extraction for MedGemma
├── processed_visium.h5ad       # Processed AnnData object
├── qc_violin.png               # QC metrics visualization
├── highly_variable_genes.png   # HVG selection
├── pca_variance.png            # PCA variance ratio
└── clustering_overview.png     # UMAP + cluster sizes
```

---

## 🔧 Troubleshooting

### Issue: Kernel not found

```bash
# Install kernel
conda activate medgemma
python -m ipykernel install --user --name medgemma --display-name "MedGemma"

# Restart Jupyter
```

### Issue: Data file not found

```bash
# Check data directory
ls ../data/sample/

# Download if missing
bash ../scripts/download_data.sh
```

### Issue: Out of memory

```python
# In notebook, use backed mode:
adata = sc.read_10x_h5(h5_file, backed='r')
```

### Issue: ImportError

```bash
# Verify environment
conda activate medgemma
python -c "import scanpy as sc; print(sc.__version__)"

# Reinstall if needed
pip install -r ../requirements.txt
```

---

## 📝 Development Guidelines

### Notebook Style

```python
# ✅ Good: One-line docstrings
def process_data(adata):
    """QC and filter Visium spots."""
    return filtered_adata

# ❌ Bad: Verbose comments
def process_data(adata):
    # This function processes the data
    # First we do QC
    # Then we filter
    ...
```

### Random Seeds

```python
# Always set seeds for reproducibility
SEED = 42
np.random.seed(SEED)
sc.tl.leiden(adata, random_state=SEED)
```

### Paths

```python
# ✅ Use relative paths from project root
PROJECT_ROOT = Path.cwd().parent
DATA_DIR = PROJECT_ROOT / "data" / "sample"

# ❌ Never hardcode absolute paths
# data_dir = "/Users/you/Desktop/data"  # NO!
```

### Before Committing

1. **Restart & Run All**: Ensure notebook runs end-to-end
2. **Clear outputs** (optional): Cell → All Output → Clear
3. **Check file sizes**: Don't commit >10MB outputs

---

## 🎯 Week 1 Milestones

- [x] Day 1-2: Scanpy baseline notebook created
- [ ] Day 3-4: Loki test (optional, time-boxed)
- [ ] Day 5-7: MedGemma integration

---

## 📚 Resources

- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [10x Visium Datasets](https://www.10xgenomics.com/datasets)
- [Spatial Transcriptomics Best Practices](https://www.sc-best-practices.org/trajectories/spatial.html)

---

**Last Updated**: 2026-01-24

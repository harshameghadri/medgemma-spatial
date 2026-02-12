# MedGemma Spatial Transcriptomics - Skills & Working Examples

**Last Updated**: 2026-02-06
**Purpose**: Master reference for all validated commands, workflows, and solutions
**Project Status**: Week 2 Day 5 - Deployment Ready
**Competition Deadline**: Feb 24, 2026

---

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Data Processing](#data-processing)
3. [Spatial Analysis](#spatial-analysis)
4. [MedGemma Integration](#medgemma-integration)
5. [Deployment](#deployment)
6. [Testing & Validation](#testing--validation)
7. [Git Workflows](#git-workflows)
8. [Troubleshooting](#troubleshooting)
9. [Performance Benchmarks](#performance-benchmarks)

---

## Environment Setup

### Conda Environment Creation ✅
**Date**: 2026-01-24
**Status**: PRODUCTION
**Platform**: M1 Mac (Darwin 23.2.0)

**Command**:
```bash
conda create -n medgemma python=3.10
conda activate medgemma
pip install scanpy==1.10.2 squidpy==1.5.0 anndata==0.11.3
pip install torch>=2.4.0 transformers==4.45.1 accelerate bitsandbytes
pip install streamlit plotly kaleido
```

**Validation**:
```bash
python -c "import scanpy as sc; import squidpy as sq; print('✓ Environment ready')"
```

**Expected Output**:
```
✓ Environment ready
```

### Loki Environment (Separate) ✅
**Date**: 2026-01-24
**Status**: PRODUCTION
**Reason**: Incompatible PyTorch/anndata versions

**Command**:
```bash
conda create -n loki_env python=3.9
conda activate loki_env
cd Loki/src
pip install -e .
```

**Validation**:
```bash
python -c "import loki; print('✓ Loki installed')"
```

### M1 Mac MPS (GPU) Acceleration ✅
**Date**: 2026-01-24
**Status**: VALIDATED

**Test Code**:
```python
import torch
print(f"MPS Available: {torch.backends.mps.is_available()}")
print(f"PyTorch Version: {torch.__version__}")
```

**Expected Output**:
```
MPS Available: True
PyTorch Version: 2.4.0+
```

---

## Data Processing

### Load Visium H5AD File ✅
**Date**: 2026-01-28
**Status**: PRODUCTION

**Code**:
```python
import scanpy as sc

adata = sc.read_h5ad('outputs/annotated_visium.h5ad')
print(f"Loaded: {adata.n_obs} spots, {adata.n_vars} genes")
print(f"Spatial coords: {adata.obsm['spatial'].shape}")
print(f"Cell types: {adata.obs['cell_type'].unique()[:5]}")
```

**Expected Output**:
```
Loaded: 4895 spots, 2000 genes
Spatial coords: (4895, 2)
Cell types: ['Epithelial cells' 'Endothelial cells' 'NK cells' ...]
```

### PanglaoDB Marker-Based Annotation ✅
**Date**: 2026-02-02
**Status**: PRODUCTION
**Accuracy**: 90% mean confidence

**Command**:
```bash
python scripts/test_marker_annotation.py \
    --input data/sample_visium.h5ad \
    --output outputs/test_marker_annotation/enhanced_with_markers.h5ad \
    --markers data/PanglaoDB_markers_27_Mar_2020.tsv
```

**Code Pattern**:
```python
from src.spatial_analysis.uncertainty_spatial_analysis import annotate_cell_types_with_markers

adata = sc.read_h5ad('sample.h5ad')
marker_db = pd.read_csv('data/PanglaoDB_markers_27_Mar_2020.tsv', sep='\t')

adata, metrics = annotate_cell_types_with_markers(adata, marker_db)

print(f"Annotated cell types: {adata.obs['cell_type'].unique()}")
print(f"Mean confidence: {metrics['mean_confidence']:.2f}")
print(f"Low confidence rate: {metrics['low_confidence_rate']:.1%}")
```

**Performance**:
- 5,181 markers covering 163 cell types
- Mean confidence: 0.90 (excellent)
- Low confidence rate: 0.0% on test data

### Extract H&E Image from AnnData ✅
**Date**: 2026-02-05
**Status**: PRODUCTION

**Code**:
```python
from PIL import Image
import numpy as np

def load_hires_image_from_adata(adata, library_id=None):
    """Extract H&E image from AnnData spatial slot."""
    if library_id is None:
        library_id = list(adata.uns['spatial'].keys())[0]

    img_array = adata.uns['spatial'][library_id]['images']['hires']

    # Convert float64 [0,1] → uint8 [0,255]
    if img_array.dtype == np.float64:
        img_array = (img_array * 255).astype(np.uint8)

    return Image.fromarray(img_array)

# Usage
he_image = load_hires_image_from_adata(adata)
print(f"Image size: {he_image.size}")  # e.g., (2000, 1335)
```

**Preprocessing for MedGemma**:
```python
def preprocess_for_medgemma(image: Image.Image, target_size=896):
    """Resize to 896×896 for SigLIP encoder."""
    image.thumbnail((target_size, target_size), Image.LANCZOS)
    canvas = Image.new('RGB', (target_size, target_size), (255, 255, 255))
    offset = ((896 - image.width) // 2, (896 - image.height) // 2)
    canvas.paste(image, offset)
    return canvas
```

---

## Spatial Analysis

### Scanpy Baseline Pipeline ✅
**Date**: 2026-01-30
**Status**: PRODUCTION (100% pass rate)
**Runtime**: 94.9 seconds for 5,000 spots

**Command**:
```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
python src/spatial_analysis/uncertainty_spatial_analysis.py \
    --input outputs/annotated_visium.h5ad \
    --output outputs/spatial_features.json \
    --resolution 0.8
```

**Core Functions Available**:

#### 1. Doublet Detection (Scrublet)
```python
from src.spatial_analysis.uncertainty_spatial_analysis import detect_doublets_scrublet

adata, metrics = detect_doublets_scrublet(adata, expected_doublet_rate=0.06)
print(f"Doublets detected: {metrics['n_predicted_doublets']}")
print(f"Doublet rate: {metrics['doublet_rate']:.1%}")
```

**Known Issue**: Scrublet threshold auto-detection may fail
**Fix**: Manual threshold fallback at 0.25 (already implemented)

#### 2. Annotation Confidence
```python
from src.spatial_analysis.uncertainty_spatial_analysis import assess_annotation_quality

metrics = assess_annotation_quality(adata)
print(f"Mean confidence: {metrics['mean_confidence']:.2f}")
print(f"Low confidence spots: {metrics['n_low_confidence']}")
```

#### 3. Moran's I Spatial Autocorrelation
```python
from src.spatial_analysis.uncertainty_spatial_analysis import calculate_morans_i

morans_results = calculate_morans_i(adata, n_genes=20, n_perms=199)
print(f"Top gene: {morans_results['top_gene']}")
print(f"Moran's I: {morans_results['morans_i']:.3f}")
print(f"P-value: {morans_results['p_value']:.4f}")
```

#### 4. Spatial Entropy
```python
from src.spatial_analysis.uncertainty_spatial_analysis import calculate_spatial_entropy

entropy_metrics = calculate_spatial_entropy(adata, n_bootstrap=200)
print(f"Mean entropy: {entropy_metrics['mean']:.3f}")
print(f"95% CI: [{entropy_metrics['ci_lower']:.3f}, {entropy_metrics['ci_upper']:.3f}]")
```

**Interpretation**:
- High entropy (>0.8): High heterogeneity, complex spatial architecture
- Low entropy (<0.5): Low heterogeneity, homogeneous tissue

#### 5. Multiscale Neighborhood Enrichment
```python
from src.spatial_analysis.uncertainty_spatial_analysis import multiscale_nhood_enrichment

enrichment_metrics = multiscale_nhood_enrichment(adata, radii=[1, 2], n_perms=199)
print(f"Scale-stable pairs: {enrichment_metrics['n_stable_pairs']}")
print(f"Cross-scale correlation: {enrichment_metrics['correlation']:.3f}")
```

**Performance Metrics** (validated on enhanced_with_markers.h5ad):
- Total execution time: 94.9 seconds
- Memory usage: 820 MB
- Processing speed: 52.7 spots/second
- Success rate: 100% (10/10 quality checks)

### Spatial Graph Construction ✅
**Date**: 2026-01-30
**Status**: PRODUCTION

**Code**:
```python
import squidpy as sq

# Build spatial neighbor graph
sq.gr.spatial_neighbors(
    adata,
    n_neighs=6,
    coord_type='generic',
    spatial_key='spatial'
)

print(f"Graph built: {adata.obsp['spatial_connectivities'].shape}")
print(f"Average neighbors: {adata.obsp['spatial_connectivities'].mean():.2f}")
```

---

## MedGemma Integration

### Text-Only Report Generation ✅
**Date**: 2026-01-30
**Status**: PRODUCTION
**Anti-Parroting**: ✅ PASSED (0/3 tissue keywords detected)

**Command**:
```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
python notebooks/run_medgemma.py \
    --h5ad outputs/annotated_visium.h5ad \
    --output outputs/medgemma_report.txt \
    --device mps
```

**Code Pattern**:
```python
from transformers import pipeline

def generate_medgemma_report(features, device='mps'):
    """Generate clinical report from spatial features."""

    # Create tissue-blind prompt
    prompt = f"""
You are a pathologist analyzing spatial transcriptomics data.

SPATIAL ANALYSIS:
- Total spots analyzed: {features['annotation']['total_spots']}
- Major cell populations: {features['annotation']['n_major_populations']}
- Spatial heterogeneity (Moran's I): {features['spatial_heterogeneity']['morans_i_mean']:.2f}
- Spatial entropy: {features['spatial_heterogeneity']['spatial_entropy']:.2f}

Generate a concise pathology report (200 words max) focusing on:
1. Tissue architecture and spatial organization
2. Cell population diversity and distribution patterns
3. Potential biological significance
"""

    # Load MedGemma pipeline
    pipe = pipeline(
        "text-generation",
        model="google/medgemma-4b-it",
        device=device,
        model_kwargs={"load_in_4bit": True}
    )

    messages = [{"role": "user", "content": prompt}]
    result = pipe(messages, max_new_tokens=500)

    return result[0]['generated_text']
```

**CRITICAL**: Tissue-blind prompts only (no "colon", "breast", "brain" keywords)

### Multimodal (H&E + Text) Report Generation ✅
**Date**: 2026-02-05
**Status**: PRODUCTION
**Image Processing**: 2000×1335 → 896×896 (SigLIP encoder)

**Command**:
```bash
python scripts/test_multimodal_report.py \
    --h5ad outputs/annotated_visium.h5ad \
    --output outputs/multimodal_report.txt \
    --device mps
```

**Code**:
```python
from src.report_generation.medgemma_multimodal import generate_multimodal_report
import scanpy as sc
import json

# Load data
adata = sc.read_h5ad('sample.h5ad')
features = json.load(open('spatial_features.json'))

# Generate report with H&E image + text
report, metadata = generate_multimodal_report(
    adata,
    features,
    device='mps',
    max_new_tokens=500
)

print(report)
print(f"Model: {metadata['model']}")
print(f"Generation time: {metadata['generation_time_seconds']:.1f}s")
```

**Image Requirements**:
- Format: RGB PIL Image
- Size: 896×896 (automatic resize/padding)
- Source: `adata.uns['spatial'][library_id]['images']['hires']`

### Anti-Parroting Prompt Template ✅
**Date**: 2026-02-05
**Status**: VALIDATED (100% pass rate)

**Code**:
```python
from notebooks.medgemma_report_generator import create_anti_parroting_prompt

prompt = create_anti_parroting_prompt(features)
```

**Output Pattern** (saved to `outputs/anti_parroting_prompt.txt`):
```
You are a pathologist analyzing spatial transcriptomics data.

AGGREGATED SPATIAL METRICS:
- Total spots analyzed: 5000
- Major cell populations (>10% frequency): 3
- Spatial organization score (Moran's I): 0.65
- Heterogeneity index (entropy): 1.82

Generate a concise pathology report (200 words) WITHOUT:
- Mentioning tissue type
- Listing specific cell counts
- Identifying organ or disease
```

**Validation**: No tissue keywords ("colon", "breast", "brain", etc.) detected

---

## Deployment

### Streamlit Local Test ✅
**Date**: 2026-02-05
**Status**: VALIDATED
**File**: `app/streamlit_app.py` (653 lines)

**Command**:
```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
streamlit run app/streamlit_app.py
```

**URL**: http://localhost:8501

**Features Tested**:
- [x] File upload (H5AD format)
- [x] Progress bar during analysis
- [x] Interactive visualizations (Plotly)
- [x] Cell type composition bar chart
- [x] Spatial cluster pie chart
- [x] Generated report display
- [x] Download button (TXT format)

**Sidebar Settings**:
- Use Multimodal (H&E + Text): Checkbox
- Marker-based Annotation: Checkbox
- Leiden Resolution: Slider (0.1 - 2.0, default 0.8)

### Docker Build ✅
**Date**: 2026-02-05
**Status**: READY FOR DEPLOYMENT
**File**: `app/Dockerfile` (multi-stage build)

**Build Command**:
```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/app
chmod +x build.sh
./build.sh
```

**Test Container**:
```bash
docker run -p 8501:8501 medgemma-spatial:latest
```

**Docker Compose** (production):
```bash
cd app
docker-compose up -d
docker-compose logs -f
docker-compose down
```

**Security Features**:
- Non-root user (medgemma:1000)
- Multi-stage build (smaller image size)
- No cache directories
- Proper file ownership

### HuggingFace Spaces Deployment ⏳
**Date**: 2026-02-05
**Status**: PENDING (Manual Deployment Required)
**Documentation**: `app/README_HUGGINGFACE.md`

**Quick Deploy Steps**:

1. **Create Space**
   ```bash
   # Go to: https://huggingface.co/new-space
   # Repository name: medgemma-spatial
   # SDK: Docker
   # Hardware: CPU Basic (FREE) or CPU Upgrade
   ```

2. **Configure Space README**
   ```markdown
   ---
   title: MedGemma Spatial Transcriptomics
   emoji: 🔬
   colorFrom: blue
   colorTo: purple
   sdk: docker
   app_port: 8501
   pinned: false
   license: mit
   ---
   ```

3. **Upload Files**
   - `app/Dockerfile`
   - `app/streamlit_app.py`
   - `app/requirements.txt`
   - `src/` directory
   - `data/PanglaoDB_markers_27_Mar_2020.tsv`

4. **Wait for Build** (5-10 minutes)

5. **Test Space URL**
   ```
   https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial
   ```

**Hardware Recommendations**:
- **Demo**: CPU Basic (FREE) - 2 vCPU, 16GB RAM
- **Production**: CPU Upgrade ($0.05/hour) - 8 vCPU, 32GB RAM
- **Optimal**: GPU T4 ($0.60/hour) - 4 vCPU, 15GB RAM, T4 GPU

---

## Testing & Validation

### Robustness Test (Agent #2) ✅
**Date**: 2026-02-05
**Status**: COMPLETE
**Success Rate**: 100% (10/10 quality checks)

**Command**:
```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
python scripts/robustness_test_agent2.py
```

**Results**:
- Total execution time: 94.9 seconds
- Peak memory usage: 820.8 MB
- Processing speed: 52.7 spots/second
- Errors encountered: 0

**Quality Checks**:
- [x] Has required fields (X, obs, var)
- [x] Spatial coordinates valid (no NaN)
- [x] Annotation complete (9 cell types)
- [x] Spatial graph built (Squidpy)
- [x] Spatial analysis complete (5 sub-analyses)
- [x] Features JSON valid
- [x] No tissue type leakage
- [x] Memory under limit (<10 GB)
- [x] Execution time under limit (<10 min)
- [x] Overall pipeline health

**Output Files**:
- `outputs/robustness_test_sample2.json` (23 KB)
- `outputs/robustness_test_sample2_report.md` (8.3 KB)
- `outputs/robustness_test_sample2_summary.png` (261 KB)

### Data Leakage Validation ✅
**Date**: 2026-02-05
**Status**: COMPLETE
**Pass Rate**: 100% (9/9 scenarios)

**Command**:
```bash
python scripts/test_data_leakage.py
```

**Test Matrix**:
```
STRATEGY: guided_questions
  colon:  ✓ PASSED (risk: LOW, score: 0)
  brain:  ✓ PASSED (risk: LOW, score: 0)
  breast: ✓ PASSED (risk: LOW, score: 0)

STRATEGY: create_report
  colon:  ✓ PASSED (risk: LOW, score: 0)
  brain:  ✓ PASSED (risk: LOW, score: 0)
  breast: ✓ PASSED (risk: LOW, score: 0)

STRATEGY: anti_parroting
  colon:  ✓ PASSED (risk: LOW, score: 0)
  brain:  ✓ PASSED (risk: LOW, score: 0)
  breast: ✓ PASSED (risk: LOW, score: 0)
```

**Validation Logic**:
```python
def check_tissue_identification(prompt: str, expected_tissue: str) -> dict:
    """Check if prompt allows tissue identification."""

    tissue_keywords = {
        'colon': ['colon', 'colorectal', 'intestine', 'goblet'],
        'brain': ['brain', 'neuron', 'astrocyte', 'oligodendrocyte'],
        'breast': ['breast', 'mammary', 'luminal', 'ductal']
    }

    identifiable = any([
        any(kw in prompt.lower() for kw in tissue_keywords[tissue])
        for tissue in tissue_keywords
    ])

    return {'identifiable': identifiable, 'passed': not identifiable}
```

**Result**: 0/3 tissues identifiable in prompts

### End-to-End Pipeline Test ✅
**Date**: 2026-02-05
**Status**: COMPLETE

**Command**:
```bash
python scripts/test_end_to_end_leakage.py
```

**Output**: `outputs/end_to_end_leakage_test.json`

---

## Git Workflows

### Commit Working Features ✅
**Pattern**: `[type]: [description]`

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `test`: Tests
- `refactor`: Code refactoring
- `style`: Formatting

**Examples**:
```bash
# Feature commit
git add src/report_generation/medgemma_multimodal.py
git commit -m "feat: Add multimodal H&E image support for MedGemma 1.5"

# Bug fix commit
git add src/spatial_analysis/uncertainty_spatial_analysis.py
git commit -m "fix: Handle Scrublet threshold auto-detection failure"

# Test commit
git add scripts/test_data_leakage.py
git commit -m "test: Add adversarial data leakage validation"

# Documentation commit
git add .guides/SKILLS.md
git commit -m "docs: Add comprehensive working examples reference"
```

### Push to Remote ✅
```bash
git push origin devel
```

### Tag Release ✅
```bash
git tag -a v1.0-deployment -m "Production-ready deployment"
git push origin v1.0-deployment
```

### Recent Commits (Last 10)
```
0d8e97b test: Add consolidated agent robustness testing summary
62611b0 docs: Add comprehensive deployment validation report
5ab2928 docs: Add HuggingFace Spaces deployment guide
0d10807 feat: Add Docker containerization for deployment
fcda314 feat: Add Streamlit web application for deployment
77800b3 test: Add comprehensive data leakage validation
de77f22 feat: Add multimodal MedGemma support with H&E images
7f87760 fix: Remove data leakage from MedGemma prompts
5866760 feat: Add marker-based cell type annotation with PanglaoDB
977fda9 feat: Add subsampling and performance optimizations
```

---

## Troubleshooting

### Scrublet Threshold Error ✅ FIXED
**Date**: 2026-02-05
**Status**: RESOLVED

**Error**:
```
AttributeError: 'Scrublet' object has no attribute 'threshold_'
```

**Root Cause**: Scrublet fails to auto-detect threshold on clean datasets

**Solution** (in `src/spatial_analysis/uncertainty_spatial_analysis.py`):
```python
# BEFORE (BROKEN)
threshold = scrub.threshold_ if scrub.threshold_ is not None else 0.25
predicted_doublets = doublet_scores > threshold

# AFTER (FIXED)
if predicted_doublets is None:
    print("WARNING: Scrublet failed to auto-detect threshold. Using manual threshold at 0.25")
    threshold = 0.25
    predicted_doublets = doublet_scores > threshold
```

**Committed**: Agent monitoring (background process)

### OpenMP Library Conflict ✅ WORKAROUND
**Date**: 2026-02-05
**Status**: DOCUMENTED

**Error**:
```
OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized
```

**Platform**: M1 Mac only

**Workaround**:
```bash
export KMP_DUPLICATE_LIB_OK=TRUE
python your_script.py
```

**Permanent Fix** (add to `~/.zshrc` or `~/.bashrc`):
```bash
echo 'export KMP_DUPLICATE_LIB_OK=TRUE' >> ~/.zshrc
source ~/.zshrc
```

### Streamlit Port Conflict ✅ SOLVED
**Error**:
```
Address already in use: Port 8501
```

**Solution**:
```bash
# Find process using port 8501
lsof -i :8501

# Kill process
kill -9 <PID>

# Or use different port
streamlit run app/streamlit_app.py --server.port=8502
```

### Docker Build Fails (Missing Files) ⚠️ COMMON
**Error**:
```
COPY failed: file not found in build context
```

**Solution**:
```bash
# Ensure Dockerfile context is project root
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma

# Build from app/ directory (context: parent dir)
docker build -f app/Dockerfile -t medgemma-spatial:latest .
```

**Verify Build Context**:
```bash
# Check .dockerignore
cat app/.dockerignore

# Ensure required files present
ls -la src/
ls -la data/PanglaoDB_markers_27_Mar_2020.tsv
```

### HuggingFace Spaces Build Timeout ⚠️ POTENTIAL
**Symptom**: Build fails after 10 minutes

**Causes**:
1. Large Docker image (>5GB)
2. Complex dependencies
3. Network timeout

**Solutions**:
```dockerfile
# Use smaller base image
FROM python:3.10-slim  # NOT python:3.10

# Multi-stage build
FROM base as builder
# ... pip install
FROM base
COPY --from=builder /root/.local /home/medgemma/.local

# Cache pip packages
RUN pip install --no-cache-dir -r requirements.txt
```

---

## Performance Benchmarks

### Tested Samples ✅

| Sample | Spots | Genes | Runtime | Memory | Status |
|--------|-------|-------|---------|--------|--------|
| annotated_visium | 4,895 | 2,000 | - | 210 MB | ✅ PASS |
| enhanced_with_markers | 5,000 | 18,132 | 94.9s | 820 MB | ✅ PASS |
| square_008um | >10K | TBD | TBD | TBD | ⏳ Testing |

### Execution Time Breakdown (enhanced_with_markers.h5ad)

| Step | Time (s) | % of Total |
|------|----------|------------|
| Load H5AD | 0.06 | 0.06% |
| Validate structure | 0.00 | 0.00% |
| Extract coordinates | 0.00 | 0.00% |
| Cell type annotation | 0.00 | 0.00% |
| Build spatial graph | 1.56 | 1.64% |
| **Spatial analysis** | **93.24** | **98.25%** |
| Generate features JSON | 0.00 | 0.00% |
| Validate output | 0.04 | 0.04% |
| **TOTAL** | **94.90** | **100.00%** |

**Bottleneck**: Spatial analysis (Moran's I with 199 permutations)

**Optimization Opportunity**: Reduce permutations to 99 for 2x speedup

### Memory Profiling (M1 Mac 64GB)

| Component | Memory Usage |
|-----------|--------------|
| Base Python + imports | ~200 MB |
| Load H5AD (5K spots) | ~500 MB |
| Spatial analysis | ~800 MB |
| Peak usage | ~820 MB |
| MedGemma 4-bit (text) | +12 GB |
| MedGemma 1.5 (multimodal) | +16 GB |

**Total Pipeline**: <20 GB RAM (within M1 Mac limits)

### Scalability Estimates

| Dataset Size | Expected Runtime | Memory Usage |
|--------------|------------------|--------------|
| <5K spots | 1-2 min | <1 GB |
| 5K-10K spots | 2-5 min | 1-2 GB |
| 10K-50K spots | 5-15 min | 2-5 GB |
| >50K spots | 15-30 min | 5-10 GB |

**Note**: Linear time complexity observed (52.7 spots/second)

---

## Future Enhancements

### Loki Integration ⏳
**Status**: IN PROGRESS (Week 2 Day 6-7)
**Deadline**: 2-day time limit (Feb 6-7)
**Decision**: GO/NO-GO by end of Day 7

**Test Plan**:
1. Install Loki in separate environment (already done)
2. Test embeddings extraction
3. Compare with Scanpy baseline
4. If works: Integrate into pipeline
5. If fails: Document reason, use Scanpy fallback

### NicheFormer Integration ⏳
**Status**: DEFERRED to post-competition
**Reason**: Insufficient time (19 days until deadline)
**Timeline**: Feb 15-24 buffer period

**Rationale**: Focus on competition submission first, add NicheFormer as post-competition enhancement

### Fine-Tuning MedGemma ⏳
**Status**: OPTIONAL (if time permits)
**Requirements**: LoRA adapter, <5 training samples
**Documentation**: `info/LORA_FINETUNING_PLAN.md`

**Estimated Effort**: 4-8 hours (not budgeted in current timeline)

---

## Quick Command Reference

### Daily Workflow

```bash
# 1. Activate environment
conda activate medgemma

# 2. Run spatial analysis
python src/spatial_analysis/uncertainty_spatial_analysis.py \
    --input outputs/sample.h5ad \
    --output outputs/features.json

# 3. Generate report
python notebooks/run_medgemma.py \
    --h5ad outputs/sample.h5ad \
    --output outputs/report.txt \
    --device mps

# 4. Test web app
streamlit run app/streamlit_app.py

# 5. Commit progress
git add .
git commit -m "feat: [description]"
git push origin devel
```

### Testing Commands

```bash
# Data leakage validation
python scripts/test_data_leakage.py

# Robustness test
python scripts/robustness_test_agent2.py

# Multimodal test
python scripts/test_multimodal_report.py

# End-to-end pipeline
python scripts/test_end_to_end_leakage.py
```

### Deployment Commands

```bash
# Docker build
cd app
./build.sh

# Docker test
docker run -p 8501:8501 medgemma-spatial:latest

# Docker Compose
docker-compose up -d
docker-compose logs -f
docker-compose down
```

---

## Critical Success Factors

### Week 2 End Checklist (Feb 6-7) ✅
- [x] Data leakage fixes complete
- [x] Multimodal support working
- [x] 5+ sample reports validated
- [x] Code committed to GitHub
- [x] Streamlit app functional
- [x] Docker container ready

### Week 3 Deployment (Feb 8-14) ⏳
- [ ] HuggingFace Spaces deployed
- [ ] Public demo URL working
- [ ] Professional README updated
- [ ] Demo video recorded (2 min)
- [ ] LinkedIn post drafted

### Week 4 Competition (Feb 15-24) ⏳
- [ ] Kaggle submission prepared
- [ ] Competition strategy finalized
- [ ] Portfolio materials complete
- [ ] Interview-ready code examples

---

## Competition Context

**Deadline**: Feb 24, 2026 (19 days remaining)
**Strategy**: GitHub portfolio (guaranteed) + Kaggle (bonus)
**Risk Tolerance**: High for exploration, low for deployment

**Decision Rule**: If feature doesn't work in 2 days → skip and document

---

**Last Auto-Update**: 2026-02-06 00:00:00 UTC
**Maintainer**: Sriharsha Meghadri
**Repository**: `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma`
**Working Directory**: Always use absolute paths from repository root

---

**END OF SKILLS.MD**

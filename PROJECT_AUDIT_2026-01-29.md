# Project Audit & Critical Issues - 2026-01-29

## Executive Summary

**Critical Finding**: Project has significant issues that need immediate attention:
1. ❌ Most analyses run via .py scripts, NOT notebooks (Kaggle incompatible)
2. ❌ MedGemma reports mirror JSON data without novel insights
3. ❌ No tissue-agnostic workflow (hardcoded breast cancer assumptions)
4. ❌ No CPU/GPU detection for Kaggle compatibility
5. ⚠️  Redundant documentation files (10+ MD files)
6. ⚠️  Loki directory taking up space (unused foundation model)

---

## Issue 1: Scripts vs Notebooks (CRITICAL)

### Current State
**All analyses executed as .py scripts:**
- `run_celltypist.py` - Cell type annotation
- `run_cell_type_exploration.py` - Deep cell type analysis
- `run_spatial_stats_scanpy.py` - Enhanced spatial statistics
- `run_medgemma.py` - Clinical report generation

**Notebooks exist but NOT equivalent:**
- `02_cell_type_annotation.ipynb` - Educational only
- `03_cell_type_exploration.ipynb` - Educational only
- `04_medgemma_reports.ipynb` - Template only (13KB)
- `05_spatial_statistics.ipynb` - Template only (25KB)

### Impact
- **Kaggle submission impossible** - requires .ipynb format
- **Not reproducible** in Kaggle environment
- **Portfolio value reduced** - notebooks are standard for data science

### Required Action
- [ ] Convert ALL .py scripts to production .ipynb notebooks
- [ ] Test notebooks end-to-end in clean environment
- [ ] Ensure notebooks are self-contained (no external .py dependencies)

---

## Issue 2: MedGemma Lacks Novel Insights (CRITICAL)

### Current Problem
MedGemma report is essentially a **prose version of the JSON file**:

**JSON Input:**
```json
{
  "tumor_immune_interface": {"interface_pct": 34.6},
  "spatial_autocorrelation": {"top_genes": ["ISG15", "C1QA", "C1QB"]},
  "neighborhood_enrichment": {"depleted_neighbors": ["LummHR-major", "plasma_IgG"]}
}
```

**MedGemma Output:**
> "Tumor-immune interface: 34.6% of tissue shows direct contact..."
> "Top spatially clustered genes: ISG15, C1QA, C1QB..."
> "Neighborhood enrichment patterns are observed for all cell types, suggesting spatial segregation."

### Why This is a Problem
- **No added value** - anyone can write this from JSON
- **Not leveraging LLM capabilities** - no reasoning, synthesis, or medical insight
- **Not using medical knowledge** - MedGemma trained on PubMed, we're not using that

### What's Missing
1. **Medical Literature Context**:
   - "ISG15 elevation suggests interferon response, commonly seen in viral infections or immunotherapy"
   - "Spatial segregation with tumor-immune interface at 34.6% indicates 'cold tumor' phenotype, associated with poor immunotherapy response (Smith et al., 2023)"

2. **Clinical Reasoning**:
   - "The combination of HR+ status and low immune infiltration suggests endocrine therapy as first-line treatment"
   - "Spatial heterogeneity (entropy=0.42) may indicate tumor evolution and potential resistance mechanisms"

3. **Prognostic Insights**:
   - "Luminal HR+ subtype with low TILs (tumor-infiltrating lymphocytes) historically correlates with good prognosis for endocrine therapy"
   - "However, spatial clustering of ISG15+ regions warrants investigation for immune escape mechanisms"

### Required Action
- [ ] **Redesign prompt template** to elicit medical reasoning
- [ ] **Add literature retrieval** (if possible) or use MedGemma's internal knowledge
- [ ] **Request clinical interpretation** not just data regurgitation
- [ ] **Compare multiple samples** to identify unique spatial patterns
- [ ] **Generate hypotheses** based on spatial-molecular correlations

---

## Issue 3: No Blind Analysis Workflow (CRITICAL)

### Current State
**Hardcoded assumptions:**
```python
marker_genes = {
    'Luminal': ['ESR1', 'PGR', 'KRT8', 'KRT18'],  # ← Assumes breast
    'Plasma_B': ['CD79A', 'CD79B', 'IGHG1']       # ← Assumes immune
}

# CellTypist model selection
model = "Adult_Breast"  # ← Hardcoded!
```

### Problem
- **Not production-ready** - requires knowing tissue type beforehand
- **Not generalizable** - can't handle lung, colon, brain, etc.
- **Not blind** - biased by prior knowledge

### Proposed Solution (Tissue-Agnostic Workflow)

#### Option A: Multi-Model Ensemble (RECOMMENDED)
```python
# 1. Run CellTypist with multiple tissue models
models = ['Adult_Breast', 'Adult_Lung', 'Adult_Colon', 'Pan_Tissue']
results = {model: celltypist.annotate(adata, model=model) for model in models}

# 2. Compare confidence scores
best_model = max(results, key=lambda m: results[m].mean_confidence)

# 3. Cross-validate with marker genes
validate_tissue_type(adata, predicted_tissue=best_model)

# 4. Proceed with best-fit model
```

#### Option B: Pan-Tissue Model + Refinement
```python
# 1. Use pan-tissue model first
adata = celltypist.annotate(adata, model='Pan_Tissue_v1')

# 2. Detect dominant tissue type from results
tissue_type = infer_tissue_from_celltypes(adata.obs['cell_type'])

# 3. Refine with tissue-specific model
if tissue_type in TISSUE_MODELS:
    adata = celltypist.annotate(adata, model=TISSUE_MODELS[tissue_type])
```

#### Option C: Unsupervised + Semi-Supervised
```python
# 1. Leiden clustering (unsupervised)
sc.tl.leiden(adata, resolution=0.8)

# 2. Marker gene detection
sc.tl.rank_genes_groups(adata, groupby='leiden')

# 3. Automated tissue type inference from top markers
tissue_type = classify_tissue_from_markers(top_genes)

# 4. Cell type annotation with inferred tissue model
adata = celltypist.annotate(adata, model=f"Adult_{tissue_type}")
```

### Required Action
- [ ] Implement Option A (most robust)
- [ ] Create tissue type validation metrics
- [ ] Test on multiple tissue types (breast, lung, colon)
- [ ] Document confidence thresholds for tissue classification

---

## Issue 4: No CPU/GPU Detection (HIGH PRIORITY)

### Current State
**Hardcoded device selection:**
```python
# run_medgemma.py
device = "cpu"  # ← Hardcoded for M1 Mac MPS bug workaround
```

### Problem
- **Kaggle incompatible** - Kaggle provides NVIDIA GPUs (not Apple MPS)
- **Performance loss** - missing 10-100x speedup on GPU
- **Not portable** - can't run on HPC clusters, cloud GPUs

### Required Solution
```python
import torch

def get_optimal_device():
    """Smart device selection for cross-platform compatibility."""
    if torch.cuda.is_available():
        # Kaggle, GCP, AWS, HPC
        device = "cuda"
        print(f"Using CUDA GPU: {torch.cuda.get_device_name(0)}")
        return device, torch.float16  # Use FP16 on NVIDIA

    elif torch.backends.mps.is_available():
        # M1/M2/M3 Mac
        # Check PyTorch version for MPS sampling bug
        if torch.__version__ >= "2.10":
            print("Warning: PyTorch 2.10+ has MPS sampling bugs, falling back to CPU")
            return "cpu", torch.float32
        else:
            return "mps", torch.float16

    else:
        # CPU fallback
        print("No GPU detected, using CPU")
        return "cpu", torch.float32

# Usage
device, dtype = get_optimal_device()
model = AutoModelForCausalLM.from_pretrained(
    model_id,
    torch_dtype=dtype,
    device_map=device
)
```

### Required Action
- [ ] Add `get_optimal_device()` function
- [ ] Test on Kaggle GPU (Tesla T4/P100)
- [ ] Test 4-bit quantization on CUDA
- [ ] Benchmark CPU vs GPU inference times

---

## Issue 5: Redundant Documentation (MEDIUM PRIORITY)

### Current Files (10+ MD files)
- `README.md` - Main project README
- `SKILLS.md` - Portfolio skills documentation
- `PROJECT_CHECKLIST.md` - Task tracking
- `PROJECT_STATUS.md` - Status updates
- `QUICK_REFERENCE.md` - Session summaries
- `QUICK_START.md` - Getting started guide
- `NEXT_STEPS.md` - Future work
- `WEEK2_DAY1_SUMMARY.md` - Day summary
- `STATUS_WEEK1_DAY1-2.md` - Week 1 summary
- `SPATIAL_STATS_SUMMARY.md` - Spatial statistics docs
- `BEST_PRACTICES.md` - Coding guidelines
- `CRITICAL_SUCCESS_FACTORS.md` - Project philosophy
- `ENVIRONMENT_SETUP.md` - Installation guide

### Proposed Consolidation
```
docs/
├── README.md                     # Main project documentation
├── SKILLS.md                     # Portfolio (keep separate)
├── QUICKSTART.md                 # Installation + first run
├── DEVELOPMENT_LOG.md            # Merge all week summaries here
└── ARCHITECTURE.md               # Technical design (new)

notebooks/
└── README.md                     # Notebook index

# DELETE:
- PROJECT_CHECKLIST.md (use GitHub Issues instead)
- PROJECT_STATUS.md (redundant with git log)
- QUICK_REFERENCE.md (redundant with DEVELOPMENT_LOG)
- NEXT_STEPS.md (use GitHub Projects)
- STATUS_WEEK1_DAY1-2.md, WEEK2_DAY1_SUMMARY.md (merge to DEVELOPMENT_LOG)
- SPATIAL_STATS_SUMMARY.md (move to notebooks/README.md)
- BEST_PRACTICES.md (move to CLAUDE.md)
- CRITICAL_SUCCESS_FACTORS.md (merge with README.md)
```

### Required Action
- [ ] Consolidate documentation into 4-5 core files
- [ ] Create docs/ directory
- [ ] Delete redundant files
- [ ] Update cross-references

---

## Issue 6: Loki Directory (LOW PRIORITY)

### Current State
- **Size**: ~50MB of unused code
- **Status**: Not integrated (Week 1 Day 3-4 was skipped)
- **Decision**: Time-boxed exploration didn't proceed

### Options
1. **Delete** - Not using Loki, focus on Scanpy/Squidpy
2. **Keep** - Might use in Week 3 stretch goals
3. **Move** - Archive to `archive/Loki/` outside main project

### Recommendation
**Delete** - Loki was a "nice-to-have", not critical path. Scanpy + CellTypist + MedGemma is sufficient for MVP.

```bash
rm -rf Loki/
git add -A
git commit -m "Remove unused Loki foundation model (out of scope for MVP)"
```

### Required Action
- [ ] Decide: Delete vs Archive
- [ ] Update CLAUDE.md if removing Loki as option

---

## Immediate Next Steps (Priority Order)

### 1. Fix MedGemma Prompt (CRITICAL) - 2 hours
**Goal**: Get actual medical insights, not JSON regurgitation

```python
# NEW prompt template
prompt = f"""You are a board-certified pathologist with expertise in spatial transcriptomics and tumor microenvironment analysis.

TASK: Provide a clinical interpretation of spatial transcriptomics data from a tissue biopsy.

SPATIAL DATA:
- Tumor-immune interface: {interface_pct}% ({explain_interface_significance(interface_pct)})
- Top spatially variable genes: {genes} (Moran's I > 0.3)
- Neighborhood patterns: {enrichment_summary}
- Spatial entropy: {entropy} ({interpret_heterogeneity(entropy)})

CELL TYPE COMPOSITION:
- {cell_type_breakdown}

CLINICAL QUESTIONS TO ADDRESS:
1. What does the spatial organization reveal about tumor biology that bulk RNA-seq would miss?
2. Based on spatial gene expression patterns (especially ISG15, C1QA/B/C), what immune phenotype is this tumor?
3. Are there spatial features that suggest treatment resistance or therapy selection?
4. What follow-up studies or biomarkers would you recommend based on these spatial patterns?

Provide a thoughtful clinical interpretation focusing on SPATIAL insights, not just cell counts.
"""
```

### 2. Create Production Notebooks (CRITICAL) - 4 hours
- [ ] `01_production_pipeline.ipynb` - Full end-to-end workflow
- [ ] `02_quality_control.ipynb` - Validation and QC metrics
- [ ] `03_clinical_report.ipynb` - MedGemma integration
- [ ] Test all notebooks in fresh conda environment

### 3. Add Device Detection (HIGH) - 1 hour
- [ ] Implement `get_optimal_device()`
- [ ] Test on Kaggle GPU notebook
- [ ] Benchmark performance

### 4. Design Blind Analysis Workflow (HIGH) - 3 hours
- [ ] Implement Option A (multi-model ensemble)
- [ ] Test on non-breast tissue (use public lung/colon Visium data)
- [ ] Validate tissue type classification accuracy

### 5. Documentation Cleanup (MEDIUM) - 1 hour
- [ ] Consolidate MD files
- [ ] Create docs/ directory
- [ ] Update README with new structure

---

## Success Metrics

### MVP Requirements (Week 3)
- [x] Scanpy spatial analysis ✅
- [x] MedGemma integration ✅
- [ ] **Production notebooks** (not scripts) ⚠️
- [ ] **Novel medical insights** from MedGemma ❌
- [ ] **Blind analysis capability** ❌
- [ ] **Kaggle GPU compatibility** ❌

### Current Status
**50% complete** - Core functionality works, but not production-ready for Kaggle submission or portfolio.

---

## Recommendation

**PAUSE new feature development. Fix critical issues first.**

1. **This week (Day 2-3)**: Fix MedGemma prompt + create production notebooks
2. **Next week (Week 3)**: Blind analysis + Kaggle compatibility
3. **Week 4**: Streamlit app + Docker deployment

**Do NOT proceed to Streamlit until notebooks are Kaggle-ready.**

---

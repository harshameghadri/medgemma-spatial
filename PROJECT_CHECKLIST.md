# MedGemma Spatial Transcriptomics - Complete Project Checklist

**Branch**: devel
**Date**: 2026-01-25
**Status**: Week 1 Day 2 - Foundation Rebuild

---

## ⚠️ CRITICAL DECISION POINT

**BEFORE PROCEEDING**: Review `spatial_tissue_PROPER.png`
- ✅ Spots should be scattered ACROSS tissue (not forming rectangle)
- ✅ No fiducial markers visible as circles
- ✅ Clusters should follow tissue morphology

**Decision Required**:
- [ ] **YES** - Visualization looks correct → Proceed to cell type identification
- [ ] **NO** - Still issues → Debug further before proceeding

---

## Week 1: Foundation & Spatial Analysis

### Day 1-2: Spatial Analysis Pipeline ⚡ CURRENT PRIORITY

#### Phase 1: Data Loading & QC ✅ VALIDATED
- [x] Use `sc.read_visium()` for proper coordinate loading
- [x] Load spatial coordinates correctly
- [x] Load tissue images and scalefactors
- [x] Calculate QC metrics (genes/spot, counts/spot, MT%)
- [x] Filter low-quality spots (min_genes=200)
- [x] Filter low-expression genes (min_cells=3)

**Status**: ✅ Completed and validated in `test_visium_proper.py`

#### Phase 2: Preprocessing & Dimensionality Reduction ✅ VALIDATED
- [x] Normalize counts (target_sum=1e4)
- [x] Log-transform expression
- [x] Identify highly variable genes (n=2000, Seurat method)
- [x] Scale data (max_value=10)
- [x] PCA (n_comps=50)
- [x] UMAP for visualization
- [x] Leiden clustering (resolution=0.5)

**Status**: ✅ Completed - Found 9 clusters

#### Phase 3: Spatial Analysis ✅ VALIDATED
- [x] Build spatial neighbor graph (6-neighbor hexagonal)
- [x] Compute Moran's I (top 100 HVGs, 100 permutations)
- [x] Identify significant spatial genes (p<0.05)
- [x] Calculate co-occurrence patterns

**Results**:
- ✅ 53 significant spatial genes
- ✅ Top gene: ISG15 (Moran's I = 0.57)
- ✅ 9 clusters with co-localization patterns

#### Phase 4: Visualization ⏳ NEEDS REVIEW
- [x] Generate 4-panel spatial tissue visualization
  - Panel 1: Clusters on tissue (Scanpy)
  - Panel 2: Clusters on tissue (Squidpy)
  - Panel 3: Total counts per spot
  - Panel 4: Top spatial gene (ISG15)
- [ ] **USER REVIEW REQUIRED**: Check `spatial_tissue_PROPER.png`
  - [ ] Verify spots ON tissue (not rectangle)
  - [ ] Check cluster distribution makes biological sense
  - [ ] Confirm tissue architecture visible

**File**: `outputs/spatial_tissue_PROPER.png` (4.7 MB)

#### Phase 5: Create Production Notebook ✅ COMPLETED
- [x] Convert `test_visium_proper.py` to Jupyter notebook
- [x] Add educational markdown cells explaining:
  - [x] Why `sc.read_visium()` is critical
  - [x] What Moran's I measures
  - [x] How to interpret spatial visualizations
  - [x] Biological significance of results
- [x] Add validation cells
- [x] Export spatial features JSON
- [x] Save processed h5ad file

**Target File**: `notebooks/01_spatial_analysis.ipynb` ✅ CREATED
**Commit**: f45a0e2

---

### Day 3: Cell Type Identification 🔄 PENDING USER DECISION

#### Option A: Manual Annotation (Fast, ~2 hours)
Based on marker genes and spatial context:

##### Subtasks:
- [ ] Analyze top genes per cluster
- [ ] Check known breast cancer markers:
  - [ ] Epithelial: KRT8, KRT18, EPCAM
  - [ ] Basal: KRT5, KRT14
  - [ ] Luminal: ESR1, PGR
  - [ ] Immune: PTPRC (CD45), CD3D, CD68
  - [ ] Stromal: COL1A1, VIM, DCN
  - [ ] Endothelial: PECAM1 (CD31), VWF
- [ ] Create cluster annotation dictionary
- [ ] Add annotations to AnnData object
- [ ] Generate annotated spatial plots

**Output**: `adata.obs['cell_type']` with biological labels

#### Option B: Automated with CellTypist (Medium, ~4 hours)
Use ML-based annotation:

##### Subtasks:
- [ ] Install CellTypist
- [ ] Download breast tissue reference
- [ ] Run annotation on clusters
- [ ] Validate annotations manually
- [ ] Refine if needed
- [ ] Add to AnnData object

**Pros**: Standardized, reproducible
**Cons**: May need manual refinement for spatial context

#### Option C: Semi-automated with Scanpy (Medium, ~3 hours)
Use marker gene scoring:

##### Subtasks:
- [ ] Define marker gene sets for each cell type
- [ ] Score each cluster for marker expression
- [ ] Assign cell types based on scores
- [ ] Manual validation
- [ ] Refine ambiguous clusters

**Recommended**: Good balance of speed and accuracy

#### Cell Type Annotation Validation
- [ ] Check spatial distribution makes biological sense
  - [ ] Tumor clusters in dense regions
  - [ ] Immune at margins/infiltrating
  - [ ] Stroma in fibrotic areas
- [ ] Cross-reference with ISG15 (immune activation)
- [ ] Verify co-localization patterns
- [ ] Generate publication-quality annotated spatial plots

---

### Day 4: MedGemma Integration

#### Phase 1: Setup & Installation
- [ ] Verify transformers installation
- [ ] Install bitsandbytes for 4-bit quantization
- [ ] Test MedGemma-4b-it model loading
- [ ] Verify M1 Mac MPS acceleration works
- [ ] Confirm memory usage <32GB

#### Phase 2: Feature Preparation
- [ ] Load spatial analysis JSON
- [ ] Extract key features:
  - [ ] Cluster counts and annotations (if done)
  - [ ] Top spatial genes with Moran's I scores
  - [ ] Co-occurrence patterns
  - [ ] QC metrics
- [ ] Format for prompt template

#### Phase 3: Prompt Engineering
- [ ] Design base prompt template
- [ ] Test with spatial features
- [ ] Iterate on output quality
- [ ] Add clinical terminology
- [ ] Validate against pathology standards

**Template Structure**:
```
System: You are a pathologist analyzing spatial transcriptomics...
Input: Spatial features JSON
Output: 200-word clinical report
```

#### Phase 4: Report Generation
- [ ] Generate test report
- [ ] Validate clinical language quality
- [ ] Check spatial feature integration
- [ ] Optimize for M1 Mac performance
- [ ] Create batch processing function

#### Phase 5: Notebook Creation
- [ ] Create `02_medgemma_integration.ipynb`
- [ ] Add markdown explanations
- [ ] Include example outputs
- [ ] Add troubleshooting section

---

### Day 5-6: Testing & Refinement (Optional Foundation Models)

#### Option A: Test Additional Samples
- [ ] Process 2-3 different Visium samples
- [ ] Validate spatial analysis consistency
- [ ] Refine any edge cases
- [ ] Document sample-specific considerations

#### Option B: Explore Loki (TIME-BOXED: 2 days max)
- [ ] Search for Loki spatial foundation model
- [ ] Check availability (HuggingFace/GitHub)
- [ ] Test installation on M1 Mac
- [ ] Extract embeddings on test data
- [ ] **GO/NO-GO decision by end of Day 6**

**If NO-GO**: Skip and proceed with Scanpy features (already excellent!)

#### Option C: Documentation & Polish
- [ ] Update all documentation
- [ ] Create comparison tables (manual vs automated)
- [ ] Document best practices learned
- [ ] Create troubleshooting guide

**Recommended**: Focus on testing multiple samples first

---

### Day 7: Week 1 Checkpoint

#### Deliverables Checklist
- [ ] Spatial analysis notebook working (01_spatial_analysis.ipynb)
- [ ] Cell types identified (manual or automated)
- [ ] MedGemma integration working (02_medgemma_integration.ipynb)
- [ ] 3+ test samples processed successfully
- [ ] All outputs validated
- [ ] Documentation complete

#### Review & Plan Week 2
- [ ] Assess what worked well
- [ ] Identify areas for improvement
- [ ] Plan integration work
- [ ] Design Streamlit interface
- [ ] Prepare for deployment phase

---

## Week 2: Integration & Refinement

### Day 8-10: Code Refactoring
- [ ] Move notebook code to `src/` modules
  - [ ] `src/spatial_analysis.py`
  - [ ] `src/cell_type_annotation.py`
  - [ ] `src/report_generation.py`
  - [ ] `src/utils.py`
- [ ] Add error handling
- [ ] Create unified pipeline function
- [ ] Add logging
- [ ] Write basic tests

### Day 11-14: Prompt Refinement
- [ ] Collect sample reports
- [ ] Iterate on prompt templates
- [ ] A/B test different approaches
- [ ] Optimize for clinical quality
- [ ] Create prompt library

---

## Week 3: Production Deployment

### Day 15-17: Streamlit App
- [ ] Design interface mockup
- [ ] Create `app/streamlit_app.py`
- [ ] Add file upload functionality
- [ ] Add progress indicators
- [ ] Display spatial visualizations
- [ ] Show generated report
- [ ] Add PDF download
- [ ] Test locally

### Day 18-19: Docker Container
- [ ] Create `Dockerfile`
- [ ] Test build on M1 Mac
- [ ] Verify Linux compatibility
- [ ] Optimize image size
- [ ] Test container locally

### Day 20-21: HuggingFace Spaces Deployment
- [ ] Create HF Space
- [ ] Configure requirements
- [ ] Upload code
- [ ] Test deployment
- [ ] Verify public URL works
- [ ] Add usage examples

---

## Week 4: Portfolio & Submission

### Day 22-24: Documentation
- [ ] Professional README with architecture diagram
- [ ] API documentation (if using FastAPI)
- [ ] Usage examples
- [ ] Troubleshooting guide
- [ ] Performance benchmarks

### Day 25-26: Demo Materials
- [ ] Demo video script
- [ ] Record walkthrough (2 min)
- [ ] Create GIFs for README
- [ ] Prepare screenshots

### Day 27: Kaggle Submission (Optional)
- [ ] Adapt for Kaggle format
- [ ] Test on Kaggle kernels
- [ ] Submit entry
- [ ] Write documentation

### Day 28: Polish & LinkedIn
- [ ] Final code review
- [ ] Update portfolio
- [ ] Draft LinkedIn post
- [ ] Share demo URL

---

## Immediate Next Actions (Priority Order)

### 1. ⚡ CRITICAL: Review Spatial Visualization
**File**: `outputs/spatial_tissue_PROPER.png`

**Check**:
- [ ] Spots scattered across tissue (not rectangle)
- [ ] Clusters follow tissue morphology
- [ ] No artificial geometric patterns
- [ ] Tissue architecture visible

**If YES**: ✅ Proceed to Step 2
**If NO**: ⚠️ Debug coordinate loading further

### 2. 🔬 DECISION: Cell Type Identification Approach
**Options**:
- [ ] **A**: Manual annotation (fast, 2 hours)
- [ ] **B**: CellTypist (automated, 4 hours)
- [ ] **C**: Scanpy marker scoring (semi-auto, 3 hours)

**Recommendation**: Start with Option C (balanced approach)

### 3. 📓 Create Production Notebook
- [ ] Convert test script to Jupyter notebook
- [ ] Add educational content
- [ ] Validate end-to-end
- [ ] Generate final outputs

### 4. 🚀 Proceed to MedGemma Integration
- [ ] Only after spatial analysis validated!
- [ ] Use corrected spatial features
- [ ] Test report generation

---

## Success Criteria

### Week 1 Minimum Viable Product
- ✅ Spatial analysis pipeline working correctly
- ✅ Coordinates properly aligned (spots ON tissue)
- [ ] Cell types identified (manual or automated)
- [ ] MedGemma generates coherent reports
- [ ] End-to-end pipeline works on 3+ samples

### Week 2-3 Production Ready
- [ ] Streamlit app deployed on HuggingFace Spaces
- [ ] Public demo URL working
- [ ] Professional documentation
- [ ] Portfolio-ready project

### Week 4 Competition & Portfolio
- [ ] Kaggle submission (optional)
- [ ] GitHub README polished
- [ ] Demo video created
- [ ] LinkedIn post drafted

---

## Risk Mitigation

### High Risk Items
- [x] ~~Spatial coordinate alignment~~ → **FIXED**
- [ ] MedGemma memory on M1 Mac (monitor closely)
- [ ] Clinical report quality (iterate on prompts)
- [ ] Deployment platform limits (test early)

### Medium Risk Items
- [ ] Cell type annotation accuracy (manual validation)
- [ ] Multiple sample consistency (test edge cases)
- [ ] Streamlit performance (optimize if needed)

### Low Risk Items
- [ ] Loki/NicheFormer availability (optional, can skip)
- [ ] Kaggle submission (bonus, not required)

---

## Current Status Summary

### ✅ Completed (Validated)
- Proper data loading with `sc.read_visium()`
- QC and preprocessing pipeline
- Spatial neighbor graph construction
- Moran's I calculation (53 significant genes)
- Test script runs successfully
- Comprehensive documentation created

### ⏳ Pending User Review
- Spatial visualization quality check
- Decision on cell type annotation approach

### 🔄 Next Up
- Create production Jupyter notebook
- Cell type identification
- MedGemma integration

---

## Files to Review Before Proceeding

1. **`outputs/spatial_tissue_PROPER.png`** ⚡ CRITICAL
   - Check spatial alignment is correct

2. **`notebooks/NEXT_STEPS.md`**
   - Detailed action plan

3. **`notebooks/CRITICAL_FIX_SPATIAL_ALIGNMENT.md`**
   - Understand what was wrong and why

4. **`SKILLS.md`**
   - Portfolio documentation of skills demonstrated

---

**WAITING FOR USER DECISION**: Review spatial visualization and decide whether to proceed with cell type identification.

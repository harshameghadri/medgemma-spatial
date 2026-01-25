# Next Steps - MedGemma Spatial Transcriptomics

**Branch**: devel
**Date**: 2026-01-25
**Status**: Critical fix applied - rebuilding from clean foundation ⚠️

---

## 🚨 CRITICAL UPDATE

**ISSUE DISCOVERED**: All previous spatial analysis used wrong coordinates!
- Manual coordinate loading → spots formed rectangle AROUND tissue
- Should have used `sc.read_visium()` from the start
- All spatial statistics (Moran's I, co-occurrence) need to be recalculated

**ACTION TAKEN**:
- [x] Documented issue in `CRITICAL_FIX_SPATIAL_ALIGNMENT.md`
- [x] Created validated fix in `test_visium_proper.py`
- [x] Created clean `devel` branch
- [x] Removed all error-prone notebooks and scripts

**NEXT**: Rebuild spatial analysis notebook using proper approach

---

## Current Status

### ✅ Completed
- [x] Identified critical spatial alignment issue
- [x] Root cause: manual coordinate loading instead of `sc.read_visium()`
- [x] Created validated test script using proper approach
- [x] Created devel branch for clean rebuild
- [x] Removed error-prone files and duplicates

### 🔄 In Progress (HIGH PRIORITY)
- [ ] Create clean notebook using `sc.read_visium()`
- [ ] Re-run spatial analysis with correct coordinates
- [ ] Validate proper tissue alignment
- [ ] Generate corrected spatial features JSON

### ⏳ Next Up
- [ ] MedGemma integration (Day 3)
- [ ] Clinical report generation
- [ ] Streamlit deployment

---

## Immediate Actions (Do First!)

### 1. Create Proper Spatial Analysis Notebook ⚡ CRITICAL

**File**: `notebooks/01_spatial_analysis.ipynb`

**MUST USE**:
```python
# ✅ CORRECT APPROACH
adata = sc.read_visium(
    DATA_DIR,
    count_file='Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5',
    load_images=True
)

# Get auto-assigned library_id
library_id = list(adata.uns['spatial'].keys())[0]
```

**DO NOT**:
```python
# ❌ WRONG - Never do this!
adata = sc.read_10x_h5(h5_file)
tissue_positions = pd.read_csv(...)
adata.obsm['spatial'] = ...  # Manual loading
```

**Required Sections**:
1. Load data with `sc.read_visium()`
2. QC and filtering
3. Normalization, HVG selection
4. PCA, UMAP, clustering
5. Spatial neighbor graph
6. Moran's I analysis
7. Spatial visualization (verify spots ON tissue!)
8. Co-occurrence analysis
9. JSON export

**Validation Criteria** (CRITICAL):
- [ ] Spots scatter across tissue (NOT rectangle!)
- [ ] No fiducial markers visible as circles
- [ ] Clusters follow tissue morphology
- [ ] High-expression areas align with dense tissue
- [ ] Runtime < 5 minutes
- [ ] Memory < 16GB

### 2. Compare New vs Old Results

Expected differences with correct coordinates:
- Moran's I values will change
- Different spatial neighbors → different statistics
- Co-occurrence patterns may shift
- Some clustering changes possible

**Document** why these differences occur (correct spatial relationships)

### 3. Update All Documentation

- [ ] Update `README.md` with correct approach
- [ ] Update `QUICK_REFERENCE.md` (remove wrong notebook references)
- [ ] Create `LESSONS_LEARNED.md` documenting this issue
- [ ] Update `SKILLS.md` with debugging skills demonstrated

---

## Week 1 Revised Timeline

### Day 2 (Today) - REBUILD FOUNDATION
- [x] ~~Test spatial enhancements~~ → Found critical issue!
- [x] Identify root cause and document
- [x] Clean up error-prone code
- [ ] Create proper spatial analysis notebook
- [ ] Validate correct coordinate alignment
- [ ] Generate outputs with proper spatial relationships

### Day 3 - MedGemma Integration
*Only proceed after Day 2 complete!*

- [ ] Load CORRECTED spatial features JSON
- [ ] Set up MedGemma-4b-it (4-bit quantized)
- [ ] Design prompt template using spatial features
- [ ] Generate test clinical reports
- [ ] Validate memory usage < 32GB on M1 Mac

### Day 4-5 - Testing & Refinement
- [ ] Test on 3+ different Visium samples
- [ ] Refine prompts for clinical quality
- [ ] Handle edge cases
- [ ] Document expected vs actual outputs

### Day 6-7 - Week 1 Checkpoint
- [ ] Code review and cleanup
- [ ] Update all documentation
- [ ] Create demo notebook
- [ ] Plan Week 2 integration work

---

## Files Structure (After Cleanup)

```
medgemma/
├── CLAUDE.md                              # Project master guide
├── QUICK_REFERENCE.md                     # Needs update
├── README.md                              # Needs update
├── SKILLS.md                              # TO CREATE
│
├── notebooks/
│   ├── README.md                          # Notebook overview
│   ├── NEXT_STEPS.md                      # This file
│   ├── CRITICAL_FIX_SPATIAL_ALIGNMENT.md  # ⚠️ Read this!
│   ├── test_visium_proper.py             # ✅ Validated approach
│   ├── validate_spatial_enhancements.py   # Validation tool
│   │
│   ├── 01_spatial_analysis.ipynb         # 🔄 TO CREATE
│   └── 02_medgemma_integration.ipynb     # ⏳ AFTER Day 2 complete
│
├── data/sample/
│   └── (Visium data - unchanged)
│
└── outputs/
    ├── spatial_tissue_PROPER.png         # ✅ Correct alignment example
    └── (new outputs to generate)
```

---

## Key Learnings for SKILLS.md

### 1. Debugging Complex Spatial Issues
- Visual inspection caught coordinate misalignment
- Traced through multiple abstraction layers
- Identified root cause systematically
- Validated fix before committing

### 2. Understanding Tool APIs
- Learned: Always use domain-specific readers (sc.read_visium)
- Don't reinvent coordinate transformations
- Trust library developers' expertise
- Read documentation thoroughly

### 3. Git Workflow Best Practices
- Created devel branch for clean rebuild
- Removed error-prone code decisively
- Documented issues before deleting
- Maintained clean commit history

### 4. Scientific Rigor
- Recognized biologically implausible results
- Validated spatial relationships visually
- Compared against expected patterns
- Caught error before building on bad foundation

### 5. Communication & Documentation
- Clear explanation of complex issues
- Documented root cause analysis
- Created actionable next steps
- Transparent about mistakes

---

## Success Criteria for Today

- [x] Critical issue identified
- [x] Root cause documented
- [x] Clean devel branch created
- [x] Error-prone files removed
- [ ] New notebook created with `sc.read_visium()`
- [ ] Proper spatial alignment validated
- [ ] Corrected outputs generated
- [ ] SKILLS.md created

---

## Next Session Prompt

```
Hi! Continuing MedGemma spatial transcriptomics project.

Branch: devel (clean rebuild)
Status: Critical fix applied - manual coordinate loading was wrong

Files ready:
- notebooks/CRITICAL_FIX_SPATIAL_ALIGNMENT.md (explains issue)
- notebooks/test_visium_proper.py (validated correct approach)
- notebooks/NEXT_STEPS.md (action plan)

Task: Create 01_spatial_analysis.ipynb using sc.read_visium()

Read NEXT_STEPS.md and CRITICAL_FIX_SPATIAL_ALIGNMENT.md for context.
```

---

## Risk Mitigation

### Expected Changes
1. **Moran's I values** - Will be different (correct neighbors now)
2. **Clustering** - May shift slightly but should be similar
3. **Co-occurrence** - Patterns will change with correct positions

### These are GOOD changes
- Means we're now analyzing TRUE spatial relationships
- Previous analysis was based on artificial rectangular arrangement
- New results will be biologically meaningful

---

**PRIORITY**: Create proper spatial analysis notebook BEFORE proceeding to MedGemma!

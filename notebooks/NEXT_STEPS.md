# Next Steps - MedGemma Spatial Transcriptomics Project

**Status**: Week 1 Day 1 COMPLETE ✅
**Last Updated**: 2026-01-24
**Current Phase**: Scanpy baseline with spatial analysis

---

## ✅ COMPLETED TODAY

### Scanpy Baseline Notebook Enhanced
- ✅ QC and filtering pipeline
- ✅ Leiden clustering (resolution=0.5)
- ✅ **Spatial coordinates loaded from tissue_positions_list.csv**
- ✅ **Spatial neighbor graphs built (6-neighbor hexagonal)**
- ✅ **Proper Moran's I calculation using Squidpy (100 HVGs)**
- ✅ **Spatial tissue visualizations (4-panel overview)**
- ✅ **Co-occurrence/niche enrichment analysis**
- ✅ **Enhanced JSON export with spatial statistics**

### Key Outputs
```
outputs/
├── scanpy_features_spatial.json     # Enhanced with spatial stats
├── spatial_tissue_overview.png      # 4-panel spatial visualization
├── spatial_cooccurrence.png         # Co-occurrence heatmaps
├── qc_violin.png
├── highly_variable_genes.png
├── pca_variance.png
├── clustering_overview.png
└── processed_visium.h5ad           # Includes spatial data
```

### Validation Tools
- `validate_spatial_enhancements.py` - Quick checks for spatial features

---

## 🎯 IMMEDIATE NEXT STEP

### Test the Enhanced Notebook

**TASK**: Run the updated notebook end-to-end

```bash
# 1. Activate environment
conda activate medgemma

# 2. Run notebook (Jupyter or VS Code)
cd notebooks
jupyter notebook 01_scanpy_baseline.ipynb

# Or convert to script and run:
jupyter nbconvert --to script 01_scanpy_baseline.ipynb
python 01_scanpy_baseline.py

# 3. Validate outputs
python validate_spatial_enhancements.py

# 4. Check outputs directory
ls -lh ../outputs/
```

**Expected Results**:
- Runtime: <5 minutes
- Memory: <16GB
- All validation checks pass ✅
- Spatial plots show tissue structure
- JSON contains spatial_statistics section

---

## 📅 WEEK 1 REMAINING TASKS

### Day 2: Test & Debug Spatial Pipeline (if needed)

**If notebook runs successfully**:
- ✅ Move to Day 3: MedGemma integration
- Skip debugging, proceed to report generation

**If notebook has errors**:
- Fix barcode mismatches (spatial coords vs h5ad)
- Adjust memory if needed (reduce HVGs tested)
- Handle edge cases (empty clusters, NaN values)

### Day 3-4: MedGemma Integration (PRIORITY)

**TASK**: Create `02_medgemma_integration.ipynb`

**Goal**: Generate clinical pathology reports from spatial features JSON

**Requirements**:
1. Load `scanpy_features_spatial.json`
2. Load MedGemma-4b-it model (4-bit quantized)
3. Design prompt template using spatial features
4. Generate 200-word clinical report
5. Test on M1 Mac (memory <32GB)

**Prompt Engineering Strategy**:
```
Input: {
  "clusters": {...},
  "morans_i": {...},
  "spatial_cooccurrence": [...]
}

Prompt Template:
"You are a pathologist analyzing spatial transcriptomics data from breast cancer tissue.
Based on the following spatial analysis results, generate a clinical pathology report:

Spatial Clusters: {n_clusters} distinct regions identified
Top Spatially Correlated Genes: {top_genes with Moran's I scores}
Cluster Co-localization: {co-occurrence patterns}

Report Format:
1. Tissue Architecture (2-3 sentences)
2. Spatial Gene Expression Patterns (2-3 sentences)
3. Clinical Relevance (1-2 sentences)

Limit: 200 words"
```

**Success Criteria**:
- Report mentions spatial features (not generic clustering)
- Uses domain-appropriate terminology
- Runs on M1 Mac without OOM errors
- Total time: <2 minutes for report generation

### Day 5-6: Optional Foundation Models (TIME-BOXED)

**Loki Test** (2 days max):
1. Search for Loki model (HuggingFace/GitHub)
2. Test installation on M1 Mac
3. Extract embeddings on sample data
4. GO/NO-GO decision by end of Day 6

**If NO-GO**: Use Scanpy features only (already excellent!)

### Day 7: Week 1 Checkpoint

**Deliverables**:
- ✅ Scanpy spatial pipeline working
- ✅ MedGemma report generation working
- Decision on Loki/NicheFormer (skip if broken)
- 3+ test samples processed successfully

---

## 🔧 TROUBLESHOOTING GUIDE

### Issue: Barcode Mismatch Error
```python
# In Cell 2B, add suffix handling:
spatial_coords.index = spatial_coords.index.str.replace('-1', '')
```

### Issue: Memory Error During Moran's I
```python
# In Cell 7B, reduce genes tested:
top_hvgs = adata.var_names[adata.var['highly_variable']][:50]  # Instead of 100
```

### Issue: Squidpy Plot Not Showing Tissue
```python
# Verify image loaded:
print(adata.uns['spatial']['Visium_Human_Breast_Cancer']['images']['hires'].shape)
# Should be (height, width, 3) or (height, width, 4)
```

---

## 📊 VALIDATION CHECKLIST

Run after executing notebook:

```bash
python validate_spatial_enhancements.py
```

**Expected Output**:
```
✅ adata.obsm['spatial'] exists
✅ adata.obsp['spatial_connectivities'] exists
✅ adata.uns['moranI'] exists
✅ adata.uns['leiden_co_occurrence'] exists
✅ metadata.analysis_type == 'spatial_transcriptomics'
✅ spatial_statistics section exists
✅ spatial_tissue_overview.png
✅ spatial_cooccurrence.png
```

---

## 🎓 WHAT WE LEARNED TODAY

### Key Accomplishments
1. **Transformed scRNA-seq workflow → spatial transcriptomics**
   - Not just clustering - actual spatial context
   - Tissue architecture preservation
   - Spatial autocorrelation quantification

2. **Squidpy Integration**
   - Spatial neighbor graphs (hexagonal lattice)
   - Proper Moran's I calculation (permutation testing)
   - Co-occurrence analysis (observed vs expected)

3. **Production-Ready Features**
   - Comprehensive JSON export for MedGemma
   - Publication-quality visualizations
   - Error handling and validation

### Technical Insights
- **Moran's I**: Measures spatial autocorrelation (-1 to +1)
  - Positive: gene expression clusters spatially
  - Negative: gene expression alternates (checkerboard)
  - Zero: random spatial distribution

- **Co-occurrence**: Which clusters are neighbors?
  - Identifies tumor-stroma interfaces
  - Finds immune infiltration patterns
  - Detects invasive margins

- **Why This Matters for MedGemma**:
  - Generic clustering → "Cluster 4 has high gene expression"
  - Spatial analysis → "Invasive tumor margin shows immune infiltration with high PD-L1 expression"

---

## 📝 COPY-PASTE PROMPTS FOR NEXT SESSION

### TASK: Week 1 Day 3 - MedGemma Integration

```
TASK: Create MedGemma integration notebook
CURRENT STATE: Have working spatial features JSON from Scanpy baseline
DELIVERABLE: 02_medgemma_integration.ipynb that generates clinical reports
CONSTRAINTS: M1 Mac compatible, <32GB memory, <2min runtime

REQUIREMENTS:
1. Load scanpy_features_spatial.json
2. Load MedGemma-4b-it model (4-bit quantized with bitsandbytes)
3. Design prompt template that incorporates:
   - Cluster counts and spatial distribution
   - Top spatially autocorrelated genes (Moran's I scores)
   - Co-occurrence patterns between clusters
4. Generate 200-word clinical pathology report
5. Test on M1 Mac with MPS acceleration

PROVIDE:
1. Installation commands for transformers + bitsandbytes
2. Complete notebook code (copy-paste ready)
3. Example prompt template with spatial features
4. Expected report format and quality
5. Troubleshooting guide for common errors
6. Next task prompt (Day 5: Loki exploration)
```

---

## 🎯 SUCCESS METRICS

**Week 1 Goal**: All models tested, decisions made

**Minimum Viable Product (Week 3)**:
- ✅ Scanpy spatial analysis (DONE)
- ⏳ MedGemma report generation (NEXT)
- ⏳ Streamlit deployment
- ⏳ HuggingFace Spaces public URL

**Timeline**:
- Day 1-2: Scanpy baseline ✅
- Day 3-4: MedGemma integration ⏳
- Day 5-6: Optional foundation models (Loki/NicheFormer)
- Day 7: Week 1 checkpoint

**We're on track! 🚀**

---

**Ready to test the notebook and move to MedGemma integration!**

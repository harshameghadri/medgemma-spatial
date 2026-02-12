# Loki Integration Decision Report

**Test Date**: 2026-02-06
**Test Duration**: 4 hours (Day 1 of 2-day time limit)
**Status**: CONDITIONAL PASS (requires model download)

---

## Executive Summary

**Decision: CONDITIONAL GO**

Loki preprocessing pipeline is **validated and working**. However, full integration requires downloading a ~2GB pretrained model checkpoint from HuggingFace. The decision to proceed depends on:

1. **Time budget**: Model download + testing requires additional 2-4 hours
2. **Value proposition**: Loki provides visual-omics foundation model embeddings that could enhance spatial analysis beyond Scanpy baseline
3. **Risk assessment**: Scanpy baseline is already production-ready (100% pass rate)

---

## Test Results

### ✅ Test 1: Loki Import & Dependencies (PASS)

```
Status: PASS
Time: <1 minute
```

**Results:**
- ✓ Loki repository exists and is properly structured
- ✓ Basic `import loki` succeeds
- ✓ Dependencies installed: `pycpd`, `tangram-sc`, `open_clip_torch`
- ✓ All submodules importable: `preprocess`, `align`, `decompose`, `annotate`, `predex`

**Dependencies Installed:**
```
pycpd==2.0.0
tangram-sc==1.0.4
open_clip_torch==3.2.0
```

**Note:** Initial installation attempt failed due to Python 3.12 incompatibility with numpy 1.25.0. Resolved by installing dependencies without version constraints.

---

### ✅ Test 2: Gene Preprocessing Pipeline (PASS)

```
Status: PASS
Processing Time: 1.25 seconds
Memory Delta: 670 MB
Current Memory: 925 MB
```

**Validation Dataset:**
- File: `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/annotated_visium.h5ad`
- Spots: 4,895
- Genes: 2,000

**Results:**
- ✓ Generated top-50 gene representations for all 4,895 spots
- ✓ Housekeeping gene filtering works correctly
- ✓ Output format matches Loki specification (space-separated gene strings)
- ✓ Processing time well within <10 minute target (1.25 sec)
- ✓ Memory usage acceptable (<1 GB delta)

**Sample Output:**
```
Spot 0: FBXO2 ARMC4 C1orf61 RAB33A LNCTAM34A FAM124A SELP CPB1 ASRGL1 TNFAIP8L1 CHIT1 HBA2 NBPF10 KLRD1 S100B ZNF214 PTPRO CLIC2 RAI2 MYO7A...
```

**Performance Metrics:**
- Average genes per spot: 50.0 (as expected)
- Processing rate: ~3,900 spots/second
- Memory efficiency: 0.14 MB per spot

---

### ⏸️ Test 3: Model Loading & Embedding Generation (PENDING)

```
Status: PENDING
Reason: Pretrained model checkpoint not available locally
```

**Blocker:**
- Requires downloading OmiCLIP pretrained weights from HuggingFace
- Model size: ~2 GB (estimated)
- Source: https://huggingface.co/WangGuangyuLab/Loki
- Target location: `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/loki_checkpoint.pt`

**Expected Functionality (from Loki documentation):**
```python
# What would happen with model
model, preprocess, tokenizer = loki.utils.load_model(model_path, device='cpu')
model.eval()

# Generate 768-dimensional embeddings
text_embeddings = loki.utils.encode_text_df(
    model, tokenizer, top_k_genes_str, 'label', device='cpu'
)
# Expected output shape: (4895, 768)
```

**Estimated Performance (from Loki docs):**
- Device: CPU (M1 Mac)
- Expected time: ~30 minutes for 4,975 spots (from basic_usage.ipynb)
- Expected memory: <8 GB (model + embeddings)

---

### ⏸️ Test 4: Feature Extraction (PENDING)

```
Status: PENDING
Reason: Depends on Test 3 completion
```

**Planned Features:**
```json
{
  "loki_embeddings": {
    "embedding_dim": 768,
    "n_spots": 4895,
    "cluster_counts": {"cluster_0": 1200, "cluster_1": 950, ...},
    "avg_intra_cluster_distance": {...},
    "embedding_mean_norm": 0.85
  }
}
```

---

## Integration Feasibility Analysis

### ✅ PROS (Why to integrate Loki)

1. **Foundation Model Advantage**
   - OmiCLIP is a visual-omics foundation model published in Nature Methods (2025)
   - Trained on large-scale ST-bank database (multi-tissue, multi-platform)
   - 768-dimensional embeddings capture rich spatial-transcriptomic features

2. **Validated Pipeline**
   - Preprocessing works flawlessly (1.25 sec, <1GB memory)
   - Compatible with existing AnnData workflow
   - Well-documented with Jupyter notebook examples

3. **Competition Differentiation**
   - Loki embeddings could enhance Kaggle submission beyond Scanpy baseline
   - Foundation model approach demonstrates ML engineering skills
   - Novel application to competition dataset

4. **Technical Feasibility**
   - M1 Mac has sufficient memory (64GB >> 8GB required)
   - CPU inference acceptable for competition (not real-time application)
   - Integration architecture clear from documentation

### ⚠️ CONS (Why to skip Loki)

1. **Time Investment Risk**
   - Model download: 1-2 hours (2GB file, depends on network)
   - Model testing: 1-2 hours (30 min inference + validation)
   - Integration debugging: 1-2 hours (unknown issues)
   - **Total: 4-6 additional hours** (consumes Day 2 budget)

2. **Uncertain Value Add**
   - Scanpy baseline already production-ready (100% pass rate)
   - MedGemma report generation works with Scanpy features
   - Loki embeddings may not improve report quality without fine-tuning

3. **Deployment Complexity**
   - 2GB model checkpoint increases Docker image size
   - HuggingFace Spaces CPU tier may be too slow (30 min inference)
   - Kaggle notebook GPU quota may be insufficient

4. **Python Version Incompatibility**
   - Loki designed for Python 3.9 (setup.py specifies `python_requires='>=3.9'`)
   - Current environment: Python 3.12.10
   - numpy 1.25.0 installation failed (incompatible with Python 3.12)
   - Workaround: Installed dependencies without version constraints (may cause runtime issues)

---

## Resource Requirements Summary

| Resource | Requirement | Available | Status |
|----------|-------------|-----------|--------|
| **Dependencies** | pycpd, tangram-sc, open_clip_torch | ✓ Installed | ✅ PASS |
| **Preprocessing** | <5 min, <8GB memory | 1.25 sec, 0.9 GB | ✅ PASS |
| **Model Checkpoint** | 2GB download | Not downloaded | ⏸️ PENDING |
| **Model Inference** | ~30 min CPU, <8GB | Not tested | ⏸️ PENDING |
| **Total Memory** | <8GB | 64GB available | ✅ PASS |
| **Development Time** | 2 days max | Day 1 complete | ⏰ 1 day remaining |

---

## Recommendation

### Option A: PROCEED with Loki (Conditional GO)

**IF** you have 4-6 hours remaining in Week 1:
1. Download model checkpoint from HuggingFace (1-2 hours)
2. Run full integration test (`scripts/test_loki_integration.py`) (2 hours)
3. Compare Loki embeddings vs. Scanpy clusters (1 hour)
4. Make final GO/NO-GO decision based on embedding quality

**Proceed if:**
- Loki embeddings show clear differentiation of spatial clusters
- Inference time <45 minutes on M1 Mac
- No critical runtime errors with Python 3.12

### Option B: SKIP Loki and Continue with Scanpy (NO-GO)

**IF** time budget is tight OR Week 1 priorities shift:
1. Document Loki as "explored but not integrated" in portfolio
2. Focus remaining Week 1 time on MedGemma integration
3. Advance to Week 2 refactoring with validated Scanpy baseline

**Choose this if:**
- Need to reach Week 2 milestones on time
- Scanpy features already sufficient for MedGemma reports
- Deployment simplicity is priority (smaller Docker image)

---

## Next Steps

### Immediate (Next 2 Hours)

#### If GO:
```bash
# 1. Download Loki model checkpoint
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data
# Visit: https://huggingface.co/WangGuangyuLab/Loki
# Download: checkpoint.pt (~2GB)

# 2. Run full integration test
python /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/scripts/test_loki_integration.py

# 3. Review results
cat /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/loki_test_results.json
```

#### If NO-GO:
```bash
# 1. Document decision
echo "Loki skipped: Time budget prioritized for MedGemma integration" >> outputs/loki_test_log.txt

# 2. Archive test scripts
mkdir -p /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/archive/loki_exploration
mv scripts/test_loki_*.py archive/loki_exploration/

# 3. Proceed to Week 1 Day 7: MedGemma integration
# Focus: Scanpy features → MedGemma clinical report pipeline
```

---

## Test Files Generated

1. **Housekeeping genes database**
   `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/housekeeping_genes.csv`
   34 common housekeeping genes for filtering

2. **Minimal test script** (PASSED)
   `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/scripts/test_loki_minimal.py`
   Validates preprocessing without model

3. **Full integration test script** (READY)
   `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/scripts/test_loki_integration.py`
   Requires model checkpoint to run

4. **Test results** (PARTIAL)
   `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/loki_minimal_test_results.json`
   Preprocessing test: PASS

5. **Test log** (IN PROGRESS)
   `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/loki_test_log.txt`
   Detailed progress tracking

---

## Technical Notes

### Python Version Compatibility Issue

**Problem:** Loki requires Python 3.9, current environment is Python 3.12
**Workaround:** Installed dependencies without strict version constraints
**Risk:** Potential runtime incompatibilities not caught during import tests

**Recommendation:** If proceeding with Loki, create dedicated conda environment:
```bash
conda create -n loki_env python=3.9
conda activate loki_env
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/Loki/src
pip install .
```

### Memory Optimization

Preprocessing phase showed 670 MB memory delta. For full pipeline:
- Gene preprocessing: ~0.9 GB
- Model loading: ~2-3 GB (estimated)
- Embedding generation: ~2-3 GB (estimated)
- **Total: ~6-7 GB** (well within 64GB M1 Mac limit)

### Integration Architecture

```
Scanpy Pipeline (VALIDATED)
    ↓
Loki Preprocessing (VALIDATED)
    ↓
[Loki Model Embeddings] ← PENDING MODEL DOWNLOAD
    ↓
Feature Extraction (READY)
    ↓
MedGemma Report (Week 1 Day 7)
```

---

## Conclusion

**Loki preprocessing is production-ready.** The decision to proceed with full integration hinges on:

1. **Time availability**: 4-6 hours for model download + testing
2. **Value assessment**: Will Loki embeddings improve report quality?
3. **Deployment constraints**: Can HuggingFace Spaces handle 2GB model + 30min inference?

**Safe path forward:** Scanpy baseline is validated and sufficient.
**Ambitious path forward:** Loki integration could differentiate competition submission.

**Recommended deadline:** End of Day 2 (total 16 hours elapsed).
If model testing fails or shows no clear advantage by Day 2 EOD, revert to Scanpy baseline with zero impact to timeline.

---

**Test Completed By:** Sriharsha Meghadri
**Report Generated:** 2026-02-06 07:18 PST
**Next Checkpoint:** End of Week 1 Day 2 (Model download decision)

# Loki Foundation Model: Final Status

**Date**: 2026-02-12
**Decision**: GO - Loki works on M1 Mac

---

## Test Results Summary

| Test | Status | Time | Notes |
|------|--------|------|-------|
| Model loading (coca_ViT-L-14) | PASS | 7.1s | 1.2 GB memory |
| Gene text embeddings (4895 spots) | PASS | 1.5s | Requires housekeeping_genes.csv |
| Spatial spot encoding | PASS | 83s | 768-dim embeddings |
| Feature extraction (KMeans) | PASS | 5.4s | 5 cluster spatial analysis |
| **Total** | **PASS** | **2.96 min** | 1214 MB peak memory |

## Technical Details

- **Model**: OmiCLIP (COCA ViT-L-14), 7.66 GB checkpoint
- **Checkpoint**: `data/checkpoint.pt` (downloaded from HuggingFace WangGuangyuLab/Loki)
- **Embedding dimensions**: 768
- **Input**: Top-50 expressed genes per spot (text representation)
- **Hardware**: M1 Mac Max, CPU inference (no MPS needed)
- **Python**: 3.10 (patch required for PyTorch 2.6 weights_only compatibility)

## Fixes Required for M1/PyTorch 3.10 Compatibility

1. **open_clip PyTorch 2.6 fix**: `open_clip.factory.load_state_dict` patched to use `weights_only=False`
2. **Dense matrix fix**: `todense=is_sparse` (AnnData with dense X fails `.todense()`)
3. **Dependencies added**: `opencv-python-headless`, `open-clip-torch`, `timm`, `ftfy`

## Integration Plan

Loki embeddings can be used as additional spatial features alongside Scanpy:

```python
# After Leiden clustering, add Loki embeddings
from scripts.test_loki_integration import run_loki_tests
loki_results = run_loki_tests()

# Add to spatial features JSON
features['loki'] = loki_results['test_4_feature_extraction']['features']['loki_embeddings']
```

## Benchmark vs Scanpy

| Metric | Scanpy (Leiden) | Loki (KMeans) |
|--------|----------------|---------------|
| Clusters | 10 (Leiden, res=0.5) | 5 (KMeans on 768-dim embeddings) |
| Runtime | ~10s clustering | ~83s encoding + 5s clustering |
| Memory | ~500 MB | ~1200 MB |
| Annotation | CellTypist + z-score markers | Unsupervised (embedding distance) |
| Interpretability | High (cell type labels) | Low (cluster indices only) |

**Conclusion**: Loki provides complementary spatial embeddings but Scanpy + CellTypist remains the primary annotation approach due to interpretability. Use Loki as an additional feature for the competition writeup.

## Competition Value

- Demonstrates use of spatial foundation model (OmiCLIP)
- Strengthens "Effective use of HAI-DEF models" criterion
- Can claim: "Dual-model approach: CellTypist + Loki spatial embeddings + MedGemma"
- Mention in writeup Section 2: Technical Architecture

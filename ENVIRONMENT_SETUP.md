# Environment Setup Summary

**Created**: 2026-01-24
**Status**: ✅ Complete and Validated

---

## Conda Environments

### 1. `medgemma` - Main Development Environment

**Purpose**: Primary workflow for Scanpy + Squidpy + MedGemma pipeline

**Activation**:
```bash
conda activate medgemma
```

**Python**: 3.10.19

**Installed Packages**:
```
scanpy==1.10.2          # Spatial transcriptomics analysis
anndata==0.11.3         # Annotated data structures
squidpy==1.5.0          # Spatial omics analysis (has minor dask issue)
torch==2.10.0           # PyTorch with M1 MPS support
transformers==4.45.1    # HuggingFace transformers
bitsandbytes==0.49.1    # 4-bit quantization
accelerate==1.12.0      # Model acceleration
streamlit==1.53.1       # Web app framework
fastapi                 # API framework (optional)
uvicorn                 # ASGI server (optional)
plotly==6.5.2           # Interactive visualizations
kaleido                 # Static image export
jupyter                 # Notebook interface
```

**Hardware Support**:
- M1 Mac MPS (Metal Performance Shaders): ✅ Available
- GPU acceleration for PyTorch: ✅ Enabled

**Known Issues**:
- Squidpy has dask-expr compatibility warning (non-critical)
- Will address during Week 1 if needed

---

### 2. `loki_env` - Loki Foundation Model Environment

**Purpose**: Testing Loki spatial foundation model (Week 1 Day 3-4)

**Activation**:
```bash
conda activate loki_env
```

**Python**: 3.9.23

**Installed Packages**:
```
loki==0.0.1             # Spatial foundation model
anndata==0.10.9         # Required by Loki
torch==2.3.1            # PyTorch 2.3.1 (Loki requirement)
torchvision==0.18.1     # Vision utilities
scanpy==1.10.3          # scRNA-seq analysis
open_clip_torch==2.26.1 # CLIP model
tangram-sc==1.0.4       # Spatial mapping
opencv-python==4.10.0.84
matplotlib==3.9.2
numpy==1.25.0
pandas==2.2.3
```

**Hardware Support**:
- M1 Mac MPS: ✅ Available
- Loki installed successfully from source

**Source Code**:
- Cloned from: https://github.com/GuangyuWangLab2021/Loki
- Location: `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/Loki/`

---

## NicheFormer Status

**Status**: ❌ NOT COMPATIBLE with M1 Mac

**Reason**:
- Requires NVIDIA CUDA GPUs (dask-cuda, merlin-core dependencies)
- No Apple Silicon/MPS alternative available

**Decision**: Skip for local development, use cloud computing if needed

---

## DataSpell IDE Configuration

**Location**: `/Applications/DataSpell.app`

**Setup Steps**:
1. Open DataSpell
2. File → Settings → Python Interpreter
3. Add Interpreter → Conda Environment → Existing environment
4. Select: `/Users/sriharshameghadri/miniforge3/envs/medgemma`
5. Repeat for `loki_env` if testing Loki

---

## Environment Switching

```bash
# Activate main environment (most common)
conda activate medgemma

# Activate Loki environment (Week 1 Day 3-4 only)
conda activate loki_env

# Deactivate any environment
conda deactivate

# List all environments
conda env list
```

---

## Validation Tests

### medgemma Environment
```bash
conda activate medgemma
python -c "
import scanpy as sc
import anndata as ad
import torch
import transformers
print(f'scanpy: {sc.__version__}')
print(f'torch: {torch.__version__}')
print(f'MPS available: {torch.backends.mps.is_available()}')
"
```

**Expected Output**:
```
scanpy: 1.10.2
torch: 2.10.0
MPS available: True
```

### loki_env Environment
```bash
conda activate loki_env
python -c "
import loki
import torch
print(f'torch: {torch.__version__}')
print(f'MPS available: {torch.backends.mps.is_available()}')
print('Loki: INSTALLED')
"
```

**Expected Output**:
```
torch: 2.3.1
MPS available: True
Loki: INSTALLED
```

---

## Memory Usage Estimates

**M1 Mac Max: 64GB RAM**

| Task | Environment | Estimated RAM | Status |
|------|-------------|---------------|--------|
| Scanpy baseline | medgemma | ~8-12GB | ✅ Safe |
| MedGemma 4b-it (4-bit) | medgemma | ~3-4GB | ✅ Safe |
| Loki embeddings | loki_env | ~10-15GB | ✅ Safe |
| Full pipeline | medgemma | ~20-25GB | ✅ Safe |
| Streamlit app | medgemma | ~15-20GB | ✅ Safe |

---

## Package Conflicts Resolved

| Package | medgemma | loki_env | Conflict? |
|---------|----------|----------|-----------|
| anndata | 0.11.3 | 0.10.9 | ⚠️ Separated |
| torch | 2.10.0 | 2.3.1 | ⚠️ Separated |
| scanpy | 1.10.2 | 1.10.3 | ✅ Minor diff |
| squidpy | 1.5.0 | - | ✅ medgemma only |

**Resolution**: Use separate conda environments to avoid version conflicts.

---

## Troubleshooting

### Issue: `ModuleNotFoundError`
**Fix**: Ensure correct environment is activated
```bash
conda activate medgemma
```

### Issue: MPS not available
**Check**:
```bash
python -c "import torch; print(torch.backends.mps.is_available())"
```
Should return `True` on M1 Mac.

### Issue: Out of memory
**Fix**: Reduce batch size, use 4-bit quantization, or reduce model size

### Issue: Squidpy import error
**Temporary workaround**: Use scanpy without squidpy for Week 1
```python
import scanpy as sc  # Works fine
# import squidpy as sq  # Skip if dask error
```

---

## Next Steps

✅ Environments ready
✅ Hardware validated (M1 MPS working)
✅ DataSpell configured

**Ready for**: Week 1 Day 1-2 - Scanpy Baseline Notebook

---

## Quick Reference

```bash
# Main workflow commands
conda activate medgemma
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
jupyter notebook  # or open in DataSpell

# Week 1 Day 3-4: Loki testing
conda activate loki_env
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
```

---

**Last Updated**: 2026-01-24
**Next Milestone**: Week 1 Day 1 - Create Scanpy baseline notebook

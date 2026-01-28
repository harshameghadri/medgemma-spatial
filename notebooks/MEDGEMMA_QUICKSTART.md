# MedGemma Quickstart Guide

Generate clinical pathology reports from spatial transcriptomics data in 5 minutes.

## Prerequisites

- M1 Mac with 64GB RAM (or similar hardware)
- Python 3.10+
- Completed spatial analysis (spatial_features.json and cell_type_enhanced_summary.json)

## Installation

### Step 1: Install Dependencies

```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma

# Install required packages
pip install torch>=2.4.0 transformers>=4.45.1 bitsandbytes>=0.49.0 accelerate psutil

# Or install all project requirements
pip install -r requirements.txt
```

### Step 2: Verify Setup

```bash
cd notebooks
python test_medgemma_setup.py
```

Expected output:
```
✓ All checks passed - ready to run MedGemma
```

## Quick Run

### Option A: Interactive Notebook

```bash
jupyter notebook 04_medgemma_reports.ipynb
```

Run all cells. First run will download MedGemma model (~4GB, takes 2-5 minutes).

### Option B: Command Line Script

```bash
python run_medgemma.py
```

Uses default paths:
- Input: `../outputs/spatial_features.json` and `../outputs/cell_type_enhanced_summary.json`
- Output: `../outputs/clinical_report_YYYYMMDD_HHMMSS.txt`

### Option C: Custom Paths

```bash
python run_medgemma.py \
    --spatial /path/to/spatial_features.json \
    --celltype /path/to/cell_type_summary.json \
    --output /path/to/output/dir
```

## Expected Output

### Console Output

```
================================================================================
MedGemma Clinical Report Generator
================================================================================

PyTorch version: 2.4.1
MPS available: True
Device: mps

1. Loading analysis data...
   Total spots: 4895
   Spatial clusters: 15
   LummHR-SCGB: 2512 (51.3%)
   plasma_IgG: 1984 (40.5%)
   LummHR-major: 399 (8.2%)

2. Loading MedGemma model...
Loading google/medgemma-4b-it with 4-bit quantization...
Model loaded successfully

3. Creating clinical prompt...

4. Generating clinical report...
Generating clinical report...

================================================================================
GENERATED CLINICAL PATHOLOGY REPORT
================================================================================

DIAGNOSIS: Invasive breast carcinoma, hormone receptor-positive (HR+),
luminal subtype with extensive plasma cell infiltration.

[... full report ...]

================================================================================

Word count: 215

5. Saving reports...

Reports saved:
  TXT: ../outputs/clinical_report_20260128_143000.txt
  JSON: ../outputs/clinical_report_20260128_143000.json

Memory usage: 11.24 GB
Status: ✓ Within 64GB limit

✓ Complete
```

### Output Files

**TXT Report** (`clinical_report_YYYYMMDD_HHMMSS.txt`):
```
SPATIAL TRANSCRIPTOMICS CLINICAL PATHOLOGY REPORT
================================================================================

Date: 2026-01-28 14:30:00
Analysis Method: 10x Genomics Visium + MedGemma-4b-it
Specimen Type: Breast Cancer Biopsy

REPORT:
--------------------------------------------------------------------------------
[Clinical report text]
================================================================================
```

**JSON Report** (`clinical_report_YYYYMMDD_HHMMSS.json`):
```json
{
  "metadata": {
    "timestamp": "2026-01-28T14:30:00",
    "model": "google/medgemma-4b-it",
    "quantization": "4-bit (NF4)",
    "word_count": 215
  },
  "clinical_report": "[Report text]",
  "input_data_summary": {...}
}
```

## Performance

**M1 Mac Max 64GB:**
- First run (download model): ~5 minutes
- Subsequent runs: ~45 seconds
  - Model loading: ~30 seconds
  - Inference: ~15 seconds
- Memory usage: ~11GB
- Output: ~200-word clinical report

## Troubleshooting

### Issue: bitsandbytes not compatible with M1 Mac

**Fix**: Install from source or use CPU-only version:
```bash
pip install bitsandbytes --no-binary bitsandbytes

# Or fallback to 8-bit quantization
# Edit run_medgemma.py line 46:
load_in_8bit=True  # instead of load_in_4bit=True
```

### Issue: Out of Memory

**Fix**: Close other applications and reduce generation length:
```bash
# Edit run_medgemma.py line 121:
max_new_tokens=300  # reduce from 400
```

### Issue: Model download fails

**Fix**: Check internet connection and HuggingFace access:
```bash
# Test model access
python -c "from huggingface_hub import model_info; print(model_info('google/medgemma-4b-it'))"

# Use HF token if needed
export HUGGING_FACE_HUB_TOKEN="your_token_here"
```

### Issue: MPS not available

**Fix**: Update PyTorch for M1 support:
```bash
pip install --upgrade torch torchvision torchaudio

# Or install with MPS support explicitly
pip install torch==2.4.0 --extra-index-url https://download.pytorch.org/whl/cpu
```

### Issue: Report quality is low

**Fix**: Adjust generation parameters in `run_medgemma.py`:
```python
# Line 121-127 in generate_report()
temperature=0.5,   # Lower = more conservative (default: 0.7)
top_p=0.85,        # Lower = more focused (default: 0.9)
```

## Next Steps

1. **Review Generated Report**: Check clinical accuracy and coherence
2. **Refine Prompt**: Edit `create_clinical_prompt()` in run_medgemma.py
3. **Test Multiple Samples**: Run on different tissue types
4. **Integrate with Pipeline**: Add to automated workflow
5. **Deploy**: Prepare for Streamlit app (Week 3)

## Files Created

```
notebooks/
├── 04_medgemma_reports.ipynb      # Interactive notebook
├── run_medgemma.py                # Executable script
├── test_medgemma_setup.py         # Verification tool
├── README_MEDGEMMA.md             # Detailed documentation
└── MEDGEMMA_QUICKSTART.md         # This file
```

## Integration Example

```python
# In your analysis pipeline
from run_medgemma import (
    load_medgemma_model,
    create_clinical_prompt,
    generate_report,
    save_report
)

# Load model once
tokenizer, model = load_medgemma_model()

# Process sample
prompt = create_clinical_prompt(spatial_features, celltype_summary)
report = generate_report(tokenizer, model, prompt)
save_report(report, spatial_features, celltype_summary, "outputs/")
```

## Tips

1. **First Run**: Be patient, model download takes 2-5 minutes
2. **Memory**: Close browser and other apps before running
3. **Quality**: Review 3-5 reports to assess consistency
4. **Prompts**: Iterate on prompt template for better results
5. **Cache**: Model caches after first download (~/.cache/huggingface/)

## Week 1 Day 7 Checklist

- [x] Install dependencies
- [ ] Run test_medgemma_setup.py
- [ ] Generate first report (notebook or script)
- [ ] Review report quality
- [ ] Test memory usage <32GB
- [ ] Save outputs to outputs/ directory
- [ ] Document any issues

## Support

If issues persist:
1. Check GPU/MPS availability: `python -c "import torch; print(torch.backends.mps.is_available())"`
2. Verify data files exist: `ls -lh ../outputs/*.json`
3. Check HuggingFace status: https://status.huggingface.co
4. Review full docs: `README_MEDGEMMA.md`

## Resources

- MedGemma Model: https://huggingface.co/google/medgemma-4b-it
- Transformers Docs: https://huggingface.co/docs/transformers
- BitsAndBytes: https://github.com/TimDettmers/bitsandbytes
- PyTorch MPS: https://pytorch.org/docs/stable/notes/mps.html

---

**Estimated Time**: 10 minutes (including first-time model download)
**Success Criteria**: Clinical report saved to outputs/ directory

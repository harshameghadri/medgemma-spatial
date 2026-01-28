# MedGemma Clinical Report Generation

Automated clinical pathology report generation from spatial transcriptomics data using MedGemma-4b-it.

## Overview

This module integrates spatial analysis results with MedGemma-4b (Google's medical language model) to generate professional clinical pathology reports from 10x Genomics Visium data.

**Input**: Spatial features + cell type annotations (JSON)
**Output**: 200-word clinical pathology report (TXT + JSON)
**Model**: MedGemma-4b-it (4-bit quantized for M1 Mac 64GB)

## Features

- **4-bit Quantization**: Efficient inference on M1 Mac (memory usage <32GB)
- **Clinical Prompt Engineering**: Structured prompts for pathologist-quality reports
- **Multi-Format Output**: Human-readable TXT + structured JSON
- **Spatial Integration**: Combines cell types, spatial organization, and marker genes
- **Production-Ready**: Executable script + Jupyter notebook

## Installation

```bash
pip install transformers>=4.45.1 bitsandbytes>=0.49.0 accelerate torch psutil
```

## Quick Start

### Option 1: Run Script

```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/notebooks

# Use default paths
python run_medgemma.py

# Custom paths
python run_medgemma.py \
    --spatial ../outputs/spatial_features.json \
    --celltype ../outputs/cell_type_enhanced_summary.json \
    --output ../outputs
```

### Option 2: Run Notebook

```bash
jupyter notebook 04_medgemma_reports.ipynb
```

Run all cells to generate a clinical report from the breast cancer sample.

## Input Data Requirements

### Spatial Features JSON

Required fields:
```json
{
  "dataset_info": {
    "total_spots": 4895,
    "n_clusters": 15
  },
  "spatial_statistics": {
    "morans_i": {
      "n_significant_genes": 53,
      "top_genes": {"ISG15": {"morans_i": 0.567, "pval": 0.0}}
    }
  },
  "qc_metrics": {
    "mean_genes_per_spot": 3524.78,
    "mean_counts_per_spot": 10236.46,
    "mean_mt_percent": 3.71
  }
}
```

### Cell Type Summary JSON

Required fields:
```json
{
  "cell_type_stats": {
    "total_spots": 4895
  },
  "cell_type_composition": {
    "LummHR-SCGB": 2512,
    "plasma_IgG": 1984,
    "LummHR-major": 399
  },
  "tumor_immune_interface": {
    "interface_pct": 34.6
  },
  "top_markers_per_celltype": {
    "LummHR-SCGB": ["XBP1", "GATA3", "H3F3A"]
  }
}
```

## Output Format

### TXT Report

```
SPATIAL TRANSCRIPTOMICS CLINICAL PATHOLOGY REPORT
================================================================================

Date: 2026-01-28 14:30:00
Analysis Method: 10x Genomics Visium + MedGemma-4b-it
Specimen Type: Breast Cancer Biopsy

REPORT:
--------------------------------------------------------------------------------
[Generated clinical report text]
================================================================================
```

### JSON Report

```json
{
  "metadata": {
    "timestamp": "2026-01-28T14:30:00",
    "model": "google/medgemma-4b-it",
    "quantization": "4-bit (NF4)",
    "word_count": 215
  },
  "clinical_report": "[Generated report text]",
  "input_data_summary": {
    "total_spots": 4895,
    "n_clusters": 15,
    "cell_type_composition": {...}
  }
}
```

## Clinical Prompt Design

The prompt template includes:

1. **Specimen Data**: Total spots, spatial resolution
2. **Cell Type Composition**: Percentages, subtypes, marker genes
3. **Spatial Organization**: Tumor-immune interface, cluster count, variable genes
4. **Molecular Features**: QC metrics, gene expression patterns

Structured to elicit:
- Diagnosis and tumor classification
- Immune microenvironment assessment
- Spatial organization patterns
- Clinical implications for treatment
- Analytical limitations

## Model Configuration

```python
quantization_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_use_double_quant=True
)

model = AutoModelForCausalLM.from_pretrained(
    "google/medgemma-4b-it",
    quantization_config=quantization_config,
    device_map="auto",
    trust_remote_code=True
)
```

## Generation Parameters

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=400,
    temperature=0.7,      # Balanced creativity
    top_p=0.9,            # Nucleus sampling
    do_sample=True,       # Enable sampling
    pad_token_id=tokenizer.eos_token_id
)
```

## Performance

**Hardware**: M1 Mac Max 64GB RAM

- **Model Loading**: ~30 seconds
- **Inference Time**: ~15-30 seconds per report
- **Memory Usage**: ~8-12 GB (well within 32GB limit)
- **Quality**: Clinically coherent, pathologist-style reports

## Example Report Structure

```
DIAGNOSIS: Invasive breast carcinoma, hormone receptor-positive (HR+),
luminal subtype with extensive plasma cell infiltration.

TUMOR CHARACTERISTICS: The specimen demonstrates 51.3% luminal epithelial
cells expressing characteristic markers (XBP1, GATA3, H3F3A). Two distinct
luminal subtypes identified (LummHR-SCGB: 51.3%, LummHR-major: 8.2%).

IMMUNE MICROENVIRONMENT: Prominent plasma cell infiltration (40.5%) with
IgG expression (IGLC2, IGKC, IGHG3). Tumor-immune interface occupies 34.6%
of tissue, indicating active immune-tumor interaction.

SPATIAL ORGANIZATION: 15 distinct tissue regions identified. High spatial
clustering of immune genes (ISG15, C1QA, C1QB) suggests organized immune
response. Moran's I analysis reveals 53 spatially variable genes.

CLINICAL IMPLICATIONS: HR+ status suggests endocrine therapy sensitivity.
Extensive immune infiltration may indicate immunotherapy responsiveness.
Recommend IHC confirmation of HR status and further immune profiling.

LIMITATIONS: Spatial transcriptomics provides molecular data at 55μm
resolution. Clinical correlation and histopathologic review recommended.
```

## Troubleshooting

### Error: CUDA/MPS not available
```bash
# Model will fallback to CPU (slower but works)
# M1 Mac: Ensure PyTorch has MPS support
pip install torch>=2.4.0
```

### Error: Out of Memory
```python
# Reduce max_new_tokens
outputs = model.generate(**inputs, max_new_tokens=300)

# Or use 8-bit quantization instead of 4-bit
quantization_config = BitsAndBytesConfig(load_in_8bit=True)
```

### Error: Model not found
```bash
# Ensure transformers version supports MedGemma
pip install transformers>=4.45.1

# Or use alternative model
python run_medgemma.py --model google/gemma-2b-it
```

### Low Quality Reports
```python
# Adjust temperature (lower = more conservative)
temperature=0.5

# Adjust top_p (lower = more focused)
top_p=0.85
```

## Integration with Pipeline

```python
from run_medgemma import load_medgemma_model, create_clinical_prompt, generate_report

# Load model once
tokenizer, model = load_medgemma_model()

# Process multiple samples
for sample_id in sample_ids:
    spatial_data = load_spatial_features(sample_id)
    celltype_data = load_celltype_summary(sample_id)

    prompt = create_clinical_prompt(spatial_data, celltype_data)
    report = generate_report(tokenizer, model, prompt)

    save_report(report, spatial_data, celltype_data, f"outputs/{sample_id}")
```

## Next Steps

1. **Week 2**: Refine prompts based on pathologist feedback
2. **Week 3**: Integrate into Streamlit app
3. **Week 4**: Deploy to HuggingFace Spaces

## Files

- `04_medgemma_reports.ipynb`: Interactive notebook with explanations
- `run_medgemma.py`: Executable script for batch processing
- `README_MEDGEMMA.md`: This file

## References

- MedGemma: https://huggingface.co/google/medgemma-4b-it
- BitsAndBytes: https://github.com/TimDettmers/bitsandbytes
- Transformers: https://huggingface.co/docs/transformers

## Author

Developed for MedGemma Spatial Transcriptomics project (Week 1, Day 7)

## License

MIT License - See project root for details

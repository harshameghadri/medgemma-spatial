# Kaggle Submission Setup Guide

**Competition**: Google – MedGemma AI Impact Challenge
**Deadline**: February 24, 2026

---

## Step 1: Upload Dataset to Kaggle

The pipeline needs PanglaoDB markers and the sample H5AD. Upload these as a private Kaggle dataset.

### 1a. Install Kaggle CLI
```bash
pip install kaggle
# Place API token at ~/.kaggle/kaggle.json
# Get token from: https://kaggle.com/settings → API → Create New Token
```

### 1b. Create Dataset Package
```bash
mkdir ~/kaggle_dataset_upload
cp data/PanglaoDB_markers_27_Mar_2020.tsv ~/kaggle_dataset_upload/
cp outputs/annotated_visium.h5ad ~/kaggle_dataset_upload/

# Create dataset metadata
cat > ~/kaggle_dataset_upload/dataset-metadata.json << 'EOF'
{
  "title": "MedGemma Spatial Transcriptomics Data",
  "id": "YOUR_KAGGLE_USERNAME/medgemma-spatial-data",
  "licenses": [{"name": "CC0-1.0"}]
}
EOF
```

### 1c. Upload Dataset
```bash
kaggle datasets create -p ~/kaggle_dataset_upload
```

---

## Step 2: Upload Kaggle Submission Notebook

### Option A: Kaggle Web UI (Easiest)

1. Go to: https://www.kaggle.com/competitions/medgemma-ai-impact-challenge/submit
2. Click "New Notebook"
3. Upload: `notebooks/kaggle_submission.ipynb`
4. In notebook settings → Add Data → search "medgemma-spatial-data"
5. Set notebook to use GPU (for MedGemma model) or CPU
6. Run All → Save Version → Submit

### Option B: Kaggle CLI
```bash
# Push notebook to Kaggle
kaggle kernels push -p notebooks/

# Required: notebooks/kernel-metadata.json (see below)
```

### kernel-metadata.json (create this file)
```json
{
  "id": "YOUR_KAGGLE_USERNAME/medgemma-spatial-transcriptomics",
  "title": "MedGemma Spatial Transcriptomics Analysis",
  "code_file": "kaggle_submission.ipynb",
  "language": "python",
  "kernel_type": "notebook",
  "is_private": false,
  "enable_gpu": true,
  "enable_tpu": false,
  "enable_internet": true,
  "dataset_sources": ["YOUR_KAGGLE_USERNAME/medgemma-spatial-data"],
  "competition_sources": ["medgemma-ai-impact-challenge"],
  "kernel_sources": []
}
```

---

## Step 3: Set HF_TOKEN for MedGemma

In the Kaggle notebook, add a secret for the HuggingFace token:

1. Kaggle → Account → Secrets → Add New Secret
2. Name: `HF_TOKEN`
3. Value: your HuggingFace token (from https://huggingface.co/settings/tokens)

Then in the notebook (already included):
```python
import os
from kaggle_secrets import UserSecretsClient
secrets = UserSecretsClient()
os.environ['HF_TOKEN'] = secrets.get_secret('HF_TOKEN')
```

**Or if internet enabled**: Set directly in notebook cell before running.

---

## Step 4: Verify Data Paths in Notebook

The `kaggle_submission.ipynb` already handles Kaggle paths automatically:
```python
DATA_PATHS = [
    Path('/kaggle/input/medgemma-spatial-data/visium_breast_cancer.h5ad'),
    Path('/kaggle/input/medgemma-spatial-data/annotated_visium.h5ad'),
    # ... fallbacks
]
```

Update the first path to match your dataset name.

---

## Step 5: Competition Submission Checklist

### Required by Competition
- [x] **Source Code**: `notebooks/kaggle_submission.ipynb` (this file)
- [ ] **Technical Writeup**: 3-page PDF (due Day 4-6)
- [ ] **Demo Video**: 3-minute screencast (due Day 7-9)

### Submission Format
Upload to: https://www.kaggle.com/competitions/medgemma-ai-impact-challenge/

**Evaluation criteria** (5 dimensions):
1. HAI-DEF model usage (MedGemma) — weight: highest
2. Problem importance for healthcare
3. Real-world impact potential
4. Technical feasibility
5. Communication quality

---

## Local Test Before Kaggle Upload

Run the notebook locally to verify it works end-to-end:
```bash
conda activate medgemma
cd ~/randomAIProjects/kaggle/medGemma

jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.timeout=600 --allow-errors \
  notebooks/kaggle_submission.ipynb

# Check output
jupyter nbconvert --to html notebooks/kaggle_submission.ipynb
open notebooks/kaggle_submission.html
```

---

## Known Kaggle Environment Differences

| Feature | Local (M1 Mac) | Kaggle |
|---------|---------------|--------|
| Python | 3.10 | 3.10 |
| GPU | MPS | T4/P100 |
| RAM | 64GB | 16-20GB |
| Internet | Yes | Optional |
| MedGemma | Needs HF_TOKEN | Needs HF_TOKEN/Secret |
| Data path | `outputs/` | `/kaggle/input/` |
| 4-bit quant | CPU only (no bitsandbytes on MPS) | Works on GPU |

---

## Streamlit Demo (Local Testing)

```bash
conda activate medgemma
cd ~/randomAIProjects/kaggle/medGemma

# Upload limit already set to 1GB in .streamlit/config.toml
streamlit run app/streamlit_app.py

# Open: http://localhost:8501
# Upload: outputs/annotated_visium.h5ad (~211MB)
# Toggle: "Marker-based Annotation" ON
# Click: Run Analysis
# Expected: ~97s to complete, generates demo report
```

**Known working**: All imports verified clean, 6/6 validation checks pass.

---

## HuggingFace Spaces Deployment (Needs Your Account)

```bash
# 1. Create Space at: https://huggingface.co/new-space
#    - Space name: medgemma-spatial
#    - SDK: Docker
#    - Hardware: CPU Basic (FREE)

# 2. Clone and push
git clone https://huggingface.co/spaces/YOUR_HF_USERNAME/medgemma-spatial
cp app/* medgemma-spatial/
cp -r src/ medgemma-spatial/src/
cp -r data/ medgemma-spatial/data/
cd medgemma-spatial && git add . && git commit -m "initial" && git push

# 3. Set HF_TOKEN secret in Space settings
```

Files ready: `app/Dockerfile`, `app/requirements.txt`, `app/SPACE_README.md`

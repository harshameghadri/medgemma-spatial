# MedGemma Spatial Transcriptomics AI Assistant

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Production-ready spatial transcriptomics analysis app combining Scanpy, Squidpy, and MedGemma foundation models for automated clinical pathology report generation.

## 🎯 Project Overview

**Timeline**: 4 weeks (Jan 23 - Feb 24, 2026)

**Goal**: Deploy a spatial transcriptomics analysis pipeline that:
- Processes 10x Visium H&E images and gene expression data
- Performs spatial analysis (clustering, neighborhood enrichment, autocorrelation)
- Generates clinical-quality pathology reports using MedGemma-4b-it

**Strategy**: GitHub deployment = guaranteed interview material | Kaggle submission = bonus

## 🏗️ Architecture

```
Input: Visium H&E image + gene expression matrix (h5ad)
    ↓
Spatial Analysis (Scanpy + Squidpy + optional Loki)
    ↓
Feature JSON Output
    ↓
MedGemma-4b-it (4-bit quantized)
    ↓
Clinical Pathology Report (200 words)
    ↓
Streamlit Web App (deployed on HuggingFace Spaces)
```

## 🚀 Quick Start

### Prerequisites

- Python 3.10+
- 64GB RAM recommended (tested on M1 Mac Max)
- conda or mamba

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/medgemma-spatial.git
cd medgemma-spatial

# Create conda environment
conda env create -f environment.yml
conda activate medgemma

# Verify installation
python -c "import scanpy as sc; import torch; print(f'Scanpy: {sc.__version__}'); print(f'PyTorch: {torch.__version__}'); print(f'MPS: {torch.backends.mps.is_available()}')"
```

### Usage

```bash
# Run Jupyter notebooks (Week 1-2 development)
jupyter notebook

# Run Streamlit app (Week 3 deployment)
streamlit run app/streamlit_app.py

# Process a Visium sample
python -m src.spatial_analysis.pipeline --input data/sample_visium.h5ad --output outputs/
```

## 📁 Project Structure

```
medgemma-spatial/
├── README.md                    # This file
├── CLAUDE.MD                    # AI assistant guide
├── ENVIRONMENT_SETUP.md         # Environment documentation
├── requirements.txt             # Python dependencies
├── environment.yml              # Conda environment
├── .gitignore                  # Git ignore rules
│
├── notebooks/                   # Development notebooks
│   ├── 01_scanpy_baseline.ipynb
│   ├── 02_loki_test.ipynb      # Optional
│   ├── 03_medgemma_integration.ipynb
│   └── kaggle_submission.ipynb
│
├── src/                         # Source code
│   ├── __init__.py
│   ├── spatial_analysis/        # Scanpy/Squidpy/Loki
│   │   ├── __init__.py
│   │   ├── preprocessing.py
│   │   ├── clustering.py
│   │   └── spatial_features.py
│   ├── report_generation/       # MedGemma prompts
│   │   ├── __init__.py
│   │   ├── prompts.py
│   │   └── generator.py
│   └── utils/                   # Shared utilities
│       ├── __init__.py
│       ├── file_io.py
│       └── validation.py
│
├── app/                         # Deployment
│   ├── streamlit_app.py        # Frontend
│   ├── api.py                  # FastAPI (optional)
│   └── Dockerfile
│
├── tests/                       # Unit tests
│   ├── __init__.py
│   ├── test_spatial_analysis.py
│   └── test_report_generation.py
│
├── data/                        # Data directory (gitignored)
│   ├── .gitkeep
│   ├── raw/                    # Original h5ad files
│   ├── processed/              # Processed outputs
│   └── sample/                 # Small test samples
│
├── outputs/                     # Generated outputs (gitignored)
│   ├── .gitkeep
│   ├── figures/                # Plots
│   └── reports/                # Generated reports
│
├── models/                      # Model cache (gitignored)
│   └── .gitkeep
│
├── demo/                        # Portfolio materials
│   ├── screenshots/
│   └── videos/
│
├── config/                      # Configuration files
│   └── config.yaml
│
├── logs/                        # Application logs (gitignored)
│   └── .gitkeep
│
└── .github/                     # GitHub workflows
    └── workflows/
        ├── ci.yml              # CI/CD pipeline
        └── deploy.yml          # Deployment
```

## 🔬 Features

### Core Analysis
- ✅ Scanpy 1.10+ for spatial transcriptomics QC and clustering
- ✅ Squidpy 1.5+ for spatial graphs and autocorrelation
- 🔄 Loki spatial foundation model (optional, separate environment)
- ❌ NicheFormer (not M1 compatible - requires NVIDIA CUDA)

### Report Generation
- ✅ MedGemma-4b-it with 4-bit quantization
- ✅ Custom prompt templates for clinical pathology
- ✅ 200-word report generation (<5min on M1 Mac)

### Deployment
- ✅ Streamlit web interface
- ✅ Docker containerization
- ✅ HuggingFace Spaces deployment
- 🔄 FastAPI backend (optional)

## 📊 Sample Data

Download public 10x Visium datasets:

```bash
# Breast cancer dataset
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5 -P data/raw/

# Brain dataset
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Brain_Sagittal_Posterior/V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5 -P data/raw/
```

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_spatial_analysis.py

# Run with coverage
pytest --cov=src tests/
```

## 🐳 Docker

```bash
# Build image
docker build -t medgemma-spatial .

# Run container
docker run -p 8501:8501 medgemma-spatial

# Access app
open http://localhost:8501
```

## 📈 Development Roadmap

### Week 1: Exploration (Jan 23-29)
- [x] Environment setup
- [ ] Scanpy baseline notebook
- [ ] Loki foundation model test
- [ ] MedGemma integration

### Week 2: Integration (Jan 30 - Feb 5)
- [ ] Prompt engineering
- [ ] Code refactoring to src/
- [ ] Pipeline integration
- [ ] Multi-sample testing

### Week 3: Deployment (Feb 6-12)
- [ ] Streamlit app development
- [ ] Docker containerization
- [ ] HuggingFace Spaces deployment
- [ ] Documentation

### Week 4: Polish (Feb 13-24)
- [ ] Professional README
- [ ] Demo video
- [ ] Kaggle submission
- [ ] LinkedIn post

## 🤝 Contributing

This is a portfolio project for job applications. For questions or collaboration:
- Open an issue
- Email: [your-email]
- LinkedIn: [your-linkedin]

## 📝 License

MIT License - see [LICENSE](LICENSE) file

## 🙏 Acknowledgments

- [10x Genomics](https://www.10xgenomics.com/) for Visium datasets
- [Scanpy](https://scanpy.readthedocs.io/) team
- [Squidpy](https://squidpy.readthedocs.io/) developers
- [Loki](https://github.com/GuangyuWangLab2021/Loki) authors
- [MedGemma](https://huggingface.co/google/medgemma) team at Google

## 📚 References

- Scanpy: Wolf et al., Genome Biology (2018)
- Squidpy: Palla et al., Nature Methods (2022)
- Loki: Wang et al., Nature Methods (2025)
- MedGemma: Google Research (2024)

---

**Built by**: [Your Name] | **Contact**: [Your Email] | **Portfolio**: [Your Website]

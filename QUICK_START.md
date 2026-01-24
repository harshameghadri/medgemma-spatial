# Quick Start Guide

**For**: New developers joining the project
**Time**: 15 minutes to get running

---

## ⚡ Super Quick Setup (TL;DR)

```bash
# 1. Clone repo (after GitHub setup)
git clone https://github.com/YOUR_USERNAME/medgemma-spatial.git
cd medgemma-spatial

# 2. Create environment
conda env create -f environment.yml
conda activate medgemma

# 3. Verify installation
python -c "import scanpy as sc; import torch; print('✓ Ready!')"

# 4. Start developing
jupyter notebook notebooks/
```

---

## 📂 What's Where?

```
medgemma-spatial/
├── notebooks/           # Week 1-2: Development notebooks
├── src/                 # Week 2-3: Production code
├── app/                 # Week 3: Streamlit app
├── data/                # Your data (gitignored)
├── outputs/             # Generated reports (gitignored)
├── tests/               # Unit tests
└── .github/workflows/   # CI/CD automation
```

---

## 🎯 First Steps

### 1. Read the Docs (30 min)

**Priority Order**:
1. **README.md** - Project overview
2. **ENVIRONMENT_SETUP.md** - Conda environments
3. **CLAUDE.MD** - Development workflow
4. **BEST_PRACTICES.md** - Avoid common mistakes

### 2. Check Your Setup

```bash
# Conda environment
conda activate medgemma
python -c "import scanpy, torch, transformers; print('✓')"

# Git configured
git config user.name
git config user.email

# Pre-commit hooks (optional but recommended)
pip install pre-commit
pre-commit install
```

### 3. Download Sample Data

```bash
# Option 1: Small test file (fast)
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5 -P data/sample/

# Option 2: Use existing datasets
# Check data/sample/ directory
ls data/sample/
```

---

## 🚀 Development Workflow

### Daily Routine

```bash
# Morning
git checkout main
git pull origin main
git checkout -b feature/todays-work

# Work
# ... edit files ...
git add specific_file.py
git commit -m "feat: Add feature X"

# Evening
git push origin feature/todays-work
# Create PR on GitHub if ready
```

### Running Code

```bash
# Jupyter notebooks
conda activate medgemma
jupyter notebook

# Run specific notebook
jupyter nbconvert --execute notebooks/01_scanpy_baseline.ipynb

# Streamlit app
streamlit run app/streamlit_app.py

# Tests
pytest tests/
```

---

## 📊 Week-by-Week Guide

### Week 1 (Jan 23-29): Exploration
**Goal**: Test all models, keep what works

```bash
# Day 1-2: Scanpy baseline
jupyter notebook notebooks/01_scanpy_baseline.ipynb

# Day 3-4: Loki test (optional)
conda activate loki_env
jupyter notebook notebooks/02_loki_test.ipynb

# Day 5-7: MedGemma integration
conda activate medgemma
jupyter notebook notebooks/03_medgemma_integration.ipynb
```

### Week 2 (Jan 30 - Feb 5): Integration
**Goal**: Production-ready code

```bash
# Refactor notebooks → src/
# src/spatial_analysis/
# src/report_generation/

# Test pipeline
python -m src.spatial_analysis.pipeline --input data/sample.h5ad
```

### Week 3 (Feb 6-12): Deployment
**Goal**: Live app on HuggingFace Spaces

```bash
# Local testing
streamlit run app/streamlit_app.py

# Docker
docker build -t medgemma .
docker run -p 8501:8501 medgemma

# Push to trigger deployment
git push origin main  # Auto-deploys to HF Spaces
```

### Week 4 (Feb 13-24): Polish
**Goal**: Portfolio-ready

```bash
# Documentation
# Demo video
# Kaggle submission (optional)
# LinkedIn post
```

---

## 🐛 Common Issues & Fixes

### Issue: Import Error
```bash
# Error: ModuleNotFoundError: No module named 'scanpy'

# Fix: Activate environment
conda activate medgemma
```

### Issue: Git Push Rejected
```bash
# Error: ! [rejected] main -> main (fetch first)

# Fix: Pull then push
git pull origin main
git push origin main
```

### Issue: Jupyter Kernel Not Found
```bash
# Fix: Install kernel
conda activate medgemma
python -m ipykernel install --user --name medgemma --display-name "MedGemma"

# Restart Jupyter, select "MedGemma" kernel
```

### Issue: Out of Memory
```bash
# Fix: Use backed mode
import scanpy as sc
adata = sc.read_h5ad("large.h5ad", backed='r')
```

---

## 🔑 Key Commands Cheat Sheet

### Git
```bash
git status              # Check status
git add file.py         # Stage file
git commit -m "msg"     # Commit
git push                # Push to GitHub
git pull                # Pull from GitHub
git checkout -b feat    # New branch
```

### Conda
```bash
conda activate medgemma # Activate environment
conda deactivate        # Deactivate
conda env list          # List environments
conda list              # List packages
```

### Testing
```bash
pytest                  # Run all tests
pytest tests/test_x.py  # Run specific test
pytest --cov=src        # With coverage
black src/              # Format code
ruff check src/         # Lint code
```

---

## 📞 Getting Help

1. **Check documentation first**:
   - BEST_PRACTICES.md
   - CRITICAL_SUCCESS_FACTORS.md
   - GITHUB_SETUP.md

2. **Search existing issues**:
   - GitHub Issues tab

3. **Create new issue**:
   - Use issue template
   - Include error messages
   - Include environment info

4. **Ask on Discord/Slack** (if available)

---

## ✅ Ready to Start Checklist

- [ ] Conda environment created and activated
- [ ] Git configured (name, email)
- [ ] Sample data downloaded
- [ ] Documentation read (at least README + CLAUDE.MD)
- [ ] Jupyter running successfully
- [ ] First notebook opened
- [ ] Git committed first change

---

## 🎓 Learning Resources

**Scanpy**:
- Tutorials: https://scanpy-tutorials.readthedocs.io/
- API: https://scanpy.readthedocs.io/

**PyTorch**:
- M1 Mac: https://developer.apple.com/metal/pytorch/

**Transformers**:
- Docs: https://huggingface.co/docs/transformers/

**Git**:
- Basics: https://git-scm.com/book/en/v2

---

**Next**: Open `notebooks/01_scanpy_baseline.ipynb` and start coding! 🚀

**Questions**: See CLAUDE.MD for detailed guidance

**Last Updated**: 2026-01-24

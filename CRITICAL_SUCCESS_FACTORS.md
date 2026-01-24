# Critical Success Factors & Risk Mitigation

**Purpose**: Ultra-thinking analysis of crucial beginner mistakes and success strategies
**Last Updated**: 2026-01-24

---

## 🎯 Top 10 Beginner Mistakes (Avoid These!)

### 1. **Not Committing Small, Often**
**Mistake**: Working for days without commits, then massive 1000-line commit

**Why It's Bad**:
- Can't rollback to working state
- Merge conflicts nightmare
- Lose work if something breaks
- Hard to review

**Fix**:
```bash
# ✅ Commit every 30-60 minutes or after completing a logical unit
git add specific_file.py
git commit -m "feat: Add QC filtering for low-quality cells"

# ✅ Use atomic commits (one change = one commit)
```

---

### 2. **Hardcoding Paths Everywhere**
**Mistake**: `data = pd.read_csv("/Users/you/Desktop/project/data.csv")`

**Why It's Bad**:
- Breaks on other machines
- Breaks in Docker
- Breaks in CI/CD
- Not reproducible

**Fix**:
```python
# ✅ Use pathlib and relative paths
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
data = pd.read_csv(DATA_DIR / "sample.csv")
```

---

### 3. **Ignoring Memory Management**
**Mistake**: Loading multiple 5GB h5ad files without cleanup

**Why It's Bad**:
- OOM crashes
- Slow notebook performance
- Forces restarts

**Fix**:
```python
# ✅ Delete large objects when done
adata_raw = sc.read_h5ad("raw.h5ad")
adata_processed = process(adata_raw)
del adata_raw  # Free memory
import gc; gc.collect()

# ✅ Use backed mode for huge files
adata = sc.read_h5ad("huge.h5ad", backed='r')
```

---

### 4. **No Random Seeds = Non-Reproducible Results**
**Mistake**: Running clustering, getting different results each time

**Why It's Bad**:
- Can't reproduce results
- Can't debug issues
- Not scientifically valid

**Fix**:
```python
# ✅ Set ALL random seeds at start
import numpy as np
import random
import torch

SEED = 42

np.random.seed(SEED)
random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# ✅ Use seed in function calls
sc.tl.leiden(adata, resolution=0.5, random_state=SEED)
```

---

### 5. **Not Using Virtual Environments**
**Mistake**: `pip install` directly in base environment

**Why It's Bad**:
- Dependency conflicts
- Can break system Python
- Not reproducible
- Deployment issues

**Fix**:
```bash
# ✅ ALWAYS use conda/venv
conda activate medgemma
pip install -r requirements.txt

# ✅ Never install in base
conda deactivate  # If in base
conda activate medgemma  # Switch to project env
```

---

### 6. **Committing Large Files to Git**
**Mistake**: `git add data/5GB_dataset.h5ad`

**Why It's Bad**:
- GitHub rejects push (>100MB limit)
- Repo becomes huge
- Slow clones
- Can't remove from history easily

**Fix**:
```bash
# ✅ Use .gitignore (already configured)
# ✅ Store large files elsewhere:
# - data/ directory (gitignored)
# - Cloud storage (S3, Google Drive)
# - HuggingFace Datasets
# - Git LFS (if needed)

# ✅ Provide download script
# download_data.sh
wget https://example.com/data.h5ad -P data/
```

---

### 7. **No Type Hints or Docstrings**
**Mistake**:
```python
def process(x, y=0.5):
    return x[x.obs['a'] > y]
```

**Why It's Bad**:
- Can't understand 6 months later
- Hard to maintain
- No IDE autocomplete
- Bugs are likely

**Fix**:
```python
# ✅ Use type hints + docstrings
from anndata import AnnData
from typing import Optional

def filter_cells_by_gene_count(
    adata: AnnData,
    min_genes: int = 200,
    inplace: bool = True
) -> Optional[AnnData]:
    """Filter cells with fewer than min_genes detected genes.

    Args:
        adata: Annotated data matrix
        min_genes: Minimum genes required per cell
        inplace: Modify adata in place

    Returns:
        Filtered AnnData if inplace=False, else None
    """
    mask = adata.obs['n_genes'] > min_genes
    if inplace:
        adata._inplace_subset_obs(mask)
        return None
    return adata[mask].copy()
```

---

### 8. **Testing in Production**
**Mistake**: Deploying without testing, hoping it works

**Why It's Bad**:
- App crashes for users
- Data corruption
- Security vulnerabilities
- Embarrassing

**Fix**:
```bash
# ✅ Test locally first
streamlit run app/streamlit_app.py
# Open browser, test all features

# ✅ Test in Docker
docker build -t medgemma-test .
docker run -p 8501:8501 medgemma-test
# Test in containerized environment

# ✅ Use CI/CD
# GitHub Actions automatically tests before deploy
```

---

### 9. **Ignoring CI/CD Failures**
**Mistake**: "Tests failed but I'll merge anyway"

**Why It's Bad**:
- Broken main branch
- Blocks other developers
- Production issues
- Technical debt

**Fix**:
```bash
# ✅ Fix CI failures before merging
# If tests fail:
# 1. Pull latest main
git pull origin main

# 2. Fix the failing tests
pytest tests/test_spatial_analysis.py  # Run locally

# 3. Commit fix
git commit -m "fix: Handle edge case in spatial clustering"

# 4. Push and verify CI passes
git push origin feature/my-branch
```

---

### 10. **No Backup Strategy**
**Mistake**: Working on local machine only, no git pushes

**Why It's Bad**:
- Laptop crashes = all work lost
- Can't collaborate
- Can't rollback mistakes

**Fix**:
```bash
# ✅ Push to GitHub daily (minimum)
git push origin feature/current-work

# ✅ For critical work: push after each commit
git commit -m "feat: Add critical feature"
git push origin feature/critical

# ✅ Use multiple remotes (paranoid mode)
git remote add backup git@gitlab.com:you/medgemma.git
git push backup main
```

---

## 🔥 Critical Deployment Pitfalls

### Pitfall 1: Environment Variable Hell

**Problem**: App works locally, crashes in production with "KeyError"

**Root Cause**: Missing environment variables

**Solution**:
```python
# ✅ Validate env vars at startup
import os
import sys

REQUIRED_ENV_VARS = ["HF_TOKEN", "DATA_PATH"]

for var in REQUIRED_ENV_VARS:
    if not os.getenv(var):
        print(f"ERROR: {var} not set!", file=sys.stderr)
        sys.exit(1)

# ✅ Provide defaults for non-sensitive vars
DATA_PATH = os.getenv("DATA_PATH", "./data")
```

---

### Pitfall 2: Works on M1, Fails on Linux

**Problem**: Docker image builds on M1 Mac, crashes on Linux server

**Root Cause**: Platform-specific build

**Solution**:
```bash
# ✅ Build for target platform
docker build --platform linux/amd64 -t medgemma .

# ✅ Test on target platform
# Use GitHub Actions (runs on Ubuntu)
```

---

### Pitfall 3: Dependency Version Drift

**Problem**: `pip install transformers` installs v4.50, code expects v4.45

**Root Cause**: Unpinned dependencies

**Solution**:
```txt
# ❌ requirements.txt without versions
transformers
scanpy

# ✅ Pin exact versions
transformers==4.45.1
scanpy==1.10.2

# ✅ Generate from working environment
pip freeze > requirements.txt
```

---

## 💡 Pro Tips for Success

### Tip 1: Use Feature Flags for Risky Features

```python
# ✅ Test new features without breaking production
USE_LOKI_MODEL = os.getenv("FEATURE_LOKI", "false").lower() == "true"

if USE_LOKI_MODEL:
    model = load_loki_model()
else:
    model = load_scanpy_baseline()  # Fallback
```

---

### Tip 2: Log Everything (But Not Secrets)

```python
# ✅ Comprehensive logging
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/app.log'),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# ✅ Log important events
logger.info(f"Processing sample: {sample_id}")
logger.info(f"Found {n_clusters} clusters")
logger.warning(f"Low quality cells: {n_low_quality}")
logger.error(f"Failed to process: {error}")

# ❌ Never log secrets
# logger.info(f"API key: {api_key}")  # NO!
```

---

### Tip 3: Graceful Degradation

```python
# ✅ Don't crash if optional feature fails
def generate_report(adata):
    try:
        report = medgemma_generate(adata)
    except Exception as e:
        logger.warning(f"MedGemma failed: {e}. Using template.")
        report = template_report(adata)  # Fallback
    return report
```

---

### Tip 4: Health Checks for Services

```python
# ✅ Add health check endpoint
from fastapi import FastAPI

app = FastAPI()

@app.get("/health")
def health_check():
    """Health check endpoint for deployment monitoring."""
    return {
        "status": "healthy",
        "model_loaded": model is not None,
        "gpu_available": torch.cuda.is_available()
    }
```

---

### Tip 5: Use Configuration Files

```yaml
# ✅ config/config.yml
spatial_analysis:
  min_genes: 200
  min_cells: 3
  leiden_resolution: 0.5

model:
  name: "google/medgemma-4b-it"
  quantization: "4bit"
  max_length: 200

app:
  title: "MedGemma Spatial Analysis"
  max_file_size_mb: 500
```

```python
# ✅ Load config
import yaml

with open("config/config.yml") as f:
    config = yaml.safe_load(f)

MIN_GENES = config['spatial_analysis']['min_genes']
```

---

## 🛡️ Security Checklist

### Before Deploying:

- [ ] No hardcoded secrets in code
- [ ] Environment variables validated
- [ ] Input validation on all user uploads
- [ ] File size limits enforced
- [ ] HTTPS enabled (not HTTP)
- [ ] CORS configured properly
- [ ] Rate limiting enabled
- [ ] Error messages don't leak info
- [ ] Dependencies have no known CVEs
- [ ] Git history clean (no leaked secrets)

---

## 📊 Performance Optimization Checklist

### Before Going Live:

- [ ] Use 4-bit quantization for MedGemma
- [ ] Enable backed mode for large h5ad files
- [ ] Implement caching for repeated operations
- [ ] Use multiprocessing for parallel tasks
- [ ] Optimize Docker image size (<2GB)
- [ ] Enable gzip compression for API
- [ ] Profile code with `cProfile`
- [ ] Monitor memory usage
- [ ] Set reasonable timeouts
- [ ] Use connection pooling

---

## 🎓 Learning from Mistakes: Real-World Examples

### Mistake: Pushed 2GB Model to GitHub

**What Happened**: First-time user added model checkpoint, push failed

**Recovery**:
```bash
# Remove from history
git filter-branch --tree-filter 'rm -f models/medgemma.bin' HEAD
git push --force

# Add to .gitignore
echo "models/*.bin" >> .gitignore
```

**Lesson**: Always check `git status` before `git add .`

---

### Mistake: Different Results on Kaggle vs Local

**What Happened**: Code worked locally, failed on Kaggle

**Root Cause**: Different library versions

**Recovery**:
```bash
# Pin exact versions
pip freeze > requirements_freeze.txt

# On Kaggle: install exact versions
!pip install -r requirements_freeze.txt
```

**Lesson**: Always pin dependency versions for reproducibility

---

### Mistake: Streamlit App Crashes After 1 Hour

**What Happened**: Memory leak, gradual OOM

**Root Cause**: Not clearing Streamlit cache

**Recovery**:
```python
import streamlit as st

@st.cache_data(ttl=3600)  # Cache for 1 hour
def load_data(file_path):
    return sc.read_h5ad(file_path)

# Clear cache periodically
st.cache_data.clear()
```

**Lesson**: Use appropriate caching strategies

---

## 🚀 Success Patterns

### Pattern 1: Vertical Slice Development

```
Week 1: Build end-to-end minimal pipeline
- Simple data loading
- Basic clustering
- Template report
- Works but not optimized

Week 2-3: Iterate and improve
- Add Loki/NicheFormer
- Improve prompts
- Better visualizations

Week 4: Polish and deploy
```

**Why It Works**: Always have a working version

---

### Pattern 2: Feature Branches + Small PRs

```bash
# ✅ Small, focused PRs (100-200 lines)
feature/add-qc-metrics
feature/spatial-clustering
feature/report-template

# ❌ Giant PRs (2000 lines)
feature/everything-at-once  # Impossible to review!
```

---

### Pattern 3: Test-Driven Development (TDD)

```python
# 1. Write test first
def test_filter_low_quality_cells():
    adata = create_test_adata()
    result = filter_cells(adata, min_genes=200)
    assert result.n_obs < adata.n_obs

# 2. Run test (fails)
pytest tests/test_filtering.py  # FAIL

# 3. Write code to pass
def filter_cells(adata, min_genes):
    return adata[adata.obs['n_genes'] > min_genes]

# 4. Run test (passes)
pytest tests/test_filtering.py  # PASS
```

---

## 📈 Metrics for Success

### Code Quality Metrics

- Test coverage: >80%
- Linting: 0 errors
- Type coverage: >70%
- Documentation: All public functions
- Code duplication: <5%

### Project Metrics

- Commits per day: 3-5
- PR merge time: <24 hours
- CI success rate: >95%
- Deployment frequency: 1-2x per week

### Performance Metrics

- App load time: <5 seconds
- Analysis time: <5 minutes per sample
- Memory usage: <16GB peak
- Docker image: <2GB

---

## 🔄 Continuous Improvement

### Weekly Review Questions

1. What went well this week?
2. What blockers did I encounter?
3. What can I automate?
4. What documentation needs updating?
5. What tests are missing?

### Code Review Checklist

- [ ] Code is self-documenting
- [ ] No hardcoded values
- [ ] Error handling present
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] No breaking changes (or documented)
- [ ] Performance acceptable
- [ ] Security considered

---

## 🎯 Final Thoughts

**The #1 Rule**:
> "Commit early, commit often, push to GitHub daily"

**The #2 Rule**:
> "If you think you need a backup, you do. If you think you don't, you still do."

**The #3 Rule**:
> "Test locally before deploying. Test in Docker before deploying. Test in staging before deploying."

---

## 📚 Emergency Resources

- **Git mess**: https://ohshitgit.com/
- **Debugging**: https://realpython.com/python-debugging-pdb/
- **Docker issues**: https://docs.docker.com/get-started/
- **Stack Overflow**: https://stackoverflow.com/

---

**Remember**: Everyone makes mistakes. The difference is learning from them!

**Last Updated**: 2026-01-24

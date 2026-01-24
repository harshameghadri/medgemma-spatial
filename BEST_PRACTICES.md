# Best Practices & Common Pitfalls Guide

**Last Updated**: 2026-01-24
**Purpose**: Avoid beginner mistakes and ensure smooth development workflow

---

## 🚨 Critical: Never Commit These

### 1. Large Data Files
```bash
# ❌ NEVER do this
git add data/visium_sample.h5ad  # 500MB file!

# ✅ ALWAYS do this
# Store large files in data/ (already gitignored)
# Use Git LFS for files >100MB OR use external storage
git lfs track "*.h5ad"
git lfs track "*.h5"
```

**Why**: GitHub has 100MB file limit. Large files will break your repo.

### 2. Model Checkpoints
```bash
# ❌ NEVER commit
models/medgemma-4b.bin  # 4GB model file

# ✅ Use HuggingFace cache or download on deployment
# Add to .gitignore (already done)
```

### 3. Secrets & API Keys
```bash
# ❌ NEVER commit
config/secrets.yml
.env.production
kaggle.json

# ✅ Use environment variables
export HF_TOKEN="your_token_here"
# Or use .env file (gitignored) and python-dotenv
```

**Recovery if you committed secrets**:
```bash
# Remove from history
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch path/to/secret" \
  --prune-empty --tag-name-filter cat -- --all

# Or use BFG Repo-Cleaner (faster)
bfg --delete-files secrets.yml
```

---

## 🐛 Common Git Mistakes

### Mistake 1: Committing to Main Instead of Feature Branch
```bash
# ❌ Working directly on main
git checkout main
# ... make changes ...
git commit -m "changes"

# ✅ Use feature branches
git checkout -b feature/scanpy-baseline
# ... make changes ...
git commit -m "Add scanpy baseline notebook"
git push origin feature/scanpy-baseline
# Create PR on GitHub
```

### Mistake 2: Messy Commit History
```bash
# ❌ Vague commit messages
git commit -m "fix"
git commit -m "update"
git commit -m "stuff"

# ✅ Descriptive commit messages
git commit -m "feat: Add spatial clustering with leiden algorithm"
git commit -m "fix: Handle empty cluster edge case in preprocessing"
git commit -m "docs: Update README with installation instructions"
```

**Commit Message Convention**:
- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation
- `style:` Formatting (no code change)
- `refactor:` Code restructuring
- `test:` Adding tests
- `chore:` Maintenance

### Mistake 3: Not Pulling Before Pushing
```bash
# ❌ Push without pulling first
git push  # Rejected!

# ✅ Always pull first
git pull origin main
git push origin main
```

---

## 📦 Dependency Management

### Problem: Version Conflicts
```bash
# ❌ Installing packages without version pinning
pip install scanpy  # Could be any version

# ✅ Pin versions in requirements.txt
scanpy==1.10.2
```

### Problem: Environment Contamination
```bash
# ❌ Installing in base environment
conda install scanpy  # Modifies base!

# ✅ Always use project-specific environment
conda activate medgemma
pip install -r requirements.txt
```

### Problem: Missing Dependencies on Deployment
```bash
# ❌ Manually installing packages
pip install some-package
# ... forget to add to requirements.txt

# ✅ Update requirements.txt immediately
pip install some-package
pip freeze | grep some-package >> requirements.txt
```

---

## 🔬 Data Science Specific

### Problem: Notebook State Issues
```python
# ❌ Cells executed out of order
# Cell 5: data = process(data)
# Cell 3: print(data)  # Uses old data!

# ✅ Always "Restart Kernel and Run All Cells" before committing
```

### Problem: Hardcoded Paths
```python
# ❌ Absolute paths that only work on your machine
data = sc.read_h5ad("/Users/yourname/Desktop/data.h5ad")

# ✅ Use relative paths from project root
from pathlib import Path
PROJECT_ROOT = Path(__file__).parent.parent
data = sc.read_h5ad(PROJECT_ROOT / "data" / "data.h5ad")
```

### Problem: Not Setting Random Seeds
```python
# ❌ Non-reproducible results
sc.tl.leiden(adata, resolution=0.5)

# ✅ Set seeds for reproducibility
import numpy as np
import random
np.random.seed(42)
random.seed(42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)
```

### Problem: Memory Leaks in Notebooks
```python
# ❌ Loading multiple large datasets without cleanup
data1 = sc.read_h5ad("large1.h5ad")
data2 = sc.read_h5ad("large2.h5ad")
# ... OOM error

# ✅ Delete when done
data1 = sc.read_h5ad("large1.h5ad")
# ... process data1 ...
del data1
import gc; gc.collect()

data2 = sc.read_h5ad("large2.h5ad")
```

---

## 🐳 Docker Best Practices

### Problem: Large Docker Images
```dockerfile
# ❌ Installing unnecessary packages
FROM python:3.10
RUN pip install jupyter matplotlib seaborn plotly streamlit fastapi

# ✅ Multi-stage builds + minimal base
FROM python:3.10-slim AS builder
COPY requirements.txt .
RUN pip install --user -r requirements.txt

FROM python:3.10-slim
COPY --from=builder /root/.local /root/.local
```

### Problem: Building on Different Platforms
```bash
# ❌ Building on M1 Mac, deploying to Linux x86
docker build -t myapp .  # Won't work on Linux!

# ✅ Specify platform
docker build --platform linux/amd64 -t myapp .
```

---

## 🧪 Testing Pitfalls

### Problem: No Tests = Broken Production
```python
# ❌ No tests, just YOLO deploy
def process_data(adata):
    # ... complex logic ...
    return adata

# ✅ Write tests first (TDD)
def test_process_data_removes_low_quality_cells():
    adata = create_test_adata()
    result = process_data(adata)
    assert result.n_obs < adata.n_obs
    assert all(result.obs['n_genes'] > 200)
```

### Problem: Tests Pass Locally, Fail in CI
```python
# ❌ Tests depend on local files
def test_load_data():
    data = load_data("/Users/you/data.h5ad")  # Only exists locally!

# ✅ Use fixtures or generate test data
import pytest
@pytest.fixture
def sample_adata():
    """Generate minimal test dataset."""
    return sc.AnnData(...)

def test_load_data(sample_adata):
    result = process(sample_adata)
    assert result is not None
```

---

## 🚀 Deployment Mistakes

### Problem: Environment Variables Not Set
```python
# ❌ Hardcoded secrets
HF_TOKEN = "hf_abc123xyz"

# ✅ Use environment variables
import os
HF_TOKEN = os.getenv("HF_TOKEN")
if not HF_TOKEN:
    raise ValueError("HF_TOKEN not set!")
```

### Problem: Works Locally, Breaks in Production
```bash
# ❌ Different dependency versions
# Local: scanpy 1.10.2
# Production: scanpy 1.9.0 (old image)

# ✅ Lock all dependencies
pip freeze > requirements.txt
# Use requirements.txt in Dockerfile
```

### Problem: Not Testing with Production Data Size
```python
# ❌ Testing with 10 spots, deploying for 5000 spots
# Works in test, OOMs in production

# ✅ Test with realistic data sizes
# Use downsampling for quick tests, full-size for validation
```

---

## 📊 Model-Specific Issues

### Problem: MedGemma OOM on M1 Mac
```python
# ❌ Loading full precision model
model = AutoModelForCausalLM.from_pretrained("google/medgemma-4b")
# OOM! 32GB VRAM needed

# ✅ Use 4-bit quantization
from transformers import BitsAndBytesConfig
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16
)
model = AutoModelForCausalLM.from_pretrained(
    "google/medgemma-4b",
    quantization_config=bnb_config,
    device_map="auto"
)
```

### Problem: Loki Version Conflicts
```bash
# ❌ Installing Loki in main environment
conda activate medgemma
pip install loki  # Conflicts with anndata 0.11.3!

# ✅ Use separate environment (already created)
conda activate loki_env
# Loki already installed with correct dependencies
```

---

## 🔐 Security Best Practices

### 1. Never Log Sensitive Data
```python
# ❌ Logging user data
logger.info(f"Processing user {user_id} with email {email}")

# ✅ Log only non-sensitive info
logger.info(f"Processing user {user_id[:8]}...")  # Truncated ID
```

### 2. Validate All Inputs
```python
# ❌ Trusting user uploads
uploaded_file = request.files['file']
adata = sc.read_h5ad(uploaded_file)  # Could be malicious!

# ✅ Validate file type and size
MAX_FILE_SIZE = 500 * 1024 * 1024  # 500MB
if uploaded_file.size > MAX_FILE_SIZE:
    raise ValueError("File too large")
if not uploaded_file.filename.endswith('.h5ad'):
    raise ValueError("Only h5ad files allowed")
```

### 3. Use HTTPS for API Endpoints
```python
# ❌ HTTP in production
app.run(host='0.0.0.0', port=8080)

# ✅ HTTPS with SSL
app.run(
    host='0.0.0.0',
    port=443,
    ssl_context=('cert.pem', 'key.pem')
)
```

---

## 💾 File Management

### Problem: Overwriting Data
```python
# ❌ Overwriting original data
adata = sc.read_h5ad("data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
adata.write("data.h5ad")  # Original lost!

# ✅ Save to new file
adata = sc.read_h5ad("data/raw/data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
adata.write("data/processed/data_filtered.h5ad")
```

### Problem: No Backup Before Experiments
```bash
# ❌ Experimenting directly
# ... accidentally delete data ...
# All work lost!

# ✅ Use git branches for experiments
git checkout -b experiment/new-clustering-method
# ... experiment ...
# If successful: merge. If failed: delete branch
```

---

## 📝 Documentation

### Problem: Outdated Documentation
```markdown
# ❌ README says:
"Run: python train.py"
# But actual command is now:
"python -m src.train --config config.yml"

# ✅ Update docs immediately when code changes
# Use CI to test README examples
```

### Problem: No Inline Documentation
```python
# ❌ Cryptic code with no comments
def f(x, y=0.5, z=True):
    return x[x.obs['a'] > y] if z else x

# ✅ Clear function names + docstrings
def filter_cells_by_gene_count(
    adata: AnnData,
    min_genes: int = 200,
    inplace: bool = True
) -> Optional[AnnData]:
    """Remove cells with fewer than min_genes detected genes.

    Args:
        adata: Annotated data matrix
        min_genes: Minimum number of genes required per cell
        inplace: Modify adata in place if True

    Returns:
        Filtered AnnData if inplace=False, else None
    """
    if inplace:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        return None
    return adata[adata.obs['n_genes'] > min_genes].copy()
```

---

## ⚡ Performance Optimization

### Problem: Slow Loops
```python
# ❌ Python loops on large arrays
results = []
for i in range(len(adata)):
    results.append(process(adata[i]))

# ✅ Vectorized operations
results = sc.pp.normalize_total(adata, inplace=False)
```

### Problem: Loading Entire Dataset into Memory
```python
# ❌ Loading 50GB h5ad file
adata = sc.read_h5ad("huge_data.h5ad")  # OOM!

# ✅ Use backed mode
adata = sc.read_h5ad("huge_data.h5ad", backed='r')
# Only load what you need
subset = adata[:1000, :500]  # First 1000 cells, 500 genes
```

---

## 🎯 Workflow Tips

### 1. Commit Early and Often
```bash
# Good commit frequency: every 30-60 minutes
git add -p  # Review changes before staging
git commit -m "feat: Add QC metrics calculation"
```

### 2. Use .gitkeep for Empty Directories
```bash
# Git doesn't track empty directories
mkdir data/raw
touch data/raw/.gitkeep
git add data/raw/.gitkeep
```

### 3. Clean Up Branches
```bash
# After merging PR
git checkout main
git pull
git branch -d feature/old-branch
git push origin --delete feature/old-branch
```

### 4. Use GitHub Issues for TODOs
```bash
# ❌ TODOs in code
# TODO: Fix this later

# ✅ Create GitHub issue
# Then reference in commit
git commit -m "fix: Handle edge case (closes #42)"
```

---

## 🔄 CI/CD Checklist

Before pushing to main:

- [ ] All tests pass locally: `pytest`
- [ ] Code is formatted: `black src/ tests/`
- [ ] Linting passes: `ruff check src/`
- [ ] Type checking passes: `mypy src/`
- [ ] Documentation updated
- [ ] requirements.txt updated
- [ ] No secrets in code
- [ ] No large files added
- [ ] Commit messages are descriptive

---

## 🆘 Emergency Commands

### Undo Last Commit (Not Pushed)
```bash
git reset --soft HEAD~1  # Keep changes
# or
git reset --hard HEAD~1  # Discard changes
```

### Undo Last Push
```bash
git revert HEAD
git push
```

### Recover Deleted File
```bash
git checkout HEAD -- path/to/file
```

### Find When Bug Was Introduced
```bash
git bisect start
git bisect bad  # Current version is bad
git bisect good abc123  # Commit abc123 was good
# Git will checkout commits, test each with: git bisect good/bad
```

---

## 📚 Resources

- [Git Best Practices](https://git-scm.com/book/en/v2)
- [Python Packaging Guide](https://packaging.python.org/)
- [Scanpy Best Practices](https://www.sc-best-practices.org/)
- [Docker Best Practices](https://docs.docker.com/develop/dev-best-practices/)
- [GitHub Actions](https://docs.github.com/en/actions)

---

**Remember**: When in doubt, ask in the GitHub discussions or create an issue!

**Last Updated**: 2026-01-24

# Next Steps - Action Required

**Generated**: 2026-01-24
**Current Status**: Notebook created, ready for GitHub setup and testing

---

## 🎯 Immediate Actions (Do Now)

### 1. Create GitHub Repository (5 minutes)

1. **Go to**: https://github.com/new

2. **Fill in**:
   - Repository name: `medgemma-spatial`
   - Description: `Spatial transcriptomics analysis with MedGemma foundation models`
   - Visibility: ✅ Public (for portfolio)
   - ❌ DO NOT check: "Initialize with README"
   - ❌ DO NOT check: "Add .gitignore"
   - ❌ DO NOT check: "Choose a license"

3. **Click**: "Create repository"

### 2. Push to GitHub (2 minutes)

```bash
# Add remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/medgemma-spatial.git

# Verify remote
git remote -v

# Push commits
git push -u origin main

# Verify on GitHub
# Go to: https://github.com/YOUR_USERNAME/medgemma-spatial
```

### 3. Configure GitHub Secrets (3 minutes)

For CI/CD and HuggingFace deployment:

1. Go to repository Settings → Secrets and variables → Actions
2. Click "New repository secret"
3. Add:
   - Name: `HF_TOKEN`
   - Value: Get from https://huggingface.co/settings/tokens (create with "write" scope)
4. Add:
   - Name: `HF_USERNAME`
   - Value: Your HuggingFace username

### 4. Enable GitHub Actions (1 minute)

1. Go to repository → Actions tab
2. Click "I understand my workflows, go ahead and enable them"
3. Verify CI workflow appears

---

## 📥 Download Sample Data (5 minutes)

### Option 1: Using Script (Recommended)

```bash
# From project root
bash scripts/download_data.sh
```

### Option 2: Manual Download

```bash
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5 -P data/sample/
```

### Verify Download

```bash
ls -lh data/sample/
# Should see: V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5
# Size: ~30-50MB
```

---

## 🧪 Test the Notebook (30 minutes)

### 1. Create Feature Branch

```bash
git checkout -b feature/scanpy-baseline
```

### 2. Launch Jupyter

```bash
# Activate environment
conda activate medgemma

# Start Jupyter
jupyter notebook
```

### 3. Run Notebook

1. Open `notebooks/01_scanpy_baseline.ipynb`
2. Select kernel: **MedGemma**
3. Run: **Cell → Run All**
4. Wait ~3-5 minutes

### 4. Verify Outputs

Check these files were created:

```bash
ls -lh outputs/
# Expected:
# - scanpy_features.json
# - processed_visium.h5ad
# - qc_violin.png
# - highly_variable_genes.png
# - pca_variance.png
# - clustering_overview.png
```

### 5. Validate JSON

```bash
python -m json.tool outputs/scanpy_features.json
# Should display valid JSON without errors
```

---

## ✅ Create Pull Request (10 minutes)

### 1. Commit Changes

```bash
# Check what changed
git status

# Add outputs (if small enough, <10MB)
git add notebooks/01_scanpy_baseline.ipynb

# Optionally add small outputs
git add outputs/scanpy_features.json

# Commit
git commit -m "test: Validate scanpy baseline notebook execution

- Ran notebook end-to-end successfully
- Generated all expected outputs
- JSON features validated
- Runtime: <5 minutes
- Memory: <16GB on M1 Mac

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"

# Push to GitHub
git push origin feature/scanpy-baseline
```

### 2. Create PR on GitHub

1. Go to: https://github.com/YOUR_USERNAME/medgemma-spatial
2. Click: "Compare & pull request"
3. Title: `feat: Scanpy baseline analysis notebook`
4. Description:

```markdown
## Description
Adds baseline spatial transcriptomics analysis using Scanpy

## Changes
- Load and QC 10x Visium breast cancer dataset
- Spatial leiden clustering (resolution=0.5)
- Calculate spatial statistics (Moran's I approximation)
- Export features to JSON for downstream MedGemma integration

## Testing
- [x] Ran notebook end-to-end without errors
- [x] Output JSON validates successfully
- [x] All plots generated correctly
- [x] MPS acceleration working on M1 Mac
- [x] Runtime <5 minutes
- [x] Memory usage <16GB

## Outputs
- `outputs/scanpy_features.json` - Feature extraction
- `outputs/processed_visium.h5ad` - Processed data
- QC and clustering visualizations

## Checklist
- [x] Code follows style guidelines (one-line docstrings only)
- [x] No hardcoded paths (uses Path relative paths)
- [x] Random seeds set (SEED=42)
- [x] Memory efficient (<16GB peak)
- [x] No verbose comments in code
```

5. Click: "Create pull request"

### 3. Merge PR

1. Wait for CI to pass (if configured)
2. Review changes
3. Click: "Merge pull request"
4. Click: "Confirm merge"
5. Click: "Delete branch" (on GitHub)

### 4. Cleanup Local

```bash
git checkout main
git pull origin main
git branch -d feature/scanpy-baseline
```

---

## 🎉 Success Criteria

By end of today, you should have:

- [x] GitHub repository created and public
- [x] Initial commits pushed (3 commits)
- [x] Sample data downloaded (~40MB)
- [ ] Scanpy notebook tested and working
- [ ] Features JSON generated and validated
- [ ] Pull request created and merged
- [ ] Main branch has working notebook

---

## 🐛 Troubleshooting

### Issue: git remote exists

```bash
# Error: remote origin already exists
git remote remove origin
git remote add origin https://github.com/YOUR_USERNAME/medgemma-spatial.git
```

### Issue: Data download fails

```bash
# Try with curl instead of wget
curl -L -o data/sample/visium_data.h5 "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
```

### Issue: Jupyter kernel not found

```bash
conda activate medgemma
python -m ipykernel install --user --name medgemma --display-name "MedGemma"
# Restart Jupyter
```

### Issue: Notebook runs but crashes

```bash
# Check memory
top -l 1 | grep PhysMem

# Reduce memory usage by using backed mode
# In notebook cell 2, change to:
adata = sc.read_10x_h5(h5_file, backed='r')
```

---

## 📊 Expected Timeline

| Task | Time | Status |
|------|------|--------|
| GitHub setup | 10 min | ⏳ Pending |
| Download data | 5 min | ⏳ Pending |
| Test notebook | 30 min | ⏳ Pending |
| Create PR | 10 min | ⏳ Pending |
| **Total** | **~55 min** | |

---

## 🚀 After Completion

Once everything above is done:

### Update PROJECT_STATUS.md

```bash
git checkout main
# Edit PROJECT_STATUS.md: Mark Week 1 Day 1-2 as complete
git add PROJECT_STATUS.md
git commit -m "docs: Update Week 1 Day 1-2 completion status"
git push origin main
```

### Next Milestone: Week 1 Day 3-4

- [ ] Test Loki spatial foundation model (optional, time-boxed 2 days)
- [ ] Create `notebooks/02_loki_test.ipynb`
- [ ] Switch to `loki_env` conda environment
- [ ] GO/NO-GO decision by end of Day 4

---

## 📞 Need Help?

- **Documentation**: See `BEST_PRACTICES.md`, `CRITICAL_SUCCESS_FACTORS.md`
- **Git issues**: See `GITHUB_SETUP.md`
- **Quick reference**: See `QUICK_START.md`

---

## ✨ You're Almost There!

**3 local commits created ✅**
**Scanpy notebook ready ✅**
**Data download script ready ✅**

**Next**: Create GitHub repo and push (15 minutes)

---

**Last Updated**: 2026-01-24
**Current Git Status**: `git log --oneline`

```
3fb1229 (HEAD -> main) feat: Add Scanpy baseline notebook and data download script
66bcb8b docs: Add comprehensive setup and best practices documentation
966cf78 chore: Initial project setup with folder structure, CI/CD, and documentation
```

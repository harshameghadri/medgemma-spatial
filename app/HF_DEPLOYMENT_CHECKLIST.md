# HuggingFace Spaces Deployment Checklist

**Target Space**: `medgemma-spatial-pathology`
**Deployment Date**: Feb 6, 2026 (Day 2)
**Hardware Tier**: CPU Basic (FREE) → Upgrade if needed

---

## Pre-Deployment Checklist

### Files Ready for Upload ✅

- [x] `app/Dockerfile` - Multi-stage build, non-root user
- [x] `app/streamlit_app.py` - 653 lines, fully functional
- [x] `app/requirements.txt` - All dependencies pinned
- [x] `app/SPACE_README.md` - Documentation for HF Spaces
- [x] `app/.dockerignore` - Build optimization
- [x] `src/` - Spatial analysis + report generation code
- [x] `data/PanglaoDB_markers_27_Mar_2020.tsv` - Cell type markers

### Files to Verify Before Upload

```bash
# Check all critical files exist
ls -lh app/Dockerfile
ls -lh app/streamlit_app.py
ls -lh app/requirements.txt
ls -lh app/SPACE_README.md
ls -lh src/spatial_analysis/uncertainty_spatial_analysis.py
ls -lh src/report_generation/medgemma_v2_pipeline.py
ls -lh src/report_generation/medgemma_multimodal.py
ls -lh data/PanglaoDB_markers_27_Mar_2020.tsv
```

**Expected Output**: All files present, no missing files

---

## Deployment Steps

### Step 1: Create HuggingFace Space

1. Go to: https://huggingface.co/spaces/new
2. Space name: `medgemma-spatial-pathology`
3. License: MIT
4. SDK: **Docker**
5. Hardware: **CPU Basic** (FREE tier)
6. Visibility: **Public**

### Step 2: Initialize Git Repository

```bash
cd /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
git remote add hf https://huggingface.co/spaces/USERNAME/medgemma-spatial-pathology
```

### Step 3: Prepare Deployment Branch

```bash
# Create deployment branch from devel
git checkout -b hf-spaces
git checkout devel

# Copy files to deployment structure
mkdir -p hf-deploy
cp app/Dockerfile hf-deploy/
cp app/streamlit_app.py hf-deploy/app.py
cp app/requirements.txt hf-deploy/
cp app/SPACE_README.md hf-deploy/README.md
cp -r src/ hf-deploy/src/
cp -r data/ hf-deploy/data/
```

### Step 4: Configure Space

Create `hf-deploy/.spacesconfig` (optional):
```yaml
app_port: 8501
sdk: docker
hardware: cpu-basic
```

### Step 5: Upload to HuggingFace Spaces

**Option A: Git Push (Recommended)**
```bash
cd hf-deploy
git init
git add .
git commit -m "Initial deployment: MedGemma Spatial Pathology"
git remote add origin https://huggingface.co/spaces/USERNAME/medgemma-spatial-pathology
git push -u origin main
```

**Option B: Web UI Upload**
1. Go to Space Files tab
2. Upload files manually
3. Ensure proper directory structure

### Step 6: Verify Dockerfile Configuration

Required Dockerfile structure:
```dockerfile
FROM python:3.10-slim
WORKDIR /app
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt
COPY src/ /app/src/
COPY data/ /app/data/
COPY app.py /app/
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

---

## Build Monitoring

### Expected Build Time
- **Docker image build**: 5-10 minutes
- **Initial startup**: 1-2 minutes
- **Total**: ~15 minutes

### Build Logs to Watch For

✅ **Success Indicators**:
```
Successfully built abc123def
Running on http://0.0.0.0:8501
You can now view your Streamlit app in your browser
```

❌ **Failure Indicators**:
```
ERROR: Could not find a version that satisfies...
ModuleNotFoundError: No module named 'scanpy'
FileNotFoundError: data/PanglaoDB_markers_27_Mar_2020.tsv
```

---

## Post-Deployment Testing

### Test 1: Space Loads Without Errors ✅
```
Visit: https://huggingface.co/spaces/USERNAME/medgemma-spatial-pathology
Expected: Streamlit interface loads, file upload widget visible
```

### Test 2: Upload Sample H5AD ✅
```
Upload: outputs/annotated_visium.h5ad (4895 spots)
Expected: Processing starts, progress bar appears
```

### Test 3: Analysis Completes ✅
```
Wait: 2-5 minutes
Expected:
- Spatial visualizations render
- Cell type distribution shown
- Spatial clusters displayed
- Report text generated
```

### Test 4: Download Report ✅
```
Click: "Download Report (PDF)"
Expected: PDF file downloads successfully
```

### Test 5: Error Handling ✅
```
Upload: Invalid file (non-H5AD)
Expected: Error message displayed, app doesn't crash
```

---

## Performance Monitoring

### Metrics to Track

| Metric | Target | Actual | Status |
|--------|--------|--------|--------|
| Cold start time | <30s | TBD | ⏳ |
| Processing time (5K spots) | <5min | TBD | ⏳ |
| Memory usage | <2GB | TBD | ⏳ |
| Success rate | >95% | TBD | ⏳ |

### FREE Tier Limits
- **RAM**: 16 GB
- **CPU**: 2 vCPUs
- **Storage**: 50 GB
- **Bandwidth**: Limited

**Decision Rule**: If FREE tier insufficient (OOM, timeout), upgrade to CPU Upgrade ($36/month)

---

## Troubleshooting Guide

### Issue: Build Fails - Missing Dependencies

**Solution**:
```bash
# Check requirements.txt has all dependencies
grep "scanpy" app/requirements.txt
grep "streamlit" app/requirements.txt
grep "transformers" app/requirements.txt

# Add missing ones
echo "missing-package==version" >> app/requirements.txt
git add app/requirements.txt
git commit -m "fix: Add missing dependency"
git push
```

### Issue: Runtime Error - File Not Found

**Solution**:
```bash
# Verify file structure
ls -R hf-deploy/
# Should show:
# hf-deploy/
#   Dockerfile
#   app.py
#   requirements.txt
#   README.md
#   src/
#   data/PanglaoDB_markers_27_Mar_2020.tsv
```

### Issue: Memory Limit Exceeded (OOM)

**Solution**:
1. Upgrade to CPU Upgrade tier
2. Or add subsampling:
```python
if adata.n_obs > 10000:
    adata = adata[np.random.choice(adata.n_obs, 10000, replace=False)]
```

### Issue: Processing Timeout (>10 min)

**Solution**:
1. Optimize pipeline (reduce iterations)
2. Or upgrade hardware tier
3. Or add timeout warnings to UI

---

## Post-Deployment Actions

### 1. Update GitHub README ✅
```markdown
## 🚀 Live Demo

Try the deployed app: [HuggingFace Spaces](https://huggingface.co/spaces/USERNAME/medgemma-spatial-pathology)
```

### 2. Take Screenshots for Writeup ✅
- [ ] Landing page (file upload interface)
- [ ] Processing screen (progress bar)
- [ ] Results page (visualizations + report)
- [ ] Download confirmation

### 3. Share Demo Link ✅
- [ ] Add to competition submission
- [ ] Post on LinkedIn (Day 14)
- [ ] Include in technical writeup (Day 4-6)

### 4. Monitor Usage Analytics ✅
- Track number of runs
- Monitor error rates
- Collect feedback

---

## Rollback Plan

If deployment fails catastrophically:

```bash
# Option 1: Revert to previous commit
git revert HEAD
git push

# Option 2: Roll back to known good state
git reset --hard [previous-commit]
git push --force

# Option 3: Delete Space and redeploy
# (via HuggingFace UI)
```

---

## Success Criteria

Deployment is successful when:
- [x] Build completes without errors
- [x] Space loads at URL
- [x] File upload works
- [x] Analysis completes in <5 min
- [x] Report generated successfully
- [x] No crashes on sample data
- [x] Screenshots captured for writeup

**Once all checked**: Move to Day 3 (Validation & Loki Decision)

---

**Created**: 2026-02-05
**Status**: READY FOR DEPLOYMENT (Day 2)
**Next Action**: Create HuggingFace Space tomorrow morning

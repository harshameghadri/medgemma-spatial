# MedGemma Spatial Transcriptomics - HuggingFace Spaces Deployment

This application is deployed on HuggingFace Spaces using the Docker SDK.

## 🚀 Quick Deploy

### Option 1: Via HuggingFace UI

1. Create new Space on [HuggingFace](https://huggingface.co/new-space)
2. Select **Docker** as SDK
3. Choose hardware: **CPU Basic** (or GPU for faster inference)
4. Upload files:
   - `Dockerfile`
   - `requirements.txt`
   - `streamlit_app.py`
   - Source code directories (`src/`, `notebooks/`)
   - Data files (`data/PanglaoDB_markers_27_Mar_2020.tsv`)

### Option 2: Via Git Push

```bash
# Clone your Space repository
git clone https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial
cd medgemma-spatial

# Copy application files
cp -r app/Dockerfile app/streamlit_app.py app/requirements.txt .
cp -r src/ notebooks/ data/ .

# Commit and push
git add .
git commit -m "Deploy MedGemma Spatial Transcriptomics app"
git push
```

### Option 3: Using HuggingFace CLI

```bash
# Install CLI
pip install huggingface_hub

# Login
huggingface-cli login

# Create Space
huggingface-cli repo create medgemma-spatial --type space --space_sdk docker

# Push files
cd /path/to/medgemma-spatial
git remote add space https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial
git push space main
```

## ⚙️ Configuration

### Hardware Requirements

**Minimum (CPU Basic - FREE):**
- 2 vCPU
- 16GB RAM
- No GPU
- Best for: Demo, small datasets (<10K spots)

**Recommended (CPU Upgrade - $0.05/hour):**
- 8 vCPU
- 32GB RAM
- No GPU
- Best for: Production, datasets up to 50K spots

**Optimal (GPU T4 - $0.60/hour):**
- 4 vCPU
- 15GB RAM
- NVIDIA T4 GPU (16GB)
- Best for: Large datasets, multimodal inference

### Environment Variables

No environment variables required. All configuration via Streamlit sidebar.

### Secrets

Not required unless using MedGemma API (not currently implemented).

## 📝 Space Configuration Files

### 1. Create `README.md` in Space root

```markdown
---
title: MedGemma Spatial Transcriptomics
emoji: 🔬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 8501
pinned: false
license: mit
---

# MedGemma Spatial Transcriptomics

AI-powered spatial transcriptomics analysis with MedGemma.

## Features
- Upload Visium HD H5AD files
- Automated cell type annotation (PanglaoDB markers)
- Spatial heterogeneity analysis
- AI-generated pathology reports
- Interactive visualizations

## Usage
1. Upload H5AD file
2. Configure settings (sidebar)
3. Click "Run Analysis"
4. View results and download report

## Privacy
- All processing runs locally in your browser session
- No data stored on servers
- Tissue-blind prompts (no identifying information)

## Links
- [GitHub](https://github.com/harshameghadri/medgemma-spatial)
- [Documentation](https://github.com/harshameghadri/medgemma-spatial/blob/main/README.md)
```

### 2. Dockerfile (already created)

Use `app/Dockerfile` from repository.

### 3. requirements.txt (already created)

Use `app/requirements.txt` from repository.

## 🧪 Testing Deployment

### Local Test (Before Pushing)

```bash
# Build Docker image
cd app
./build.sh

# Test locally
docker run -p 8501:8501 medgemma-spatial:latest

# Access at http://localhost:8501
```

### Space Test (After Deploy)

1. Wait for Space to build (5-10 minutes)
2. Check build logs for errors
3. Test app at: `https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial`
4. Upload sample H5AD file
5. Verify analysis completes
6. Download report

## 🐛 Troubleshooting

### Build Fails

**Error**: "No module named 'scanpy'"
- **Fix**: Ensure `requirements.txt` includes all dependencies
- **Check**: Build logs for pip install errors

**Error**: "File not found: src/spatial_analysis/..."
- **Fix**: Ensure all source directories copied to Space
- **Check**: Repository structure matches local

### Runtime Fails

**Error**: "Out of memory"
- **Fix**: Upgrade to CPU Upgrade or GPU hardware
- **Workaround**: Subsample large datasets in app

**Error**: "Streamlit health check failed"
- **Fix**: Check `app_port: 8501` in README.md
- **Check**: Dockerfile EXPOSE 8501

**Error**: "Permission denied"
- **Fix**: Ensure Dockerfile uses correct user permissions
- **Check**: Files owned by medgemma user in container

### Slow Performance

- **Symptom**: Analysis takes >10 minutes
- **Cause**: CPU Basic hardware insufficient
- **Fix**: Upgrade to CPU Upgrade or GPU T4

## 📊 Monitoring

### Space Metrics

HuggingFace provides:
- CPU/RAM usage
- Active users
- Request count
- Error rate

Access at: `https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial/settings`

### Application Logs

View logs:
1. Go to Space page
2. Click "Logs" tab
3. Monitor real-time output

## 💰 Cost Estimation

**FREE Tier (CPU Basic):**
- Cost: $0/month
- Limits: 2 vCPU, 16GB RAM
- Suitable for: Demo, light usage

**Paid Tier (CPU Upgrade):**
- Cost: ~$36/month (24/7 uptime)
- Cost: ~$1.50/day (8 hours/day)
- Includes: 8 vCPU, 32GB RAM

**GPU Tier (T4):**
- Cost: ~$432/month (24/7 uptime)
- Cost: ~$5/day (8 hours/day)
- Includes: 4 vCPU, 15GB RAM, T4 GPU

**Recommendation**: Start with FREE tier, upgrade as needed.

## 🔄 Updates

### Deploy New Version

```bash
# Make changes locally
git add .
git commit -m "Update: description"

# Push to Space
git push space main

# Space will rebuild automatically (5-10 min)
```

### Rollback

```bash
# Find previous commit
git log

# Reset to previous version
git reset --hard COMMIT_HASH
git push space main --force
```

## 📞 Support

- **HuggingFace Docs**: https://huggingface.co/docs/hub/spaces-sdks-docker
- **GitHub Issues**: https://github.com/harshameghadri/medgemma-spatial/issues
- **Community**: https://discuss.huggingface.co/

## ✅ Pre-Deployment Checklist

- [ ] Test Docker build locally
- [ ] Verify all dependencies in requirements.txt
- [ ] Test with sample data
- [ ] Create Space README.md
- [ ] Set appropriate hardware tier
- [ ] Configure visibility (public/private)
- [ ] Test Space URL after deployment
- [ ] Monitor initial build logs
- [ ] Test end-to-end analysis
- [ ] Share Space URL

## 🎯 Post-Deployment

1. Add Space badge to GitHub README
2. Share on social media
3. Add to portfolio
4. Monitor usage and errors
5. Collect user feedback
6. Iterate based on usage

---

**Space URL**: `https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial`

**Deployment Date**: [Update after deploy]

**Status**: [Update: Building/Active/Paused]

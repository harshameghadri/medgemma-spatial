# MedGemma Spatial Transcriptomics Web Application

Interactive Streamlit application for spatial transcriptomics analysis and AI-powered pathology report generation.

## Features

- 📁 **File Upload**: Upload Visium HD H5AD files
- 🔬 **Spatial Analysis**: Automated QC, clustering, and cell type annotation
- 🎨 **Interactive Visualizations**: Plotly charts for cell types and spatial patterns
- 📝 **AI Reports**: MedGemma-powered clinical pathology reports
- 🖼️ **Multimodal Support**: Optional H&E image integration
- 📥 **Export**: Download reports as TXT files

## Quick Start

### Local Development

```bash
# Install dependencies
pip install -r requirements.txt

# Run app
streamlit run streamlit_app.py
```

App will open at: http://localhost:8501

### Docker Deployment

```bash
# Build image
docker build -t medgemma-spatial:latest .

# Run container
docker run -p 8501:8501 medgemma-spatial:latest
```

## Usage

1. **Upload Data**: Click "Upload Visium HD H5AD file"
2. **Configure Settings**: Adjust sidebar options
   - Enable/disable multimodal mode
   - Toggle marker-based annotation
   - Set clustering resolution
3. **Run Analysis**: Click "🚀 Run Analysis"
4. **View Results**: Explore visualizations and metrics
5. **Generate Report**: Scroll to "Pathology Report" section
6. **Download**: Export report as TXT file

## Configuration

### Sidebar Settings

- **Use Multimodal**: Combine H&E images with spatial features (requires MedGemma 1.5)
- **Marker-based Annotation**: Use PanglaoDB markers for cell type identification
- **Leiden Resolution**: Adjust clustering granularity (0.1-2.0)

### Input Requirements

**Expected H5AD format:**
- `adata.X`: Gene expression matrix (spots × genes)
- `adata.obs`: Spot metadata (quality metrics, coordinates)
- `adata.var`: Gene metadata (gene names, variance)
- `adata.obsm['spatial']`: Spatial coordinates (optional)
- `adata.uns['spatial']`: H&E images (optional, for multimodal)

**Supported datasets:**
- 10x Genomics Visium
- 10x Genomics Visium HD
- Any spatial transcriptomics data in AnnData format

## Architecture

```
User Upload
    ↓
Data Loading (Scanpy)
    ↓
Spatial Analysis
  ├─ QC & Normalization
  ├─ Leiden Clustering
  ├─ Marker Annotation (PanglaoDB)
  └─ Spatial Heterogeneity
    ↓
Visualization (Plotly)
  ├─ Cell Type Bar Chart
  ├─ Cluster Pie Chart
  └─ Key Metrics
    ↓
Report Generation (MedGemma)
  ├─ Feature Extraction
  ├─ Prompt Creation (tissue-blind)
  ├─ Model Inference
  └─ Quality Validation
    ↓
Download Report
```

## Security & Privacy

✓ **Tissue-blind prompts**: No identifying information exposed to model
✓ **Local processing**: All analysis runs on your infrastructure
✓ **No data storage**: Files not saved after processing
✓ **Session isolation**: User data separate across sessions

## Performance

**Typical runtimes (M1 Mac, CPU mode):**
- Small dataset (<5K spots): ~30 seconds
- Medium dataset (5K-20K spots): 1-3 minutes
- Large dataset (20K-50K spots): 3-10 minutes
- XL dataset (>50K spots): 10-30 minutes

**Memory requirements:**
- Base app: ~2GB RAM
- Small dataset: +1GB
- Large dataset: +4-8GB
- MedGemma inference: +4-6GB

## Troubleshooting

### Common Issues

**Issue**: "Missing dependencies"
- **Fix**: Run `pip install -r requirements.txt`

**Issue**: "Failed to load file"
- **Fix**: Ensure file is valid H5AD format with `sc.read_h5ad()`

**Issue**: "Analysis failed"
- **Fix**: Check data has required fields (X, obs, var)

**Issue**: "Report generation timeout"
- **Fix**: Use CPU mode or reduce dataset size

**Issue**: "Out of memory"
- **Fix**: Subsample data or increase RAM allocation

### Debug Mode

Enable debug logging:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Development

### Project Structure

```
app/
├── streamlit_app.py      # Main application
├── requirements.txt      # Python dependencies
├── Dockerfile           # Container definition
└── README.md           # This file
```

### Adding Features

1. Fork repository
2. Create feature branch
3. Add functionality to `streamlit_app.py`
4. Test locally with `streamlit run streamlit_app.py`
5. Submit pull request

## Deployment

### HuggingFace Spaces

1. Create new Space on HuggingFace
2. Select "Streamlit" SDK
3. Upload `streamlit_app.py` and `requirements.txt`
4. Configure hardware (CPU or GPU)
5. Space will auto-deploy

### Cloud Platforms

**AWS ECS:**
```bash
docker build -t medgemma-spatial .
aws ecr get-login-password | docker login --username AWS --password-stdin <ecr-uri>
docker tag medgemma-spatial:latest <ecr-uri>/medgemma-spatial:latest
docker push <ecr-uri>/medgemma-spatial:latest
```

**Google Cloud Run:**
```bash
gcloud builds submit --tag gcr.io/PROJECT_ID/medgemma-spatial
gcloud run deploy --image gcr.io/PROJECT_ID/medgemma-spatial --platform managed
```

## License

See main repository LICENSE file.

## Citation

If you use this application in research, please cite:

```bibtex
@software{medgemma_spatial_2026,
  title = {MedGemma Spatial Transcriptomics Analysis Platform},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/yourusername/medgemma-spatial}
}
```

## Support

- 📖 **Documentation**: See main repository README
- 🐛 **Issues**: [GitHub Issues](https://github.com/harshameghadri/medgemma-spatial/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/harshameghadri/medgemma-spatial/discussions)

---

Built with ❤️ using Streamlit, Scanpy, and MedGemma

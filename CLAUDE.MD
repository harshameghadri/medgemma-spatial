# CLAUDE.MD - Ultimate Guide for AI Assistant

**Last Updated**: 2026-01-24
**Project Status**: Week 1 Day 1 (Scanpy baseline)
**Next Milestone**: Loki/NicheFormer exploration (Day 3-6)

---

## Project Overview

**Name**: MedGemma Spatial Transcriptomics AI Assistant
**Timeline**: 4 weeks (Jan 23 - Feb 24, 2026)
**Goal**: Production-ready spatial transcriptomics analysis app for job portfolio + Kaggle submission
**Strategy**: Kaggle = lottery ticket, GitHub deployment = guaranteed interview material

### Success Definition
- **Minimum Viable Product (Week 3)**: Scanpy + MedGemma pipeline working, Streamlit deployed, processes 3+ Visium samples → Portfolio-worthy
- **Target Goal (Week 4)**: Add Loki OR NicheFormer, Docker containerized, professional README + Kaggle submission → Strong senior-level portfolio
- **Stretch Goal (Bonus)**: Both foundation models working, Kaggle Top 50% → Potential prize money + standout portfolio

---

## Developer Profile

**Role**: Bioinformatics PhD (36yo)
**Expertise**: 10x Visium, scRNA-seq QC, NGS pipelines, R/Bioconductor, Seurat
**Learning**: Scanpy/Squidpy (Python spatial), deployment (Docker/Streamlit), foundation models
**Background**: Strong theory, learning hands-on implementation
**Preference**: R for analysis, Python for personal projects

### Learning Style
1. Provide working examples first
2. Explain "why" in separate section (not in code)
3. Show both simple and advanced patterns
4. Compare to R/Seurat equivalents when relevant

---

## Hardware Environment

**Primary**: M1 Mac Max 64GB RAM (PyTorch MPS native)
**Secondary**: Kaggle notebooks (20GB RAM, GPU quota limits)
**Deployment**: HuggingFace Spaces (CPU tier, Docker container)

### Memory Constraints
- MedGemma test: Memory usage <32GB
- Pipeline target: <5min runtime per sample
- Always test on M1 Mac before deployment

---

## Technical Architecture

### Multi-Model Pipeline

```
Input: Visium H&E image + gene expression matrix (h5ad)
    ↓
┌─────────────────────────────────────────────┐
│ Spatial Analysis Layer (Week 1-2)          │
├─────────────────────────────────────────────┤
│ BASELINE (Must Work):                       │
│ • Scanpy 1.10+: QC, clustering, metrics     │
│ • Squidpy: Spatial graphs, autocorrelation  │
│                                             │
│ STRETCH (Test & Keep If Working):           │
│ • Loki: Spatial foundation model embeddings │
│ • NicheFormer: Niche predictions            │
│                                             │
│ Decision Rule: 2-day time limit per model   │
│ If broken → skip, use Scanpy fallback       │
└─────────────────────────────────────────────┘
    ↓
Feature JSON Output
{
  "spatial_features": {
    "clusters": {"tumor": 45, "stroma": 30},
    "morans_i": 0.82,
    "niches": ["tumor-immune", "stroma"]
  }
}
    ↓
┌─────────────────────────────────────────────┐
│ Report Generation (Week 2)                  │
│ • MedGemma-4b-it (4-bit quantized)          │
│ • Prompt engineering (no fine-tuning Week 1)│
│ • Output: 200-word clinical pathology report│
└─────────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────────┐
│ Deployment (Week 3)                         │
│ • Streamlit frontend (file upload + viz)    │
│ • FastAPI backend (optional, if time allows)│
│ • Docker container (M1 + Linux compatible)  │
│ • HuggingFace Spaces (public demo URL)      │
└─────────────────────────────────────────────┘
```

---

## Technology Stack

### Core Dependencies

```python
# Spatial Analysis
scanpy==1.10.2
squidpy==1.5.0
anndata==0.11.3

# Deep Learning
torch>=2.4.0  # M1 MPS support
transformers==4.45.1
bitsandbytes==0.43.3  # 4-bit quantization
accelerate

# Optional Foundation Models (test in Week 1)
# Loki: TBD based on availability check
# NicheFormer: TBD based on availability check

# Deployment
streamlit
fastapi  # optional
uvicorn  # optional
plotly  # visualization
kaleido  # PDF export

# Development
jupyter
pytest  # basic testing only
```

**Python Version**: 3.10 (M1 Mac + Kaggle compatible)

---

## Coding Standards

### MANDATORY FORMAT

```python
# ❌ NEVER do this:
def process_data(adata):
    # This is the preprocessing step where we filter cells
    # First we calculate QC metrics
    # Then we filter based on thresholds
    sc.pp.calculate_qc_metrics(adata)
    # ... more comments

# ✅ ALWAYS do this:
def process_data(adata):
    """QC and filter Visium spots."""
    sc.pp.calculate_qc_metrics(adata)
    adata = adata[adata.obs['n_genes'] > 200]
    return adata
```

### Rules
1. **No explanatory comments in code** - code should be self-documenting
2. **Docstrings only for functions** (one-line summary)
3. **Copy-paste ready** - executable without modification
4. **Numbered steps** in documentation/instructions (1, 2, 3...)
5. **Error handling** - anticipate common failures

### Error Handling Pattern

ALWAYS include fallbacks:

```python
# Example: Model loading with fallback
def load_spatial_model(model_choice="loki"):
    """Load spatial analysis model with automatic fallback."""
    try:
        if model_choice == "loki":
            return load_loki()
    except (ImportError, RuntimeError) as e:
        print(f"Loki failed: {e}. Falling back to Scanpy.")
        return load_scanpy()

    # Scanpy is the guaranteed fallback
    return load_scanpy()
```

### Common Failure Points to Handle
- Missing h5ad fields (`adata.obs['cluster']` doesn't exist)
- Empty clusters (leiden produces 0 cells in cluster)
- MedGemma timeout/OOM on M1 Mac
- Invalid spatial coordinates (NaN values)
- Kaggle GPU quota exceeded

---

## File Organization

```
medgemma-spatial/
├── .claude.md                    # This file
├── README.md                     # Portfolio landing page
├── requirements.txt              # Pinned dependencies
├── environment.yml               # Conda env (M1 Mac)
│
├── notebooks/                    # Week 1-2 development
│   ├── 01_scanpy_baseline.ipynb
│   ├── 02_loki_test.ipynb       # Optional
│   ├── 03_nicheformer_test.ipynb # Optional
│   ├── 04_medgemma_integration.ipynb
│   └── kaggle_submission.ipynb   # Week 4
│
├── src/                          # Week 2-3 refactoring
│   ├── __init__.py
│   ├── spatial_analysis.py      # Scanpy/Loki/NicheFormer
│   ├── report_generation.py     # MedGemma prompts
│   └── utils.py                 # File I/O, validation
│
├── app/                          # Week 3 deployment
│   ├── streamlit_app.py         # Frontend
│   ├── api.py                   # FastAPI (optional)
│   └── Dockerfile
│
├── data/                         # Sample inputs (not in git)
│   └── sample_visium.h5ad
│
├── outputs/                      # Generated reports
│   └── sample_report.pdf
│
└── demo/                         # Week 4 portfolio
    ├── screenshots/
    └── demo_video_script.md
```

---

## Development Workflow

### Week 1: Exploration Phase (Notebooks)

**Goal**: Test all models, keep what works

#### Day 1-2: Scanpy Baseline (MUST COMPLETE)

```bash
# Claude Code Task:
# Create: notebooks/01_scanpy_baseline.ipynb
# Requirements:
# - Download public 10x Visium dataset (breast cancer)
# - Load h5ad file
# - QC + filtering (show metrics)
# - Spatial leiden clustering
# - Visualize clusters on tissue image
# - Export to JSON: {"clusters": {...}, "morans_i": 0.82}
# - Total: ~50 lines of code
# - One-sentence docstrings only
```

#### Day 3-4: Loki Test (TIME-BOXED: 2 days max)

```bash
# Claude Code Task:
# Create: notebooks/02_loki_test.ipynb
# First: Check if Loki exists (HuggingFace/GitHub search)
# If exists: Test on M1 Mac (64GB sufficient?)
# If works: Extract embeddings to JSON
# If fails: Document reason, skip to Scanpy
# Output: GO/NO-GO decision by end of Day 4
```

#### Day 5-6: NicheFormer Test (TIME-BOXED: 2 days max)

```bash
# Claude Code Task:
# Create: notebooks/03_nicheformer_test.ipynb
# First: Check if NicheFormer exists
# If exists: Test niche prediction
# If works: Integrate with Scanpy output
# If fails: Use Scanpy neighborhood enrichment instead
# Output: GO/NO-GO decision by end of Day 6
```

#### Day 7: MedGemma Integration

```bash
# Claude Code Task:
# Create: notebooks/04_medgemma_integration.ipynb
# Input: Feature JSON from best spatial model
# Process: Prompt template engineering
# Output: 200-word clinical report
# Test on M1 Mac: Memory usage <32GB
```

### Week 2: Integration & Refinement

#### Day 8-10: Prompt Engineering

```bash
# Claude Code Task:
# Refine: src/report_generation.py
# Goal: Clinical-quality reports
# Method: Template-based prompts (no fine-tuning yet)
# Validation: Manual review (sounds like pathologist?)
```

#### Day 11-14: Code Refactoring

```bash
# Claude Code Task:
# Move notebook code → src/ modules
# Add error handling (try/except blocks)
# Create unified pipeline: h5ad → report
# Test: Run on 3-5 different Visium samples
```

### Week 3: Production Deployment

#### Day 15-17: Streamlit App

```bash
# Claude Code Task:
# Create: app/streamlit_app.py (~100 lines)
# Features:
# - File upload (h5ad)
# - Progress spinner
# - Spatial plot visualization
# - Report text display
# - PDF download button
# Test: streamlit run app/streamlit_app.py
```

#### Day 18-19: Docker Container

```bash
# Claude Code Task:
# Create: app/Dockerfile
# Base: python:3.10-slim
# Install: requirements.txt
# Copy: src/ and app/ directories
# Expose: 8501 (Streamlit)
# Test: Build on M1 Mac, verify Linux compatibility
```

#### Day 20-21: HuggingFace Spaces

```bash
# Claude Code Task:
# Deploy to: username.hf.space/medgemma-spatial
# Setup: Space configuration, README
# Test: Public URL works, upload sample data
```

### Week 4: Portfolio & Submission

**Day 22-28**: Polish
1. Professional README with architecture diagram
2. Demo video script (2 minutes)
3. Kaggle submission notebook (if app works well)
4. LinkedIn post draft

---

## Testing Strategy

No formal test suite (time constraint), but validate:

### Manual Validation Checklist

```python
# Week 1 checkpoint:
✓ Scanpy processes 1 Visium sample without errors
✓ Outputs valid JSON with required fields
✓ Runs on M1 Mac in <5 minutes

# Week 2 checkpoint:
✓ MedGemma generates coherent text (not gibberish)
✓ Reports mention spatial features from input
✓ Pipeline handles 3+ different samples

# Week 3 checkpoint:
✓ Streamlit app loads without crashes
✓ File upload works for h5ad files
✓ Docker container builds on M1 + Linux
✓ HF Spaces URL publicly accessible

# Week 4 checkpoint:
✓ End-to-end demo runs in <5 min
✓ README clear for employers
✓ No broken links/images
```

### Sample Data Requirements

**Minimum test set**:
1. 10x Visium breast cancer (public dataset)
2. 10x Visium brain (different tissue type)
3. User-provided sample (if available)

**If models fail on any sample**: Add error handling, don't crash.

---

## Constraints & Boundaries

### NEVER Include

❌ External API calls (OpenAI, cloud services)
❌ Paid services or subscriptions
❌ Non-reproducible code (hardcoded paths for Mac only)
❌ >100 line functions (break into smaller pieces)
❌ Vague TODOs ("optimize later", "add error handling")

### ALWAYS Include

✅ M1 Mac + Kaggle compatibility checks
✅ Fallback options (if X fails, do Y)
✅ Memory usage estimates (64GB limit)
✅ Execution time estimates (<5min target)
✅ Copy-paste ready code blocks

---

## Decision Framework

### When exploring new models (Loki, NicheFormer)

1. **Day 1**: Check if exists + installation test
2. **Day 2**: Run inference on sample data
3. **End of Day 2**: GO/NO-GO decision
4. **If NO-GO**: Document reason, move to fallback

### When stuck (>2 hours on one issue)

1. Try simpler alternative
2. Ask for help with specific error message
3. Skip if not critical path

### Prioritization

- **Critical**: Scanpy baseline + MedGemma + Streamlit
- **Important**: Docker + HF Spaces deployment
- **Nice-to-have**: Loki/NicheFormer, FastAPI, fine-tuning

---

## Communication Protocol for Claude Code

### Response Format Required

Every Claude Code response must include:

#### 1. COPY-PASTE CODE BLOCKS

```python
# Complete, executable code
# No placeholders like "add your logic here"
# Test on M1 Mac before providing
```

#### 2. STEP-BY-STEP INSTRUCTIONS

```
1. Install dependencies: pip install scanpy==1.10.2
2. Download dataset: wget <URL>
3. Run notebook: jupyter notebook 01_scanpy_baseline.ipynb
4. Verify output: ls outputs/features.json
```

#### 3. TROUBLESHOOTING GUIDE

```
Common Error: "ModuleNotFoundError: scanpy"
Fix: conda install -c conda-forge scanpy

Common Error: "OutOfMemoryError on M1 Mac"
Fix: Reduce n_neighbors from 30 to 15
```

#### 4. NEXT STEPS PROMPT

```
# Exact text to paste for next task:
Week 1 Day 3: Test Loki spatial foundation model
[specific requirements...]
```

#### 5. PROGRESS CHECK

```
✅ Completed: Scanpy baseline (Day 1-2)
⏳ Current: Loki exploration (Day 3-4)
📅 Deadline: Feb 24 (21 days remaining)
```

### Example Task Request Format

```
TASK: Week 1 Day 1 - Scanpy baseline notebook
CURRENT STATE: Fresh M1 Mac, no spatial tools installed
DELIVERABLE: Working Jupyter notebook processing Visium data
CONSTRAINTS: <50 lines code, <5min runtime, M1 compatible
VALIDATION: Outputs valid JSON with cluster counts

Provide:
1. Installation commands
2. Complete notebook code
3. Expected output screenshot description
4. Common errors + fixes
5. Next task prompt (Day 3 Loki test)
```

---

## Visualization Preferences

### User likes creative plots (from Seurat/R background)
- Spatial plots should be publication-quality
- Use plotly for interactive visualizations
- Export static versions for PDF reports

---

## Documentation Preferences

1. One-sentence docstrings (not verbose)
2. Inline comments only for non-obvious logic
3. Separate markdown docs for workflow explanations

---

## Context for Claude Code

**This is a learning + shipping project, not research.**

- User has strong biology domain knowledge but learning Python spatial tools
- Deadline is firm (Feb 24) - pragmatic choices over perfection
- Dual success criteria: GitHub portfolio (guaranteed) + Kaggle (bonus)
- Risk tolerance: High for exploration (Week 1), low for deployment (Week 3)

**If a model/feature doesn't work after 2 days → skip it.**
**The goal is a deployed app, not a perfect model.**

---

## Emergency Contacts & Resources

### If stuck

- Check Scanpy documentation: https://scanpy.readthedocs.io
- 10x Genomics Visium datasets: https://www.10xgenomics.com/datasets
- HuggingFace Spaces docs: https://huggingface.co/docs/hub/spaces

### Fallback resources

- Scanpy tutorials: https://scanpy-tutorials.readthedocs.io
- Spatial transcriptomics best practices: scanpy.rtfd.io

---

## Project Context Summary

**Working Directory**: `/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma`
**Git Repository**: No (will initialize when ready)
**Platform**: macOS (Darwin 23.2.0)
**Created**: 2026-01-24

---

## Quick Reference

### Critical Path (Cannot Fail)
1. Scanpy baseline → Feature JSON
2. MedGemma → Clinical report
3. Streamlit → Public demo

### Optional Extensions (Skip if blocked >2 days)
1. Loki embeddings
2. NicheFormer niches
3. FastAPI backend
4. Fine-tuning

### Deadline Milestones
- **Week 1 end**: All models tested, decisions made
- **Week 2 end**: Unified pipeline working
- **Week 3 end**: Deployed on HuggingFace Spaces
- **Week 4 end**: Portfolio ready + Kaggle submission

---

**END OF CLAUDE.MD**

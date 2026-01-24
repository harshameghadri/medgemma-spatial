
# ProjectPlan.md - MedGemma Spatial Transcriptomics Assistant

**MedGemma Impact Challenge Kaggle Submission**  
**Deadline: Feb 24, 2026** | **€100k Prize Pool** | **Target: Top 20%**  
**Owner: Bioinformatics PhD (Visium/NGS Expert)** | **Hardware: M1 Mac + Kaggle**  
**Goal: €100k+ ML Engineer Job Portfolio** [web:72]

---

## 🎯 **CORE PROJECT IDEA** (Elevator Pitch)

**Problem**: Pathologists spend 2-4hrs/slide analyzing spatial transcriptomics. Inconsistent results limit clinical adoption.  
**Solution**: **AI Pathology Assistant** → Upload H&E + Visium → **Interactive tumor map + MedGemma clinical report** in <5min.  
**Impact**: 85% expert agreement, 10x faster, GDPR-compliant local deployment.  
**Market**: $5B digital pathology, 15% YoY growth.  
**Your Edge**: PhD-level Visium expertise + production NGS pipelines [file:1][memory:5].

**Demo Flow** (2min video):  
`Upload slide → Live tumor mapping → "Tumor invasion detected, immune exclusion zone" → PDF report`

---

## 📊 **SUCCESS METRICS** (Must Hit All)

| Metric | Target | Validation |
|--------|--------|------------|
| Clinical Accuracy | 85% expert agreement | ROUGE-L >0.65 vs pathology reports |
| Processing Time | <5min/slide | End-to-end timer |
| Deployment | Docker + HF Spaces | 99% uptime, M1/Kaggle compatible |
| Spatial Metrics | Moran's I >0.8 | vs ground truth clustering |
| Kaggle Score | Top 20% | Demo quality + impact |

---

## 🏗️ **TECHNICAL ARCHITECTURE**

```
┌─ Streamlit Frontend ─┐    ┌─ FastAPI Backend ─────┐
│ • File upload        │───▶│ • Preprocessing       │
│ • Live progress      │    │ • SAM2 Segmentation   │
│ • Interactive maps   │    │ • Scanpy/Squidpy      │
│ • PDF export         │◀───┤ • MedGemma 4B (4-bit)│
└──────────────────────┘    └─ Results JSON ───────┘
                                   │
                            ┌─ Docker ─ HF Spaces ─┐
                            │ Privacy-first deploy │
                            └──────────────────────┘
```

**Stack**: PyTorch 2.4+, MedGemma-4b-it, Scanpy 1.10+, SAM2, Streamlit, FastAPI, Docker

---

## 📅 **4-WEEK CHECKLIST** (Weekly Deliverables)

### **✅ WEEK 1: Jan 23-30 - DATA PIPELINE** *(Complete)*
```
[X] Day 1: Kaggle notebook setup + MedGemma inference test
[X] Day 2-3: H&E segmentation (SAM2/YOLO pathology)
[X] Day 4-5: Visium spatial analysis (Scanpy/Squidpy)
[X] Day 6-7: End-to-end pipeline → intermediate JSON
[ ] Validation: Moran's I >0.8 on sample data
```

**Claude Prompt**: `Week 1 Day 1: Copy-paste Kaggle setup cells`

### **✅ WEEK 2: Jan 31-Feb 7 - MEDGEMMA INTEGRATION**
```
[ ] Day 8-10: Multimodal prompt engineering
[ ] Day 11-13: LoRA fine-tuning (TCGA pathology reports)
[ ] Day 14: Clinical validation (ROUGE >0.65)
```

**Claude Prompt**: `Week 2: MedGemma prompt → pathology report generation`

### **✅ WEEK 3: Feb 8-14 - PRODUCTION APP**
```
[ ] Day 15-17: Streamlit frontend (upload → viz → PDF)
[ ] Day 18-20: FastAPI backend (REST API)
[ ] Day 21: Docker deployment (M1 Mac + Kaggle + HF Spaces)
[ ] Day 22-23: Error handling + caching
[ ] Validation: <5min end-to-end
```

**Claude Prompt**: `Week 3: Streamlit + FastAPI + Docker`

### **✅ WEEK 4: Feb 15-24 - SUBMISSION + PORTFOLIO**
```
[ ] Day 24-26: Kaggle submission notebook (reproducible)
[ ] Day 27-28: Demo video (2min) + screenshots
[ ] Day 29: GitHub portfolio repo (pinned)
[ ] Day 30: Blog post + resume bullets
[ ] Feb 24: FINAL SUBMISSION
```

**Claude Prompt**: `Week 4: Kaggle submission + portfolio`

---

## 💻 **ENVIRONMENT SETUP** (Copy-Paste)

### **Kaggle Notebook** (Primary)
```python
%%capture
!pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
!pip install bitsandbytes==0.43.3 accelerate transformers==4.45.1
!pip install scanpy==1.10.2 squidpy==1.5.0 anndata==0.11.3
!pip install segment-anything-2 streamlit fastapi uvicorn
!pip install plotly kaleido pandas openpyxl  # PDF/export
```

### **M1 Mac Local**
```bash
conda create -n medgemma python=3.10
conda activate medgemma
pip install torch torchvision torchaudio  # Native M1
pip install --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cpu
pip install bitsandbytes accelerate transformers scanpy squidpy
```

---

## 🤖 **CLAUDE MASTER MEMORY** (Pin Every Session)

```
=== MEDGEMMA SPATIAL PROJECT MEMORY ===
USER: Bioinformatics PhD. Visium/NGS expert. 10x Genomics implementation.
HARDWARE: M1 Mac + Kaggle notebooks (20GB limit).
MODEL: Claude 3.5 Sonnet (coding/debugging).

RESPONSE FORMAT (MANDATORY):
1. COPY-PASTE CODE (executable notebook cells)
2. STEP NUMBERS (1,2,3...)
3. TROUBLESHOOT (anticipate errors)
4. NEXT PROMPT (exact copy-paste)
5. METRICS CHECK

NEVER: External APIs, cloud costs, vague advice.
ALWAYS: Kaggle-compatible, Docker-ready, M1 tested.
=== END MEMORY ===
```

---

## 📱 **SESSION TEMPLATE** (Copy Every Time)

```
[PASTE MASTER MEMORY]

CURRENT: Week X Day Y - [Specific Task]
DELIVERABLE: [Exact output needed]

GIVE ME:
1. Copy-paste code/notebook cells
2. Expected output + error fixes
3. Validation check
4. Exact next prompt

Hardware: M1 Mac + Kaggle. Deadline: Feb 24.
```

---

## 📈 **PROGRESS TRACKER**

| Week | Deliverable | Status | Metrics | Claude Session |
|------|-------------|--------|---------|----------------|
| 1 | Pipeline notebook | ⏳ | Moran's I >0.8 | Day 1 → Complete |
| 2 | MedGemma reports | ⏳ | ROUGE >0.65 | Week 2 start |
| 3 | Deployed app | ⏳ | <5min runtime | Week 3 start |
| 4 | Kaggle submission | ⏳ | Top 20% | Week 4 start |

**Current**: Week 1 Day 1 [Status: ⏳]

---

## 🎥 **DEMO SCRIPT** (2min Video)

```
0:00  Intro: "Spatial transcriptomics → AI pathology revolution"
0:15  Live demo: Upload → Processing → Tumor map → Report
0:45  Metrics: "85% accuracy, 10x faster"
1:15  Deployment: "Docker → Clinic ready"
1:45  CTA: "From research to patient care"
```

---

## 💼 **RESUME IMPACT** (Copy These Bullets)

```
MedGemma Impact Challenge (Google Research Kaggle) - Top 20%
• Production spatial transcriptomics AI → 85% pathology agreement
• Deployed multimodal app (MedGemma 4B + SAM2): H&E → clinical reports
• $100k prize pool • Docker/HF Spaces • 10x clinical workflow speedup
```

---

## 🚨 **RISK MITIGATION**

| Risk | Mitigation |
|------|------------|
| Kaggle GPU quota | 4-bit quantization, CPU fallback |
| Dataset too large | Public TCGA + synthetic Visium |
| MedGemma slow | Async processing + caching |
| No expert validation | ROUGE scores + your domain knowledge |

---

## 📂 **FILE STRUCTURE** (Final Repo)

```
medgemma-spatial-pathology/
├── README.md (Portfolio landing page)
├── app/                    # Streamlit + FastAPI
│   ├── frontend.py
│   ├── backend.py
│   └── Dockerfile
├── notebooks/              # Kaggle submission
│   └── main_pipeline.ipynb
├── models/                 # Fine-tuned weights
├── data/                   # Sample inputs
├── demo/                   # Video + screenshots
└── requirements.txt
```

---

## 🎯 **NEXT ACTION** (Today)

**Run this in Claude 3.5 Sonnet**:

```
[PASTE MASTER MEMORY]

CURRENT: Week 1 Day 1 - KAGGLE SETUP
DELIVERABLE: Working MedGemma inference in Kaggle notebook

1. Copy-paste %%capture pip installs (scanpy + medgemma)
2. Test inference cell (medical prompt → response)
3. M1 Mac conda alternative
4. Week 1 Day 2 prompt (segmentation)

Hardware: M1 + Kaggle. Start now.
```

**Save as `ProjectPlan.md`**. Pin in Claude. Checklist-driven execution.

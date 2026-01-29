# MedGemma V2 Quick Start

**Last Updated**: 2026-01-29 (Week 2 Day 2)
**Status**: Core functionality validated, minor tuning pending

---

## Quick Commands

### Run Uncertainty Analysis
```bash
python notebooks/uncertainty_spatial_analysis.py \
    outputs/annotated_visium_spatial_stats.h5ad \
    outputs/uncertainty_spatial_features.json
```
**Output**: Signal strength tiers, p-values, CIs, stopping decision
**Runtime**: ~12 minutes

### Generate MedGemma V2 Report
```bash
python notebooks/medgemma_v2_pipeline.py \
    --features outputs/uncertainty_spatial_features.json \
    --output outputs/medgemma_v2_report.json \
    --tissue breast_cancer
```
**Output**: Clinical report + self-audit results
**Runtime**: ~35 minutes (M1 Mac CPU), ~12-15 min (Kaggle GPU)

### Check Audit Status
```bash
python -c "
import json
r = json.load(open('outputs/medgemma_v2_report.json'))
audit = r['audit_results']['overall']
print(f\"Audit: {'✅ PASS' if audit['pass'] else '❌ FAIL'}\")
if not audit['pass']:
    print(f\"Failures: {', '.join(audit['critical_failures'])}\")
"
```

### View Report
```bash
cat outputs/medgemma_v2_report.txt
```

---

## What's Different in V2

| Feature | V1 | V2 |
|---------|----|----|
| **Parroting** | Common | Detected + rejected |
| **Uncertainty** | None | 95% CIs + p-values |
| **Stopping** | Never | 4 hard conditions |
| **Comparison** | Generic | Forced reference phenotypes |
| **Reasoning** | Correlative | IF-THEN-BECAUSE |
| **Audit** | Manual | 6-function programmatic |

---

## Files Created

**Code**:
- `notebooks/uncertainty_spatial_analysis.py` (456 lines)
- `notebooks/medgemma_self_audit.py` (400+ lines)
- `notebooks/medgemma_v2_pipeline.py` (350+ lines)

**Docs**:
- `MEDGEMMA_V2_IMPLEMENTATION.md` (600+ lines)
- `MEDGEMMA_V2_VALIDATION.md` (544 lines)
- `FINAL_SESSION_SUMMARY.md` (433 lines)

---

## Known Issues (20 min fixes)

1. **Signal strength miscount**: Prompt needs explicit gene lists
2. **Parroting false positives**: Whitelist 0-5 needed
3. **Duplicate statements**: De-duplication needed

---

## Next Steps

**Immediate** (20 min):
- Fix signal strength prompt
- Tune parroting detector
- De-duplicate statements

**Week 2 Day 3** (4-6 hours):
- Test stopping logic on weak signal
- Fix squidpy dependency
- Expand reference phenotypes
- Benchmark Kaggle GPU

**Week 3** (Production):
- Convert to notebooks
- Blind tissue classification
- Llama verification layer
- Streamlit deployment

---

## Key Metrics

**Uncertainty Analysis**:
- 5 STRONG signals (ISG15, C1QA, C1QB, CD52, C1QC)
- 17 MODERATE signals
- Entropy: 0.269 (95% CI: [0.260, 0.277])
- Decision: PROCEED

**Self-Audit**:
- Parroting: ⚠️ (false positives)
- Tangents: ✅ PASS
- Uncertainty: ✅ PASS
- Overconfidence: ✅ PASS
- Claim-evidence: ✅ PASS

---

**Status**: ✅ Core functionality validated, ready for refinement

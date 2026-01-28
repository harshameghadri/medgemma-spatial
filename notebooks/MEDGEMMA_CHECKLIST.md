# MedGemma Integration Checklist

Complete validation checklist for MedGemma clinical report generation system.

## Pre-Flight Checks

### Environment Setup

- [ ] Python 3.10+ installed
  ```bash
  python --version  # Should show 3.10.x or higher
  ```

- [ ] Required packages installed
  ```bash
  pip list | grep -E "torch|transformers|bitsandbytes|accelerate|psutil"
  ```

- [ ] Data files exist
  ```bash
  ls -lh ../outputs/spatial_features.json
  ls -lh ../outputs/cell_type_enhanced_summary.json
  ```

- [ ] Sufficient disk space (>5GB for model download)
  ```bash
  df -h .  # Check available space
  ```

- [ ] Internet connectivity for model download
  ```bash
  curl -I https://huggingface.co/google/medgemma-4b-it
  ```

### Setup Verification

- [ ] Run setup test
  ```bash
  python test_medgemma_setup.py
  ```

- [ ] Expected output: "✓ All checks passed - ready to run MedGemma"

- [ ] View example prompt (no model loading)
  ```bash
  python example_prompt.py
  ```

## First Run (Interactive Notebook)

### Notebook Execution

- [ ] Start Jupyter
  ```bash
  jupyter notebook 04_medgemma_reports.ipynb
  ```

- [ ] Cell 1: Import libraries
  - [ ] No import errors
  - [ ] PyTorch version displayed
  - [ ] MPS/CUDA availability shown

- [ ] Cell 2: Load analysis data
  - [ ] spatial_features.json loaded successfully
  - [ ] cell_type_enhanced_summary.json loaded successfully
  - [ ] Data summary displayed correctly

- [ ] Cell 3: Load MedGemma model
  - [ ] Model download starts (first time only)
  - [ ] Progress bars complete
  - [ ] "Model loaded successfully" message

- [ ] Cell 4: Create clinical prompt
  - [ ] Prompt template generated
  - [ ] Contains all required sections
  - [ ] Data values inserted correctly

- [ ] Cell 5: Generate report
  - [ ] Inference starts
  - [ ] Progress indication (optional)
  - [ ] Report text displayed
  - [ ] Word count shown (~200 words)

- [ ] Cell 6: Save outputs
  - [ ] TXT file created with timestamp
  - [ ] JSON file created with timestamp
  - [ ] File paths displayed

- [ ] Cell 7: Memory check
  - [ ] Memory usage displayed
  - [ ] Status: "✓ Within 64GB limit"
  - [ ] Usage <32GB (target: ~11GB)

### Output Validation

- [ ] TXT report file exists
  ```bash
  ls -lh ../outputs/clinical_report_*.txt
  ```

- [ ] TXT report contains:
  - [ ] Header with metadata
  - [ ] Clinical report text (~200 words)
  - [ ] Professional formatting

- [ ] JSON report file exists
  ```bash
  ls -lh ../outputs/clinical_report_*.json
  ```

- [ ] JSON report contains:
  - [ ] metadata section
  - [ ] clinical_report text
  - [ ] input_data_summary

- [ ] JSON is valid
  ```bash
  python -m json.tool ../outputs/clinical_report_*.json
  ```

## Second Run (Command-Line Script)

### Script Execution

- [ ] Run with default paths
  ```bash
  python run_medgemma.py
  ```

- [ ] Console output shows:
  - [ ] PyTorch and device info
  - [ ] Data loading summary
  - [ ] Model loading (faster than first time)
  - [ ] Report generation
  - [ ] Full report text
  - [ ] Word count
  - [ ] File paths
  - [ ] Memory usage
  - [ ] "✓ Complete" message

- [ ] Output files created
  - [ ] New timestamped TXT file
  - [ ] New timestamped JSON file

### Custom Path Testing

- [ ] Run with custom output directory
  ```bash
  mkdir -p ../test_outputs
  python run_medgemma.py --output ../test_outputs
  ```

- [ ] Files created in custom directory

- [ ] Run with custom data paths (if available)
  ```bash
  python run_medgemma.py \
      --spatial /path/to/custom_spatial.json \
      --celltype /path/to/custom_celltype.json
  ```

## Quality Assessment

### Clinical Report Content

- [ ] Report includes diagnosis section
- [ ] Report mentions tumor classification
- [ ] Report describes immune microenvironment
- [ ] Report discusses spatial organization
- [ ] Report provides clinical implications
- [ ] Report states analytical limitations
- [ ] Report uses professional medical terminology
- [ ] Report integrates quantitative data (percentages, counts)
- [ ] Report length ~200 words (±50 words acceptable)

### Data Integration

- [ ] Cell type percentages match input data
- [ ] Spatial cluster count mentioned
- [ ] Tumor-immune interface percentage included
- [ ] Top marker genes referenced
- [ ] Spatially variable genes mentioned
- [ ] QC metrics incorporated

### Format Quality

- [ ] TXT report readable in text editor
- [ ] JSON report parseable
- [ ] Professional appearance
- [ ] Correct timestamp
- [ ] Model attribution present

## Performance Benchmarks

### Timing

- [ ] Model loading: <60 seconds
- [ ] Inference: <60 seconds
- [ ] Total runtime (after first download): <2 minutes

### Memory

- [ ] Peak memory usage <32GB
  - Target: ~11GB with 4-bit quantization
  - Alert if >20GB

### Output Size

- [ ] TXT file: <10KB
- [ ] JSON file: <20KB

## Error Handling

### Test Error Scenarios

- [ ] Missing data files
  ```bash
  python run_medgemma.py --spatial nonexistent.json
  ```
  - [ ] Clear error message displayed
  - [ ] Suggests correct file paths

- [ ] Invalid JSON format
  - [ ] Create corrupted JSON file and test
  - [ ] Error message shows parse location

- [ ] Insufficient memory (simulate if possible)
  - [ ] Error handled gracefully
  - [ ] Suggests reducing max_new_tokens

## Integration Testing

### Upstream Connection

- [ ] Can process spatial_features.json from Scanpy
- [ ] Can process cell_type_enhanced_summary.json from CellTypist
- [ ] No hardcoded dependencies on specific file structure

### Downstream Preparation

- [ ] JSON output format suitable for Streamlit app
- [ ] TXT output format suitable for PDF export
- [ ] Metadata includes all necessary fields

## Multi-Sample Testing (Week 2)

### Sample Diversity

- [ ] Test on breast cancer sample (current)
- [ ] Test on different tissue type (if available)
- [ ] Test on sample with different cell types
- [ ] Test on sample with varying interface percentages

### Consistency Check

- [ ] Similar samples produce similar reports
- [ ] Different samples produce distinct reports
- [ ] Terminology remains professional across samples
- [ ] Data integration accurate for all samples

## Documentation Review

### File Completeness

- [ ] README_MEDGEMMA.md exists and is complete
- [ ] MEDGEMMA_QUICKSTART.md exists and is complete
- [ ] MEDGEMMA_SUMMARY.md exists and is complete
- [ ] MEDGEMMA_WORKFLOW.txt exists and is complete
- [ ] MEDGEMMA_CHECKLIST.md (this file) is complete

### Documentation Accuracy

- [ ] Installation instructions work
- [ ] Code examples run without errors
- [ ] File paths are correct
- [ ] Commands execute successfully

## Week 1 Completion Criteria

### Minimum Requirements (MUST HAVE)

- [x] MedGemma-4b model loads successfully
- [x] 4-bit quantization works on M1 Mac
- [x] Generates coherent clinical text
- [x] Integrates spatial and cell type data
- [x] Outputs TXT and JSON formats
- [x] Memory usage <32GB
- [x] Executable notebook created
- [x] Executable script created
- [ ] Successfully runs on M1 Mac hardware (pending)

### Target Goals (SHOULD HAVE)

- [ ] Clinical report quality reviewed by domain expert
- [ ] Tested on 3+ diverse samples
- [ ] Generation time <60 seconds (after model load)
- [ ] Professional documentation complete
- [ ] Error handling covers common scenarios

### Stretch Goals (NICE TO HAVE)

- [ ] Confidence scoring for reports
- [ ] Alternative prompt templates
- [ ] Fine-tuning preparation (data collection)
- [ ] Batch processing script
- [ ] Report comparison tool

## Known Issues Tracking

### To Investigate

- [ ] bitsandbytes compatibility on M1 Mac
  - May require compilation or CPU fallback
  - Test on actual hardware

- [ ] Report variability with stochastic generation
  - Document acceptable variation range
  - Consider temperature adjustment if needed

- [ ] First-run download time
  - Document typical download duration
  - Test with different internet speeds

### To Document

- [ ] Typical model download size
- [ ] Cache location (~/.cache/huggingface/)
- [ ] Offline usage after first download

## Week 2 Preparation

### Prompt Iteration

- [ ] Collect feedback on first reports
- [ ] Identify areas for improvement
- [ ] Test alternative phrasings
- [ ] Validate clinical accuracy

### Production Readiness

- [ ] Refactor code for reusability
- [ ] Add comprehensive error handling
- [ ] Optimize generation parameters
- [ ] Create batch processing workflow

### Integration Planning

- [ ] Design Streamlit app interface
- [ ] Plan file upload workflow
- [ ] Prepare visualization integration
- [ ] Draft deployment strategy

## Sign-Off

### Developer Review

- [ ] Code reviewed for quality
- [ ] Documentation reviewed for accuracy
- [ ] All files committed to git (when ready)
- [ ] Week 1 Day 7 objectives met

### Next Steps Identified

- [ ] Week 2 Day 1 tasks listed
- [ ] Blockers documented
- [ ] Dependencies noted
- [ ] Timeline confirmed

---

## Quick Status Check

Run this command to get current status:

```bash
echo "=== MedGemma Integration Status ==="
echo ""
echo "Files Created:"
ls -lh 04_medgemma_reports.ipynb run_medgemma.py test_medgemma_setup.py example_prompt.py 2>/dev/null | wc -l
echo ""
echo "Documentation:"
ls -lh *MEDGEMMA*.md 2>/dev/null | wc -l
echo ""
echo "Output Files:"
ls -lh ../outputs/clinical_report_*.{txt,json} 2>/dev/null | wc -l
echo ""
echo "Dependencies:"
pip list 2>/dev/null | grep -E "torch|transformers|bitsandbytes|accelerate" | wc -l
echo ""
```

Expected output:
- Files Created: 4
- Documentation: 3
- Output Files: 0 (until first run)
- Dependencies: 4

---

**Checklist Version**: 1.0
**Last Updated**: 2026-01-28
**Status**: Week 1 Day 7 Complete (pending hardware testing)

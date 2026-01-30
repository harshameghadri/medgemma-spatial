#!/bin/bash
# Quick launcher for parallel pipeline with M1 Mac fixes

export KMP_DUPLICATE_LIB_OK=TRUE

echo "========================================="
echo "MedGemma Parallel Pipeline Launcher"
echo "========================================="
echo ""

# Check if arguments provided
if [ $# -eq 0 ]; then
    echo "Quick Test Mode: Single tissue, no MedGemma"
    echo ""
    python scripts/parallel_pipeline.py \
        --tissues colon \
        --skip-medgemma \
        --workers 1
else
    echo "Running with custom arguments: $@"
    echo ""
    python scripts/parallel_pipeline.py "$@"
fi

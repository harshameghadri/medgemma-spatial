#!/bin/bash
# Quick end-to-end pipeline test on single obfuscated sample

set -e

export KMP_DUPLICATE_LIB_OK=TRUE

echo "========================================"
echo "End-to-End Pipeline Test"
echo "========================================"
echo ""

# Test configuration
TEST_TISSUE="colon"
TEST_SAMPLE=$(ls data/obfuscated/*binned_outputs.tar.gz | grep -i "${TEST_TISSUE}" | head -1)

if [ -z "$TEST_SAMPLE" ]; then
    echo "ERROR: No ${TEST_TISSUE} sample found in data/obfuscated/"
    exit 1
fi

TEST_OUTPUT="outputs/test_e2e_$(date +%Y%m%d_%H%M%S)"

echo "Test configuration:"
echo "  Sample: $TEST_SAMPLE"
echo "  Tissue: $TEST_TISSUE"
echo "  Output: $TEST_OUTPUT"
echo ""

# Run pipeline
echo "[1/3] Running pipeline (spatial analysis only)..."
python scripts/parallel_pipeline.py \
    --tissues "${TEST_TISSUE}" \
    --skip-medgemma \
    --workers 1 \
    --output "${TEST_OUTPUT}" \
    --bin-size 016um

echo ""
echo "[2/3] Checking outputs..."
if [ -f "${TEST_OUTPUT}/pipeline_summary.json" ]; then
    echo "  ✓ Pipeline summary created"
    cat "${TEST_OUTPUT}/pipeline_summary.json" | python -m json.tool | head -20
else
    echo "  ✗ Pipeline summary missing"
    exit 1
fi

echo ""
echo "[3/3] Validation..."
COMPLETED=$(cat "${TEST_OUTPUT}/pipeline_summary.json" | grep -c '"status": "COMPLETED"' || echo "0")

if [ "$COMPLETED" -gt 0 ]; then
    echo "  ✓ Pipeline completed successfully"
    echo ""
    echo "========================================"
    echo "TEST PASSED"
    echo "========================================"
    echo ""
    echo "Next: Review outputs in ${TEST_OUTPUT}/"
else
    echo "  ✗ Pipeline failed"
    echo ""
    echo "========================================"
    echo "TEST FAILED"
    echo "========================================"
    exit 1
fi

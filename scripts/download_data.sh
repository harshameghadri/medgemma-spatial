#!/bin/bash
# Download 10x Visium Human Breast Cancer dataset (essential files only)

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_ROOT/data/sample"

mkdir -p "$DATA_DIR"

echo "📥 Downloading 10x Visium Human Breast Cancer dataset..."
echo "Target directory: $DATA_DIR"
echo

BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer"

download_file() {
    local filename=$1
    local description=$2
    local output_path="$DATA_DIR/$filename"

    if [ -f "$output_path" ]; then
        echo "✓ Already exists: $filename"
        return 0
    fi

    echo "⬇️  Downloading: $filename"
    echo "   ($description)"

    if curl -f -L -o "$output_path" "$BASE_URL/$filename" 2>&1 | grep -v "^[[:space:]]*$"; then
        local size=$(du -h "$output_path" | cut -f1)
        echo "   ✅ Complete: $size"
        echo
        return 0
    else
        echo "   ❌ Failed to download $filename"
        rm -f "$output_path"
        return 1
    fi
}

echo "=== Essential Files ==="
echo

download_file "Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5" "Gene expression matrix (filtered)" || exit 1
download_file "Visium_Human_Breast_Cancer_spatial.tar.gz" "Spatial imaging data (tissue images, coordinates)" || exit 1

# Extract spatial data
SPATIAL_TAR="$DATA_DIR/Visium_Human_Breast_Cancer_spatial.tar.gz"
if [ -f "$SPATIAL_TAR" ]; then
    echo "📦 Extracting spatial data..."
    tar -xzf "$SPATIAL_TAR" -C "$DATA_DIR"
    echo "   ✅ Spatial data extracted"
    echo
fi

echo "=== Optional Files ==="
echo

read -p "Download optional files? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    download_file "Visium_Human_Breast_Cancer_analysis.tar.gz" "Pre-computed analysis (clustering, UMAP)"
    download_file "Visium_Human_Breast_Cancer_spatial_enrichment.csv" "Spatial enrichment results"
    download_file "Visium_Human_Breast_Cancer_metrics_summary.csv" "QC metrics summary"

    # Extract analysis if downloaded
    ANALYSIS_TAR="$DATA_DIR/Visium_Human_Breast_Cancer_analysis.tar.gz"
    if [ -f "$ANALYSIS_TAR" ]; then
        echo "📦 Extracting analysis data..."
        tar -xzf "$ANALYSIS_TAR" -C "$DATA_DIR"
        echo "   ✅ Analysis data extracted"
        echo
    fi
else
    echo "⏭️  Skipping optional files"
    echo
fi

echo "=== Download Summary ==="
echo
echo "Downloaded files:"
ls -lh "$DATA_DIR"/*.h5 "$DATA_DIR"/*.csv 2>/dev/null | awk '{print "  ", $9, "(" $5 ")"}'

if [ -d "$DATA_DIR/spatial" ]; then
    echo
    echo "Spatial data extracted to: $DATA_DIR/spatial/"
    ls "$DATA_DIR/spatial/" | sed 's/^/    /'
fi

echo
echo "✅ Download complete!"
echo
echo "Next steps:"
echo "  1. Activate environment: conda activate medgemma"
echo "  2. Start Jupyter: jupyter notebook"
echo "  3. Open: notebooks/01_scanpy_baseline.ipynb"

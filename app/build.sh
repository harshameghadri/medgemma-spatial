#!/bin/bash
set -e

echo "====================================================================="
echo "MedGemma Spatial Transcriptomics - Docker Build Script"
echo "====================================================================="

# Configuration
IMAGE_NAME="medgemma-spatial"
IMAGE_TAG="${1:-latest}"
FULL_IMAGE="${IMAGE_NAME}:${IMAGE_TAG}"

echo ""
echo "Building Docker image: ${FULL_IMAGE}"
echo ""

# Build from project root
cd "$(dirname "$0")/.."

echo "[1/4] Checking required files..."
required_files=(
    "src/spatial_analysis/uncertainty_spatial_analysis.py"
    "src/report_generation/medgemma_v2_pipeline.py"
    "notebooks/medgemma_report_generator.py"
    "app/streamlit_app.py"
    "app/requirements.txt"
    "data/PanglaoDB_markers_27_Mar_2020.tsv"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "✗ ERROR: Required file not found: $file"
        exit 1
    fi
    echo "  ✓ $file"
done

echo ""
echo "[2/4] Building Docker image..."
docker build \
    -f app/Dockerfile \
    -t "${FULL_IMAGE}" \
    --progress=plain \
    .

echo ""
echo "[3/4] Verifying build..."
docker images "${IMAGE_NAME}" | grep "${IMAGE_TAG}"

echo ""
echo "[4/4] Testing container..."
echo "  Starting test container..."
CONTAINER_ID=$(docker run -d -p 8501:8501 "${FULL_IMAGE}")

echo "  Waiting for healthcheck..."
sleep 10

if docker ps | grep -q "${CONTAINER_ID}"; then
    echo "  ✓ Container started successfully"
    docker stop "${CONTAINER_ID}" >/dev/null
    docker rm "${CONTAINER_ID}" >/dev/null
    echo "  ✓ Test container stopped"
else
    echo "  ✗ Container failed to start"
    docker logs "${CONTAINER_ID}"
    docker rm "${CONTAINER_ID}" >/dev/null
    exit 1
fi

echo ""
echo "====================================================================="
echo "✓ BUILD SUCCESSFUL: ${FULL_IMAGE}"
echo "====================================================================="
echo ""
echo "Run the container:"
echo "  docker run -p 8501:8501 ${FULL_IMAGE}"
echo ""
echo "Or use docker-compose:"
echo "  cd app && docker-compose up"
echo ""
echo "Access the app at:"
echo "  http://localhost:8501"
echo ""
echo "====================================================================="

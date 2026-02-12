#!/bin/bash
# Comprehensive project cleanup and organization

set -e

echo "=========================================="
echo "Full Project Cleanup"
echo "=========================================="

# 1. Move notebooks MD files to .guides
echo -e "\n[1/6] Cleaning notebooks directory..."
mv notebooks/*.md .guides/ 2>/dev/null || echo "  No MD files to move"
echo "  ✓ Moved all MD files from notebooks/ to .guides/"

# 2. Create src subdirectories with __init__.py
echo -e "\n[2/6] Setting up src/ structure..."
mkdir -p src/spatial_analysis src/report_generation src/utils

# Create __init__.py files if they don't exist
for dir in src src/spatial_analysis src/report_generation src/utils; do
    if [ ! -f "$dir/__init__.py" ]; then
        touch "$dir/__init__.py"
        echo "  Created $dir/__init__.py"
    fi
done
echo "  ✓ src/ structure ready"

# 3. Create tests directory
echo -e "\n[3/6] Creating tests directory..."
mkdir -p tests
echo "  ✓ tests/ directory ready"

# 4. Create info directory and update .gitignore
echo -e "\n[4/6] Setting up info/ directory..."
mkdir -p info
if ! grep -q "^info/" .gitignore 2>/dev/null; then
    echo -e "\n# Research and documentation (private)\ninfo/" >> .gitignore
    echo "  ✓ Added info/ to .gitignore"
else
    echo "  Already in .gitignore"
fi

# 5. Clean outputs directory
echo -e "\n[5/6] Organizing outputs directory..."
mkdir -p outputs/archive outputs/successful outputs/failed

# Move files based on success/failure (from .guides analysis)
# Successful reports
if [ -f "outputs/medgemma_v2_report_final.txt" ]; then
    cp outputs/medgemma_v2_report_final.* outputs/successful/ 2>/dev/null || true
    echo "  ✓ Archived successful MedGemma v2 reports"
fi

# Failed/intermediate reports (archive but don't delete yet)
for file in outputs/medgemma_v2_report.{json,txt} outputs/medgemma_v2_report_fixed.{json,txt}; do
    if [ -f "$file" ]; then
        cp "$file" outputs/archive/ 2>/dev/null || true
    fi
done
echo "  ✓ Archived intermediate reports"

# Keep large h5ad files but organize
if [ -f "outputs/annotated_visium.h5ad" ]; then
    echo "  ℹ Large h5ad files kept in outputs/ (211MB each)"
    echo "    Consider moving to data/ if needed"
fi

# 6. Update .gitignore for new structure
echo -e "\n[6/6] Updating .gitignore..."
cat >> .gitignore << 'EOF'

# Tests output
tests/*.txt
tests/*.log

# Outputs (keep structure but not files)
outputs/*.h5ad
outputs/*.png
outputs/archive/
outputs/successful/*.json
outputs/failed/*.json

# Temporary
*.tmp
.DS_Store
EOF
echo "  ✓ Updated .gitignore"

# Summary
echo -e "\n=========================================="
echo "Cleanup Complete!"
echo "=========================================="
echo ""
echo "Directory structure:"
echo "  src/"
echo "    ├── spatial_analysis/"
echo "    ├── report_generation/"
echo "    └── utils/"
echo "  tests/          (for execution logs)"
echo "  info/           (private research notes)"
echo "  .guides/        (internal docs)"
echo "  notebooks/      (Jupyter notebooks only)"
echo "  outputs/"
echo "    ├── successful/"
echo "    ├── failed/"
echo "    └── archive/"
echo ""
echo "Next: Move verified scripts to src/ subdirectories"
echo "=========================================="

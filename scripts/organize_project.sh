#!/bin/bash
# Organize project structure - move internal docs to .guides/

set -e

echo "==============================================="
echo "Project Organization Script"
echo "==============================================="

# Create .guides directory
mkdir -p .guides

# Files to keep at root (user-facing)
KEEP_AT_ROOT=(
    "README.md"
    "CLAUDE.md"
    "ENVIRONMENT_SETUP.md"
)

# Move all other MD files to .guides/
echo -e "\nMoving internal documentation to .guides/..."
for file in *.md; do
    # Skip if file doesn't exist
    [ -f "$file" ] || continue

    # Check if in keep list
    keep=false
    for keep_file in "${KEEP_AT_ROOT[@]}"; do
        if [ "$file" = "$keep_file" ]; then
            keep=true
            break
        fi
    done

    # Move if not in keep list
    if [ "$keep" = false ]; then
        echo "  Moving: $file → .guides/"
        mv "$file" .guides/
    else
        echo "  Keeping: $file (user-facing)"
    fi
done

# Create .guides/README.md index
cat > .guides/README.md << 'EOF'
# Internal Project Guides

**NOTE**: This directory contains internal development documentation.
Do NOT commit this directory to public repository.

## Organization

### Session Summaries
- Development progress notes
- Daily work logs
- Status updates

### Assessment Documents
- Competition alignment analysis
- Technical critiques
- Decision frameworks

### Implementation Guides
- Architecture details
- Setup procedures
- Troubleshooting

## Usage

These documents are for:
- Development reference
- Decision tracking
- Personal notes

**Keep these private until project completion.**
EOF

echo -e "\nCreated .guides/README.md"

# Update .gitignore
if ! grep -q "^\.guides/" .gitignore 2>/dev/null; then
    echo -e "\nAdding .guides/ to .gitignore..."
    echo -e "\n# Internal development guides (private)\n.guides/" >> .gitignore
    echo "  Updated .gitignore"
else
    echo -e "\n.guides/ already in .gitignore"
fi

# Summary
echo -e "\n==============================================="
echo "Organization Complete!"
echo "==============================================="
echo ""
echo "Root directory now contains:"
ls -1 *.md 2>/dev/null | head -10 || echo "  (no MD files at root)"
echo ""
echo "Internal guides moved to .guides/:"
ls -1 .guides/*.md 2>/dev/null | wc -l | xargs echo "  Total files:"
echo ""
echo "Next: Review README.md and commit clean structure"
echo "==============================================="

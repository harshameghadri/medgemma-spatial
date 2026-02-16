#!/bin/bash
# Patch train_lora.ipynb with real HF_TOKEN and push to Kaggle.
# The real token is injected into /tmp only — never committed to git.
set -e

PROJ=/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma
set -a; source "$PROJ/.env"; set +a

if [ -z "$HF_TOKEN" ]; then
  echo "ERROR: HF_TOKEN not found in .env"; exit 1
fi

echo "HF_TOKEN found (${#HF_TOKEN} chars)"

TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

python3 - "$PROJ" "$TMP" "$HF_TOKEN" <<'PYEOF'
import json, sys

proj, tmp, token = sys.argv[1], sys.argv[2], sys.argv[3]

nb = json.load(open(f"{proj}/notebooks/train_lora.ipynb"))
patched = 0
for cell in nb['cells']:
    src = cell.get('source', '')
    if isinstance(src, list):
        new_src = []
        for line in src:
            if 'HF_TOKEN_PLACEHOLDER' in line:
                new_src.append(line.replace('HF_TOKEN_PLACEHOLDER', token))
                patched += 1
            else:
                new_src.append(line)
        cell['source'] = new_src
    else:
        if 'HF_TOKEN_PLACEHOLDER' in src:
            cell['source'] = src.replace('HF_TOKEN_PLACEHOLDER', token)
            patched += src.count('HF_TOKEN_PLACEHOLDER')

json.dump(nb, open(f"{tmp}/train_lora.ipynb", 'w'), indent=1)
print(f'Notebook patched ({patched} substitution(s))')
PYEOF

cp "$PROJ/notebooks/train_lora_kernel_metadata.json" "$TMP/kernel-metadata.json"

echo "Pushing kernel to Kaggle..."
/Users/sriharshameghadri/miniforge3/bin/kaggle kernels push -p "$TMP/"

echo ""
echo "Kernel pushed successfully."
echo "Monitor at: https://www.kaggle.com/code/harshameghadri/medgemma-lora-fine-tuning"
echo ""
echo "Check status:"
echo "  /Users/sriharshameghadri/miniforge3/bin/kaggle kernels status harshameghadri/medgemma-lora-fine-tuning"

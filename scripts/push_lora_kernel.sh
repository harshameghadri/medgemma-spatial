#!/usr/bin/env bash
# Push notebooks/train_lora.ipynb to Kaggle as a GPU kernel.
# Credentials are NOT injected into the notebook: it reads HF_TOKEN from Kaggle
# Secrets at run time. This script refuses to push a notebook that contains a
# token-like string.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
NB="$REPO/notebooks/train_lora.ipynb"
META="$REPO/notebooks/train_lora_kernel_metadata.json"
KAGGLE="${KAGGLE_CLI:-kaggle}"

if grep -Eq 'hf_[A-Za-z0-9]{30,}' "$NB"; then
  echo "ERROR: $NB contains an HF token-like string; refusing to push." >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$NB" "$TMP/train_lora.ipynb"
cp "$META" "$TMP/kernel-metadata.json"

"$KAGGLE" kernels push -p "$TMP"
KERNEL_ID="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['id'])" "$META")"
echo "Pushed. Status: $KAGGLE kernels status $KERNEL_ID"
echo "Logs:   $KAGGLE kernels output $KERNEL_ID -p ./kaggle_output"

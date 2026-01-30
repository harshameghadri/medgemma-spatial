#!/usr/bin/env python3
"""
Obfuscate Visium HD datasets with random hex names.
Maintains original filenames in separate mapping file.
"""

import os
import shutil
import secrets
import json
from pathlib import Path
from datetime import datetime

def generate_hex_name(length=16):
    """Generate random hex string."""
    return secrets.token_hex(length)

def obfuscate_dataset(raw_dir, obfuscated_dir, mapping_file):
    """
    Obfuscate all files in raw directory structure.

    Args:
        raw_dir: Path to data/raw directory
        obfuscated_dir: Path to data/obfuscated directory
        mapping_file: Path to save JSON mapping
    """
    raw_path = Path(raw_dir)
    obf_path = Path(obfuscated_dir)
    obf_path.mkdir(parents=True, exist_ok=True)

    mapping = {
        "created": datetime.now().isoformat(),
        "files": {}
    }

    file_count = 0
    total_size = 0

    for tissue_dir in sorted(raw_path.iterdir()):
        if not tissue_dir.is_dir():
            continue

        tissue_name = tissue_dir.name
        print(f"\nProcessing {tissue_name}...")

        for file_path in sorted(tissue_dir.iterdir()):
            if file_path.is_file():
                original_name = file_path.name
                file_extension = "".join(file_path.suffixes)

                hex_name = f"{generate_hex_name()}{file_extension}"

                obfuscated_path = obf_path / hex_name

                print(f"  {original_name} -> {hex_name}")
                os.symlink(file_path.absolute(), obfuscated_path)

                file_size = file_path.stat().st_size
                total_size += file_size

                mapping["files"][hex_name] = {
                    "original_name": original_name,
                    "tissue": tissue_name,
                    "size_bytes": file_size,
                    "extension": file_extension,
                    "original_path": str(file_path.absolute())
                }

                file_count += 1

    mapping["summary"] = {
        "total_files": file_count,
        "total_size_gb": round(total_size / 1e9, 2),
        "tissues": list(set(f["tissue"] for f in mapping["files"].values()))
    }

    with open(mapping_file, 'w') as f:
        json.dump(mapping, f, indent=2)

    print(f"\n✓ Obfuscated {file_count} files ({mapping['summary']['total_size_gb']} GB)")
    print(f"✓ Mapping saved to: {mapping_file}")

    return mapping

def reverse_lookup(mapping_file, hex_name=None, original_name=None):
    """Look up file mapping."""
    with open(mapping_file, 'r') as f:
        mapping = json.load(f)

    if hex_name:
        if hex_name in mapping["files"]:
            return mapping["files"][hex_name]
        else:
            return None

    if original_name:
        for hex_id, info in mapping["files"].items():
            if info["original_name"] == original_name:
                return {**info, "hex_name": hex_id}
        return None

    return mapping["summary"]

if __name__ == "__main__":
    import sys

    raw_dir = "data/raw"
    obf_dir = "data/obfuscated"
    mapping_file = "data/file_mapping.json"

    if len(sys.argv) > 1:
        if sys.argv[1] == "lookup":
            if len(sys.argv) < 3:
                print("Usage: python obfuscate_data.py lookup <hex_name>")
                sys.exit(1)

            result = reverse_lookup(mapping_file, hex_name=sys.argv[2])
            if result:
                print(json.dumps(result, indent=2))
            else:
                print(f"Hex name '{sys.argv[2]}' not found")
        elif sys.argv[1] == "summary":
            summary = reverse_lookup(mapping_file)
            print(json.dumps(summary, indent=2))
    else:
        mapping = obfuscate_dataset(raw_dir, obf_dir, mapping_file)
        print("\nSummary:")
        print(json.dumps(mapping["summary"], indent=2))

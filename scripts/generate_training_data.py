#!/usr/bin/env python3
"""Generate synthetic LoRA fine-tuning data for MedGemma spatial transcriptomics reports."""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

PATTERNS = [
    "hot_immune",
    "cold_immune_desert",
    "stromal_fibrotic",
    "tls_forming",
    "mixed_dual",
    "high_heterogeneous",
    "myeloid_dominant",
]

EXAMPLES_PER_PATTERN = 20


def _total_spots(rng):
    return int(rng.integers(3000, 20001))


def _cell_counts(rng, total, fractions):
    keys = list(fractions.keys())
    raw = np.array([fractions[k] for k in keys], dtype=float)
    raw /= raw.sum()
    counts = (raw * total).astype(int)
    counts[np.argmax(counts)] += total - counts.sum()
    return {k: int(v) for k, v in zip(keys, counts)}


def generate_hot_immune(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.20, 0.30),
        "Stromal_CAF": rng.uniform(0.05, 0.12),
        "T_cells": rng.uniform(0.18, 0.26),
        "B_cells": rng.uniform(0.04, 0.08),
        "NK_cells": rng.uniform(0.07, 0.12),
        "Macrophages": rng.uniform(0.04, 0.08),
        "Dendritic_cells": rng.uniform(0.02, 0.06),
        "Plasma_cells": rng.uniform(0.01, 0.04),
        "Neutrophils": rng.uniform(0.01, 0.03),
        "Endothelial": rng.uniform(0.02, 0.05),
        "Unknown": rng.uniform(0.01, 0.04),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(8, 15)),
            "n_cell_types": int(rng.integers(10, 18)),
            "mean_confidence": float(rng.uniform(0.78, 0.95)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.51, 0.75)),
            "entropy_mean": float(rng.uniform(1.4, 2.5)),
            "n_enriched_pairs": int(rng.integers(100, 201)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.8, 1.6))},
    }


def generate_cold_immune_desert(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.45, 0.60),
        "Stromal_CAF": rng.uniform(0.20, 0.30),
        "T_cells": rng.uniform(0.01, 0.03),
        "B_cells": rng.uniform(0.005, 0.015),
        "NK_cells": rng.uniform(0.003, 0.010),
        "Macrophages": rng.uniform(0.01, 0.025),
        "Dendritic_cells": rng.uniform(0.002, 0.008),
        "Plasma_cells": rng.uniform(0.002, 0.007),
        "Neutrophils": rng.uniform(0.002, 0.007),
        "Endothelial": rng.uniform(0.03, 0.07),
        "Unknown": rng.uniform(0.02, 0.05),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(5, 10)),
            "n_cell_types": int(rng.integers(8, 14)),
            "mean_confidence": float(rng.uniform(0.65, 0.82)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.05, 0.148)),
            "entropy_mean": float(rng.uniform(0.1, 0.55)),
            "n_enriched_pairs": int(rng.integers(0, 25)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.1, 0.5))},
    }


def generate_stromal_fibrotic(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.15, 0.25),
        "Stromal_CAF": rng.uniform(0.40, 0.55),
        "T_cells": rng.uniform(0.04, 0.09),
        "B_cells": rng.uniform(0.01, 0.04),
        "NK_cells": rng.uniform(0.01, 0.03),
        "Macrophages": rng.uniform(0.03, 0.07),
        "Dendritic_cells": rng.uniform(0.01, 0.03),
        "Plasma_cells": rng.uniform(0.01, 0.03),
        "Neutrophils": rng.uniform(0.01, 0.03),
        "Endothelial": rng.uniform(0.02, 0.05),
        "Unknown": rng.uniform(0.01, 0.04),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(6, 12)),
            "n_cell_types": int(rng.integers(9, 16)),
            "mean_confidence": float(rng.uniform(0.70, 0.88)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.20, 0.40)),
            "entropy_mean": float(rng.uniform(0.5, 1.2)),
            "n_enriched_pairs": int(rng.integers(20, 80)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.3, 0.9))},
    }


def generate_tls_forming(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.20, 0.30),
        "Stromal_CAF": rng.uniform(0.08, 0.16),
        "T_cells": rng.uniform(0.12, 0.20),
        "B_cells": rng.uniform(0.10, 0.18),
        "NK_cells": rng.uniform(0.03, 0.07),
        "Macrophages": rng.uniform(0.04, 0.08),
        "Dendritic_cells": rng.uniform(0.05, 0.10),
        "Plasma_cells": rng.uniform(0.04, 0.09),
        "Neutrophils": rng.uniform(0.01, 0.03),
        "Endothelial": rng.uniform(0.02, 0.05),
        "Unknown": rng.uniform(0.01, 0.04),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(7, 13)),
            "n_cell_types": int(rng.integers(10, 18)),
            "mean_confidence": float(rng.uniform(0.75, 0.93)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.42, 0.68)),
            "entropy_mean": float(rng.uniform(1.2, 2.1)),
            "n_enriched_pairs": int(rng.integers(80, 160)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.6, 1.3))},
    }


def generate_mixed_dual(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.26, 0.34),
        "Stromal_CAF": rng.uniform(0.10, 0.18),
        "T_cells": rng.uniform(0.14, 0.20),
        "B_cells": rng.uniform(0.05, 0.10),
        "NK_cells": rng.uniform(0.04, 0.08),
        "Macrophages": rng.uniform(0.05, 0.10),
        "Dendritic_cells": rng.uniform(0.02, 0.05),
        "Plasma_cells": rng.uniform(0.02, 0.05),
        "Neutrophils": rng.uniform(0.01, 0.04),
        "Endothelial": rng.uniform(0.02, 0.05),
        "Unknown": rng.uniform(0.01, 0.04),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(7, 13)),
            "n_cell_types": int(rng.integers(10, 17)),
            "mean_confidence": float(rng.uniform(0.72, 0.90)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.25, 0.50)),
            "entropy_mean": float(rng.uniform(0.9, 1.6)),
            "n_enriched_pairs": int(rng.integers(40, 120)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.5, 1.1))},
    }


def generate_high_heterogeneous(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.15, 0.25),
        "Stromal_CAF": rng.uniform(0.10, 0.20),
        "T_cells": rng.uniform(0.08, 0.15),
        "B_cells": rng.uniform(0.05, 0.12),
        "NK_cells": rng.uniform(0.04, 0.09),
        "Macrophages": rng.uniform(0.07, 0.14),
        "Dendritic_cells": rng.uniform(0.04, 0.09),
        "Plasma_cells": rng.uniform(0.03, 0.07),
        "Neutrophils": rng.uniform(0.03, 0.07),
        "Endothelial": rng.uniform(0.03, 0.07),
        "Unknown": rng.uniform(0.04, 0.10),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(10, 16)),
            "n_cell_types": int(rng.integers(14, 22)),
            "mean_confidence": float(rng.uniform(0.65, 0.78)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.05, 0.199)),
            "entropy_mean": float(rng.uniform(1.51, 2.5)),
            "n_enriched_pairs": int(rng.integers(10, 50)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(1.2, 2.5))},
    }


def generate_myeloid_dominant(rng):
    total = _total_spots(rng)
    fracs = {
        "Epithelial": rng.uniform(0.20, 0.30),
        "Stromal_CAF": rng.uniform(0.08, 0.14),
        "T_cells": rng.uniform(0.05, 0.10),
        "B_cells": rng.uniform(0.02, 0.06),
        "NK_cells": rng.uniform(0.02, 0.05),
        "Macrophages": rng.uniform(0.20, 0.30),
        "Dendritic_cells": rng.uniform(0.03, 0.07),
        "Plasma_cells": rng.uniform(0.01, 0.04),
        "Neutrophils": rng.uniform(0.10, 0.15),
        "Endothelial": rng.uniform(0.02, 0.05),
        "Unknown": rng.uniform(0.02, 0.05),
    }
    return {
        "annotation": {
            "n_spots": total,
            "n_clusters": int(rng.integers(6, 12)),
            "n_cell_types": int(rng.integers(9, 16)),
            "mean_confidence": float(rng.uniform(0.70, 0.88)),
            "cell_type_counts": _cell_counts(rng, total, fracs),
        },
        "spatial_heterogeneity": {
            "morans_i_mean": float(rng.uniform(0.28, 0.55)),
            "entropy_mean": float(rng.uniform(0.8, 1.5)),
            "n_enriched_pairs": int(rng.integers(40, 130)),
        },
        "uncertainty": {"mean_prediction_entropy": float(rng.uniform(0.5, 1.2))},
    }


PATTERN_GENERATORS = {
    "hot_immune": generate_hot_immune,
    "cold_immune_desert": generate_cold_immune_desert,
    "stromal_fibrotic": generate_stromal_fibrotic,
    "tls_forming": generate_tls_forming,
    "mixed_dual": generate_mixed_dual,
    "high_heterogeneous": generate_high_heterogeneous,
    "myeloid_dominant": generate_myeloid_dominant,
}

SYSTEM_PROMPT = (
    "You are an expert computational pathologist specializing in spatial transcriptomics. "
    "Generate clinical pathology reports that identify tissue mechanisms of action, cite "
    "relevant literature, and recommend therapeutic strategies. Never list raw cell counts. "
    "Always synthesize biological meaning."
)


def build_user_prompt(features):
    return (
        "You are analyzing Visium HD spatial transcriptomics data from human tissue.\n\n"
        "SPATIAL FEATURES:\n"
        f"{json.dumps(features, indent=2)}\n\n"
        "Generate a 280-350 word clinical pathology report that:\n"
        "1. Identifies the tissue microenvironment mechanism of action (MoA) based on "
        "cellular spatial patterns\n"
        "2. References 1-2 relevant published findings (cite author, journal, year)\n"
        "3. Recommends a specific therapeutic or mitigation strategy\n"
        "4. Notes measurement uncertainty and suggests validation approach\n\n"
        "Write in formal clinical pathology style. Do NOT list raw cell counts or "
        "enumerate statistics."
    )


def format_chatml(user_prompt, completion):
    text = (
        f"<start_of_turn>user\n{user_prompt}<end_of_turn>\n"
        f"<start_of_turn>model\n{completion}<end_of_turn>"
    )
    return {"text": text}


def call_claude(client, user_prompt, max_retries=3, backoff=5.0):
    """Call Claude API with exponential-backoff retry."""
    for attempt in range(1, max_retries + 1):
        try:
            message = client.messages.create(
                model="claude-sonnet-4-5",
                max_tokens=1024,
                system=SYSTEM_PROMPT,
                messages=[{"role": "user", "content": user_prompt}],
            )
            return message.content[0].text
        except Exception as exc:
            if attempt == max_retries:
                raise RuntimeError(
                    f"Claude API failed after {max_retries} attempts: {exc}"
                ) from exc
            print(
                f"  [retry {attempt}/{max_retries}] API error: {exc}. "
                f"Waiting {backoff}s...",
                file=sys.stderr,
            )
            time.sleep(backoff)
    raise RuntimeError("Unexpected retry loop exit")


def build_example_list(n, seed=42):
    """Build list of (pattern, features) for n examples cycling through all patterns."""
    rng = np.random.default_rng(seed)
    pattern_cycle = (PATTERNS * ((n // len(PATTERNS)) + 1))[:n]
    return [(pat, PATTERN_GENERATORS[pat](rng)) for pat in pattern_cycle]


def compute_split(n):
    """Number of train examples (remainder go to eval, preserving 120/20 ratio)."""
    if n == 140:
        return 120
    return max(1, n - max(1, round(n * 20 / 140)))


def _write_jsonl(path, records):
    """Write list of dicts as newline-delimited JSON."""
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")


def generate_dataset(n, output_dir, dry_run=False, seed=42):
    """Generate n training examples and write train/eval JSONL splits."""
    examples = build_example_list(n, seed=seed)

    if dry_run:
        print(f"\n{'='*60}")
        print(f"DRY RUN -- showing {min(3, n)} examples (no API calls)")
        print(f"{'='*60}\n")
        for i, (pattern, features) in enumerate(examples[:3]):
            print(f"--- Example {i+1} | pattern: {pattern} ---")
            print("Feature JSON:")
            print(json.dumps(features, indent=2))
            print("\nUser prompt (first 400 chars):")
            prompt = build_user_prompt(features)
            print(prompt[:400] + "...")
            target = "train.jsonl" if i < compute_split(n) else "eval.jsonl"
            print(f"\nChatML record would be appended to: {output_dir}/{target}\n")
        train_n = compute_split(n)
        print(f"Full run would produce {train_n} train + {n - train_n} eval examples.")
        return

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("ERROR: ANTHROPIC_API_KEY environment variable is not set.", file=sys.stderr)
        sys.exit(1)

    try:
        import anthropic
    except ImportError:
        print(
            "ERROR: 'anthropic' package not installed. Run: pip install anthropic",
            file=sys.stderr,
        )
        sys.exit(1)

    client = anthropic.Anthropic(api_key=api_key)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / ".gitkeep").touch()

    train_n = compute_split(n)
    train_path = output_dir / "train.jsonl"
    eval_path = output_dir / "eval.jsonl"
    train_records = []
    eval_records = []

    for idx, (pattern, features) in enumerate(examples):
        split = "train" if idx < train_n else "eval"
        example_num = (idx % EXAMPLES_PER_PATTERN) + 1
        print(f"[{idx+1}/{n}] {pattern} example {example_num} -> {split}...")

        user_prompt = build_user_prompt(features)
        completion = call_claude(client, user_prompt)
        record = format_chatml(user_prompt, completion)

        if split == "train":
            train_records.append(record)
        else:
            eval_records.append(record)

    _write_jsonl(train_path, train_records)
    _write_jsonl(eval_path, eval_records)

    print(f"\nDone.")
    print(f"  Train: {len(train_records)} examples -> {train_path}")
    print(f"  Eval:  {len(eval_records)} examples -> {eval_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate synthetic LoRA training data for MedGemma spatial reports."
    )
    parser.add_argument(
        "--n",
        type=int,
        default=140,
        help="Total number of examples to generate (default: 140).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print 3 example feature JSONs and prompts without calling the API.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/lora_training"),
        help="Directory to write train.jsonl and eval.jsonl (default: data/lora_training/).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.n < 1:
        print("ERROR: --n must be >= 1.", file=sys.stderr)
        sys.exit(1)

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        repo_root = Path(__file__).parent.parent
        output_dir = repo_root / output_dir

    generate_dataset(
        n=args.n,
        output_dir=output_dir,
        dry_run=args.dry_run,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()

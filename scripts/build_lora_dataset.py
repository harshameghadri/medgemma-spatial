#!/usr/bin/env python3
"""Re-render the existing LoRA data with the serving prompt (no API calls).

The original data (``generate_training_data.py``) used a raw JSON dump as the
user turn, while the app prompts MedGemma with
``prompt_builder.generate_medgemma_prompt``. That train/serve skew means the
adapter never sees the prompt it is used with. The features are deterministic
(``build_example_list(n, seed)``), so they can be regenerated, verified against
the stored prompt (with float tolerance), and re-rendered in TRL prompt-completion format:

    {"prompt": [{"role": "system", ...}, {"role": "user", ...}],
     "completion": [{"role": "assistant", ...}],
     "pattern": "..."}

Usage:
    python scripts/build_lora_dataset.py --in-dir data/lora_training \
        --out-dir data/lora_training_v2
"""

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from generate_training_data import build_example_list, build_user_prompt  # noqa: E402
from src.report_generation.medgemma_multimodal import SYSTEM_PROMPT  # noqa: E402
from src.report_generation.prompt_builder import generate_medgemma_prompt  # noqa: E402

USER_TAG = "<start_of_turn>user\n"
MODEL_TAG = "<start_of_turn>model\n"
END_TAG = "<end_of_turn>"
FEATURES_START = "SPATIAL FEATURES:\n"
FEATURES_END = "\n\nGenerate a"


def split_chatml(text: str):
    """Return (user_prompt, completion) from a stored Gemma-format record."""
    user = text.split(USER_TAG, 1)[1].split(END_TAG, 1)[0]
    completion = text.split(MODEL_TAG, 1)[1]
    if completion.endswith(END_TAG):
        completion = completion[: -len(END_TAG)]
    return user, completion.strip()


def to_prompt_completion(features: dict, completion: str, pattern: str) -> dict:
    """Serving-identical prompt plus the stored completion."""
    return {
        "prompt": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": generate_medgemma_prompt(features)},
        ],
        "completion": [{"role": "assistant", "content": completion}],
        "pattern": pattern,
    }


def extract_features(user_prompt: str) -> dict:
    """Parse the feature JSON embedded in an original training prompt."""
    try:
        body = user_prompt.split(FEATURES_START, 1)[1].split(FEATURES_END, 1)[0]
    except IndexError as e:
        raise ValueError("prompt does not contain a SPATIAL FEATURES block") from e
    return json.loads(body)


def features_match(a, b, rel_tol: float = 1e-9) -> bool:
    """Structural equality with float tolerance (numpy RNG differs in the last ULP across builds)."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(features_match(a[k], b[k], rel_tol) for k in a)
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=rel_tol, abs_tol=1e-12)
    return a == b


def rebuild(records: list, seed: int = 42) -> list:
    """Parse each stored prompt's features, cross-check them against the seeded
    generator (which also recovers the pattern label), and re-render."""
    examples = build_example_list(len(records), seed=seed)
    out = []
    for i, ((pattern, regenerated), record) in enumerate(zip(examples, records)):
        user, completion = split_chatml(record["text"])
        stored = extract_features(user)
        if not features_match(stored, regenerated):
            raise ValueError(
                f"Record {i}: stored features do not match build_example_list(seed={seed}); "
                "the data was not produced by the seeded generator"
            )
        if user.replace(json.dumps(stored, indent=2), "") != build_user_prompt({}).replace("{}", ""):
            raise ValueError(f"Record {i}: prompt template differs from build_user_prompt")
        out.append(to_prompt_completion(stored, completion, pattern))
    return out


def read_jsonl(path: Path) -> list:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_jsonl(path: Path, rows: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in-dir", type=Path, default=REPO / "data/lora_training")
    parser.add_argument("--out-dir", type=Path, default=REPO / "data/lora_training_v2")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    train = read_jsonl(args.in_dir / "train.jsonl")
    evals = read_jsonl(args.in_dir / "eval.jsonl")
    rows = rebuild(train + evals, seed=args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(args.out_dir / "train.jsonl", rows[: len(train)])
    write_jsonl(args.out_dir / "eval.jsonl", rows[len(train):])
    print(f"Verified {len(rows)} records against seed {args.seed}; wrote "
          f"{len(train)} train / {len(evals)} eval to {args.out_dir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Rewrite an eval set with grounded prompts; keeps the original prompt for scoring.

Each row gains ``source_prompt`` (the original user message, which grounding
checks score against) and ``added`` (the inserted evidence and examples, for the
separate copying metric).

    python scripts/build_grounded_eval.py --eval data/lora_v4/eval.jsonl --out data/grounded_v1/eval.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.report_generation.grounding import ground_prompt  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--eval", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-evidence", action="store_true")
    ap.add_argument("--no-examples", action="store_true")
    args = ap.parse_args(argv)
    out = []
    for line in open(args.eval):
        row = json.loads(line)
        source = row["prompt"][-1]["content"]
        grounded, added = ground_prompt(source, evidence=not args.no_evidence, examples=not args.no_examples)
        row["prompt"] = row["prompt"][:-1] + [{"role": "user", "content": grounded}]
        out.append({**row, "source_prompt": source, "added": added})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("".join(json.dumps(r) + "\n" for r in out))
    print(f"{len(out)} rows -> {args.out}; {sum(bool(r['added']) for r in out)} grounded")


if __name__ == "__main__":
    main()

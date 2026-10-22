#!/usr/bin/env python3
"""Build the retrieval corpus from reviewer-picked training reports.

Keeps each picked draft (``build_distill_set.py`` output) that the biology check
does not flag, with the parsed data summary used for retrieval. Training rows
only: the eval split is never read.

    python scripts/build_approved_corpus.py --train data/distill_v1/train.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.report_generation.biology_check import check_report  # noqa: E402
from src.report_generation.grounding import APPROVED, parse_summary  # noqa: E402


def build(rows):
    """(kept corpus rows, ids dropped by the biology check)."""
    kept, dropped = [], []
    for r in rows:
        prompt, text = r["prompt"][-1]["content"], r["completion"][0]["content"]
        if check_report(text, prompt, r.get("pattern")):
            dropped.append(r["source"])
            continue
        kept.append({"id": r["source"], "pattern": r.get("pattern"), "summary": parse_summary(prompt), "text": text})
    return kept, dropped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--train", required=True)
    ap.add_argument("--out", default=str(APPROVED))
    args = ap.parse_args(argv)
    kept, dropped = build([json.loads(line) for line in open(args.train)])
    Path(args.out).write_text("".join(json.dumps(r) + "\n" for r in kept))
    print(f"{len(kept)} kept, {len(dropped)} dropped by the biology check: {', '.join(dropped)}")


if __name__ == "__main__":
    main()

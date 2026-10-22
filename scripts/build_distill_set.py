#!/usr/bin/env python3
"""Turn reviewed teacher drafts into a prompt-completion training set.

``decisions.json`` maps a card id (``r007``) to ``{"pick": 0 | 1 | "none", "note": ...}``,
as exported from the review page. A pick indexes that prompt's passing drafts in
variant order, the same order the review cards showed. Only picked drafts are kept,
and a picked draft must still pass the automatic gates.

    python scripts/build_distill_set.py --drafts data/distill_v1/drafts.jsonl \\
        --decisions data/distill_v1/decisions.json --train data/lora_v4/train.jsonl \\
        --out data/distill_v1/train.jsonl
"""

import argparse
import collections
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from draft_targets import gates  # noqa: E402


def passing_by_row(drafts):
    """Passing drafts grouped by training row, in variant order."""
    rows = collections.defaultdict(list)
    for d in sorted(drafts, key=lambda d: (d["row"], d["variant"])):
        if d["passes"]:
            rows[d["row"]].append(d)
    return rows


def build(drafts, decisions, train_rows):
    """Picked (prompt, completion) rows plus a per-pattern tally."""
    passing = passing_by_row(drafts)
    out, kept = [], collections.Counter()
    for card, dec in sorted(decisions.items()):
        if dec["pick"] == "none":
            continue
        row = int(card.lstrip("r"))
        d = passing[row][dec["pick"]]
        prompt = train_rows[row]["prompt"]
        if not gates(d["text"], prompt[-1]["content"])["passes"]:
            raise ValueError(f"{card}: picked draft no longer passes the gates")
        out.append({"prompt": prompt, "completion": [{"role": "assistant", "content": d["text"]}],
                    "pattern": train_rows[row].get("pattern"), "source": f"{card}.v{d['variant']}"})
        kept[train_rows[row].get("pattern")] += 1
    return out, kept


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drafts", required=True)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--train", required=True, help="the train.jsonl the drafts were made from")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    drafts = [json.loads(line) for line in open(args.drafts)]
    decisions = json.loads(Path(args.decisions).read_text())
    train_rows = [json.loads(line) for line in open(args.train)]
    out, kept = build(drafts, decisions, train_rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.writelines(json.dumps(r) + "\n" for r in out)
    rejected = sum(d["pick"] == "none" for d in decisions.values())
    print(f"{len(out)} picked, {rejected} rejected -> {args.out}")
    print(dict(sorted(kept.items())))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Draft LoRA targets with a local teacher model, then gate them with the grounding checks.

Only the train split is drafted; eval prompts are never shown to the teacher.
Each prompt gets one greedy draft and ``--samples`` sampled drafts. A draft
passes the automatic gates when it has no citation, no number or cell type
absent from the prompt, all four sections (after heading normalisation), a
length within REPORT_WORDS (+20) and no 5+ word run copied from the prompt.
Passing drafts still need human review (``scripts/review_drafts.py``) before
they become training data.

    python scripts/draft_targets.py --teacher unsloth/Meta-Llama-3.1-8B-Instruct-bnb-4bit \\
        --train data/lora_v4/train.jsonl --out data/distill_v1/drafts.jsonl
"""

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from compare_adapter import load_model  # noqa: E402

from src.report_generation.faithfulness import normalise_headings, score_report, strip_thinking  # noqa: E402
from src.report_generation.prompt_builder import REPORT_WORDS  # noqa: E402


def gates(text: str, prompt: str) -> dict:
    """Automatic pass/fail checks for one draft."""
    s = score_report(text, prompt)
    lo, hi = REPORT_WORDS
    checks = {
        "no_citation": s["citation_sentences"] == 0,
        "numbers_grounded": not s["unsupported_numbers"],
        "cell_types_grounded": not s["unsupported_cell_types"],
        "four_sections": len(s["sections"]) == 4,
        "length": lo <= s["words"] <= hi + 20,
        "no_prompt_copy": not s["copied_spans"],
    }
    return {"checks": checks, "passes": all(checks.values()), "score": s}


def draft(tok, model, messages, sample: bool, seed: int, max_new_tokens: int) -> str:
    """One teacher draft, thinking stripped and headings normalised."""
    import torch

    ids = tok.apply_chat_template(messages, add_generation_prompt=True,
                                  return_tensors="pt", return_dict=True).to(model.device)
    kwargs = {"do_sample": False}
    if sample:
        torch.manual_seed(seed)
        kwargs = {"do_sample": True, "temperature": 0.7, "top_p": 0.9}
    with torch.inference_mode():
        out = model.generate(**ids, max_new_tokens=max_new_tokens, **kwargs)
    text = tok.decode(out[0][ids["input_ids"].shape[-1]:], skip_special_tokens=True)
    return normalise_headings(strip_thinking(text))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--teacher", required=True)
    ap.add_argument("--train", required=True, help="train.jsonl (prompt-completion rows)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--samples", type=int, default=1, help="sampled drafts per prompt, besides the greedy one")
    ap.add_argument("--max-new-tokens", type=int, default=450)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    if "eval" in Path(args.train).name:
        sys.exit("refusing to draft from an eval file")
    rows = [json.loads(line) for line in open(args.train)]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        done = {(d["row"], d["variant"]) for d in map(json.loads, out.open())}
    tok, model = load_model(args.teacher, None)
    with out.open("a") as fh:
        for i, row in enumerate(rows):
            prompt = row["prompt"][-1]["content"]
            for v in range(args.samples + 1):
                if (i, v) in done:
                    continue
                t = time.time()
                text = draft(tok, model, row["prompt"], v > 0, args.seed + 1000 * i + v, args.max_new_tokens)
                g = gates(text, prompt)
                rec = {"row": i, "variant": v, "pattern": row.get("pattern"), "teacher": args.teacher,
                       "seconds": round(time.time() - t, 1), "text": text, **g}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                failed = [k for k, ok in g["checks"].items() if not ok]
                print(f"[{i}.{v}] {str(row.get('pattern')):<20} {'PASS' if g['passes'] else 'fail ' + ','.join(failed)}",
                      flush=True)


if __name__ == "__main__":
    main()

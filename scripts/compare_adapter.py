#!/usr/bin/env python3
"""Compare base MedGemma vs a LoRA adapter on held-out prompts.

For each eval example, generates a report with the adapter disabled and
enabled, then records length, section coverage, MoA term hits, citations and
**unsupported numbers** (numeric claims absent from the prompt) - the main
failure mode of adapter v1 (see docs/REVIEW.md, L9).

    python scripts/compare_adapter.py --adapter runs/adapter_v1 \
        --eval data/lora_training_v2/eval.jsonl --n 6
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

CITATION = re.compile(r"[A-Z][a-zA-Z\-]+ et al\.?,? \(?[^)]{0,60}?(?:19|20)\d\d")
NUMBER = re.compile(r"\d+(?:\.\d+)?%?")
SECTIONS = ("MICROENVIRONMENT", "MECHANISM", "THERAPEUTIC", "VALIDATION")


def unsupported_numbers(text: str, prompt: str) -> list:
    """Numbers in the report that do not appear in the prompt (ignores 1-5, years, CD markers)."""
    allowed = set(NUMBER.findall(prompt))
    text = re.sub(r"\bCD\d+[a-z]?\b", "", text)
    return [n for n in NUMBER.findall(text)
            if n not in allowed and not re.fullmatch(r"[1-5]|(19|20)\d\d", n)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model-id", default="google/medgemma-4b-it")
    p.add_argument("--adapter", required=True)
    p.add_argument("--eval", required=True)
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--max-new-tokens", type=int, default=500)
    args = p.parse_args(argv)

    import torch
    from peft import PeftModel
    from transformers import AutoModelForImageTextToText, AutoTokenizer, BitsAndBytesConfig

    from src.report_generation.prompt_builder import evaluate_report_quality

    rows = [json.loads(line) for line in open(args.eval)][: args.n]
    tok = AutoTokenizer.from_pretrained(args.model_id)
    base = AutoModelForImageTextToText.from_pretrained(
        args.model_id, device_map={"": 0}, dtype=torch.bfloat16,
        quantization_config=BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_quant_type="nf4",
                                               bnb_4bit_compute_dtype=torch.bfloat16))
    model = PeftModel.from_pretrained(base, args.adapter).eval()

    def generate(row, use_adapter):
        ids = tok.apply_chat_template(row["prompt"], add_generation_prompt=True,
                                      return_tensors="pt", return_dict=True).to(model.device)
        with torch.inference_mode():
            if use_adapter:
                out = model.generate(**ids, max_new_tokens=args.max_new_tokens, do_sample=False)
            else:
                with model.disable_adapter():
                    out = model.generate(**ids, max_new_tokens=args.max_new_tokens, do_sample=False)
        return tok.decode(out[0][ids["input_ids"].shape[-1]:], skip_special_tokens=True).strip()

    results, totals = [], {"base": [0, 0], "adapter": [0, 0]}
    for i, row in enumerate(rows):
        prompt = row["prompt"][-1]["content"]
        rec = {"pattern": row.get("pattern"), "reference": row["completion"][0]["content"]}
        for name, flag in (("base", False), ("adapter", True)):
            t = time.time()
            text = generate(row, flag)
            bad = unsupported_numbers(text, prompt)
            cites = CITATION.findall(text)
            rec[name] = {
                "text": text, "seconds": round(time.time() - t, 1), "words": len(text.split()),
                "sections": sum(h in text.upper() for h in SECTIONS),
                "moa_hits": evaluate_report_quality(text, {})["moa_hits"],
                "citations": cites, "unsupported_numbers": bad,
            }
            totals[name][0] += len(bad)
            totals[name][1] += len(cites)
            print(f"[{i}] {str(row.get('pattern')):<20} {name:<7} words={rec[name]['words']:<4} "
                  f"moa={rec[name]['moa_hits']} cites={len(cites)} unsupported={bad}", flush=True)
        results.append(rec)

    out = Path(args.adapter) / "comparison.json"
    out.write_text(json.dumps({"totals": totals, "examples": results}, indent=2))
    print(f"TOTAL unsupported numbers / citations: {totals}\nsaved {out}")


if __name__ == "__main__":
    main()

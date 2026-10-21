#!/usr/bin/env python3
"""Compare base MedGemma vs a LoRA adapter on held-out prompts.

For each eval example, generates a report with the adapter disabled and
enabled and scores grounding against the prompt (``faithfulness.score_report``):
length, sections, MoA terms, citations, unsupported numbers and unsupported
cell types. ``--pubmed`` also checks every citation against PubMed.

    python scripts/compare_adapter.py --adapter runs/adapter_v2 --eval data/lora_v3/eval.jsonl
    python scripts/compare_adapter.py --eval data/lora_v3/eval.jsonl --out runs/base.json  # base only
    python scripts/compare_adapter.py --rescore runs/adapter_v1/comparison.json --eval ...  # no GPU
"""

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.report_generation.faithfulness import WORD, HEADING, score_report  # noqa: E402
from src.report_generation.prompt_builder import REPORT_WORDS, evaluate_report_quality  # noqa: E402

PASS_CRITERIA = (
    "citations: 0 sentences with a citation marker",
    "unsupported numbers: total <= base",
    "unsupported cell types: total <= base",
    "length: share in REPORT_WORDS range > base",
    "sections: share with all 4 >= base",
    "MoA terms: mean >= base - 1",
)
REPORTED = (
    "moa_novel_mean: MoA terms outside prompt-copied spans",
    "reports_copying_prompt: reports with a 5+ word run copied from the prompt",
)


def score(text: str, prompt: str, pubmed=None) -> dict:
    rec = score_report(text, prompt)
    rec["moa_hits"] = evaluate_report_quality(text, {})["moa_hits"]
    rec["moa_novel"] = evaluate_report_quality(own_words(text, rec["copied_spans"]), {})["moa_hits"]
    if pubmed is not None:
        from src.report_generation.pubmed_check import check_citations
        rec["pubmed"] = [{k: c.get(k) for k in ("text", "status", "pmid")} for c in check_citations(text, pubmed)]
    return rec


def own_words(text: str, spans: list) -> str:
    """Report body (headings dropped, lower-cased words) with prompt-copied spans cut out."""
    body = " ".join(w.lower() for line in text.split("\n") if not HEADING.match(line) for w in WORD.findall(line))
    for span in sorted(spans, key=len, reverse=True):
        body = body.replace(span, " | ")
    return body


def summarise(results: list, name: str) -> dict:
    recs = [r[name] for r in results if name in r]
    lo, hi = REPORT_WORDS
    return {
        "n": len(recs),
        "words_median": statistics.median(r["words"] for r in recs),
        "in_word_range": sum(lo <= r["words"] <= hi + 20 for r in recs),
        "all_sections": sum(len(r["sections"]) == 4 for r in recs),
        "moa_mean": round(statistics.mean(r["moa_hits"] for r in recs), 2),
        "moa_novel_mean": round(statistics.mean(r["moa_novel"] for r in recs), 2),
        "reports_copying_prompt": sum(bool(r["copied_spans"]) for r in recs),
        "copied_words": sum(len(x.split()) for r in recs for x in r["copied_spans"]),
        "reports_with_citation": sum(r["citation_sentences"] > 0 for r in recs),
        "unsupported_numbers": sum(len(r["unsupported_numbers"]) for r in recs),
        "reports_with_unsupported_number": sum(bool(r["unsupported_numbers"]) for r in recs),
        "unsupported_cell_types": sum(len(r["unsupported_cell_types"]) for r in recs),
        "pubmed_verified": sum(c["status"] == "verified" for r in recs for c in r.get("pubmed", [])),
        "pubmed_unverified": sum(c["status"] != "verified" for r in recs for c in r.get("pubmed", [])),
    }


def verdict(summary: dict) -> dict:
    """Apply PASS_CRITERIA (fixed before training) to an adapter-vs-base summary."""
    a, b = summary.get("adapter"), summary["base"]
    if a is None:
        return {}
    checks = {
        "citations": a["reports_with_citation"] == 0,
        "unsupported_numbers": a["unsupported_numbers"] <= b["unsupported_numbers"],
        "unsupported_cell_types": a["unsupported_cell_types"] <= b["unsupported_cell_types"],
        "length": a["in_word_range"] > b["in_word_range"],
        "sections": a["all_sections"] >= b["all_sections"],
        "moa": a["moa_mean"] >= b["moa_mean"] - 1,
    }
    return {**checks, "pass": all(checks.values())}


def load_model(model_id: str, adapter):
    """4-bit model on GPU 0: Gemma 3 vision-language checkpoints or any text-only causal LM."""
    import torch
    from transformers import (AutoConfig, AutoModelForCausalLM, AutoModelForImageTextToText, AutoTokenizer,
                              BitsAndBytesConfig)

    tok = AutoTokenizer.from_pretrained(model_id)
    config = AutoConfig.from_pretrained(model_id)
    cls = AutoModelForImageTextToText if hasattr(config, "vision_config") else AutoModelForCausalLM
    kwargs = {"device_map": {"": 0}, "dtype": torch.bfloat16}
    if not getattr(config, "quantization_config", None):
        kwargs["quantization_config"] = BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_quant_type="nf4",
                                                           bnb_4bit_compute_dtype=torch.bfloat16)
    model = cls.from_pretrained(model_id, **kwargs)
    if adapter:
        from peft import PeftModel
        model = PeftModel.from_pretrained(model, adapter)
    return tok, model.eval()


def generate(tok, model, row, use_adapter: bool, max_new_tokens: int) -> str:
    import contextlib

    import torch

    ids = tok.apply_chat_template(row["prompt"], add_generation_prompt=True,
                                  return_tensors="pt", return_dict=True).to(model.device)
    ctx = contextlib.nullcontext() if use_adapter or not hasattr(model, "disable_adapter") \
        else model.disable_adapter()
    with torch.inference_mode(), ctx:
        out = model.generate(**ids, max_new_tokens=max_new_tokens, do_sample=False)
    return tok.decode(out[0][ids["input_ids"].shape[-1]:], skip_special_tokens=True).strip()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model-id", default="google/medgemma-4b-it")
    p.add_argument("--adapter", default=None, help="omit to score the base model only")
    p.add_argument("--eval", required=True)
    p.add_argument("--n", type=int, default=0, help="0 = all eval rows")
    p.add_argument("--max-new-tokens", type=int, default=450)
    p.add_argument("--out", default=None)
    p.add_argument("--rescore", default=None, help="re-score texts in an existing comparison.json")
    p.add_argument("--pubmed", action="store_true", help="verify citations against PubMed")
    args = p.parse_args(argv)

    rows = [json.loads(line) for line in open(args.eval)]
    rows = rows[: args.n] if args.n else rows
    pubmed = None
    if args.pubmed:
        from src.report_generation.pubmed_check import PubMedChecker
        pubmed = PubMedChecker(str(Path(args.eval).with_name("pubmed_cache.json")))

    if args.rescore:
        old = json.loads(Path(args.rescore).read_text())
        old = old.get("examples", old)
        texts = [{k: v["text"] for k, v in rec.items() if isinstance(v, dict) and "text" in v} for rec in old]
        prompts = [r["prompt"][-1]["content"] for r in rows]
        if len(texts) != len(prompts):
            sys.exit(f"--rescore has {len(texts)} examples but --eval has {len(prompts)}")
        gen = None
    else:
        tok, model = load_model(args.model_id, args.adapter)
        names = ("base", "adapter") if args.adapter else ("base",)

        def gen(row):
            out = {}
            for name in names:
                t = time.time()
                out[name] = (generate(tok, model, row, name == "adapter", args.max_new_tokens),
                             round(time.time() - t, 1))
            return out

    results = []
    for i, row in enumerate(rows):
        prompt = row["prompt"][-1]["content"]
        rec = {"pattern": row.get("pattern"), "reference": row["completion"][0]["content"]}
        generated = {k: (v, None) for k, v in texts[i].items()} if gen is None else gen(row)
        for name, (text, secs) in generated.items():
            rec[name] = {"text": text, "seconds": secs, **score(text, prompt, pubmed)}
            r = rec[name]
            print(f"[{i}] {str(row.get('pattern')):<20} {name:<7} words={r['words']:<4} "
                  f"sections={len(r['sections'])} moa={r['moa_hits']} cites={r['citation_sentences']} "
                  f"num={r['unsupported_numbers']} ct={r['unsupported_cell_types']}", flush=True)
        results.append(rec)

    summary = {name: summarise(results, name) for name in ("base", "adapter") if name in results[0]}
    summary["verdict"] = verdict(summary)
    out = Path(args.out or (Path(args.adapter) / "comparison.json" if args.adapter else "comparison.json"))
    out.write_text(json.dumps({"summary": summary, "criteria": PASS_CRITERIA, "reported": REPORTED, "examples": results}, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"saved {out}")


if __name__ == "__main__":
    main()

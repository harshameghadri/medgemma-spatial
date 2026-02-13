#!/usr/bin/env python3
"""Benchmark: base MedGemma 4B-it vs fine-tuned LoRA adapter on spatial pathology reports."""

import argparse
import json
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.report_generation.prompt_builder import generate_medgemma_prompt, evaluate_report_quality


SYNTHETIC_SAMPLES = [
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 1420, "T cells": 380, "Macrophages": 210,
                "Fibroblasts": 150, "Endothelial cells": 90,
            },
            "n_spots": 2250,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.72},
        "uncertainty": {"mean_prediction_entropy": 1.2},
        "pattern": "immune_hot_tumor",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 1850, "Fibroblasts": 480, "Endothelial cells": 120,
                "T cells": 80, "Macrophages": 70,
            },
            "n_spots": 2600,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.61},
        "uncertainty": {"mean_prediction_entropy": 0.9},
        "pattern": "immune_excluded",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "B cells": 920, "T cells": 760, "Macrophages": 410,
                "Plasma cells": 290, "Dendritic cells": 180, "NK cells": 140,
            },
            "n_spots": 2700,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.55},
        "uncertainty": {"mean_prediction_entropy": 1.8},
        "pattern": "tertiary_lymphoid_structure",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 2100, "Fibroblasts": 630, "Myofibroblasts": 310,
                "Endothelial cells": 160,
            },
            "n_spots": 3200,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.82},
        "uncertainty": {"mean_prediction_entropy": 0.7},
        "pattern": "desmoplastic_stroma",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 950, "Necrotic cells": 820, "Macrophages": 490,
                "Neutrophils": 260, "T cells": 180,
            },
            "n_spots": 2700,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.44},
        "uncertainty": {"mean_prediction_entropy": 2.1},
        "pattern": "necrotic_core",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Epithelial cells": 1680, "Fibroblasts": 520, "Smooth muscle cells": 310,
                "Endothelial cells": 190, "T cells": 140,
            },
            "n_spots": 2840,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.67},
        "uncertainty": {"mean_prediction_entropy": 1.0},
        "pattern": "normal_adjacent_tissue",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 1100, "T cells": 620, "Macrophages": 580,
                "NK cells": 310, "B cells": 230, "Dendritic cells": 160,
            },
            "n_spots": 3000,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.38},
        "uncertainty": {"mean_prediction_entropy": 1.9},
        "pattern": "inflamed_tumor_border",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 2400, "Endothelial cells": 390, "Pericytes": 210,
                "Fibroblasts": 180, "T cells": 120,
            },
            "n_spots": 3300,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.76},
        "uncertainty": {"mean_prediction_entropy": 0.8},
        "pattern": "angiogenic_front",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Macrophages": 1050, "Tumor cells": 780, "T cells": 560,
                "Fibroblasts": 340, "NK cells": 170,
            },
            "n_spots": 2900,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.49},
        "uncertainty": {"mean_prediction_entropy": 1.6},
        "pattern": "myeloid_dominated_tme",
    },
    {
        "annotation": {
            "cell_type_counts": {
                "Tumor cells": 1730, "T cells": 490, "B cells": 320,
                "Plasma cells": 280, "Macrophages": 230, "Fibroblasts": 150,
            },
            "n_spots": 3200,
        },
        "spatial_heterogeneity": {"morans_i_mean": 0.58},
        "uncertainty": {"mean_prediction_entropy": 1.3},
        "pattern": "mixed_immune_infiltrate",
    },
]

MOA_TERMS = ["pathway", "signaling", "mechanism", "via", "mediated by", "driven by", "-dependent", "axis"]
RESEARCH_TERMS = ["et al", "consistent with", "described by", "shown by", "reported by"]
MITIGATION_TERMS = ["treatment", "therapeutic", "inhibit", "blockade", "agent", "inhibitor", "checkpoint", "immunotherapy", "surgery", "radiation"]
_AUTHOR_YEAR_RE = re.compile(r"\b[A-Z][a-z]+ \d{4}\b")


def _count_term_matches(text, terms):
    """Count how many terms appear in text (case-insensitive)."""
    text_lower = text.lower()
    return sum(1 for t in terms if t.lower() in text_lower)


def score_report(report, features):
    """Compute all 5 benchmark metrics for a single report."""
    moa_count = _count_term_matches(report, MOA_TERMS)
    research_count = (_count_term_matches(report, RESEARCH_TERMS) + len(_AUTHOR_YEAR_RE.findall(report)))
    mitigation_count = _count_term_matches(report, MITIGATION_TERMS)
    quality = evaluate_report_quality(report, features)
    wc = len(report.split())
    return {
        "moa_score": {"count": moa_count, "pass": moa_count >= 3},
        "research_score": {"count": research_count, "pass": research_count >= 1},
        "mitigation_score": {"count": mitigation_count, "pass": mitigation_count >= 2},
        "anti_parroting": {"parroting_risk": quality["parroting_risk"], "pass": quality["parroting_risk"] == "LOW"},
        "word_count": {"count": wc, "pass": 250 <= wc <= 400},
        "metrics_passed": sum([moa_count >= 3, research_count >= 1, mitigation_count >= 2, quality["parroting_risk"] == "LOW", 250 <= wc <= 400]),
    }


def _detect_device():
    """Return the best available torch device string."""
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda"
        if torch.backends.mps.is_available():
            return "mps"
    except ImportError:
        pass
    return "cpu"


def load_base_model(model_id, hf_token):
    """Load base MedGemma with 4-bit quantization on CUDA, float32 on CPU/MPS."""
    from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig
    import torch
    device = _detect_device()
    print(f"  Device: {device}")
    tokenizer = AutoTokenizer.from_pretrained(model_id, token=hf_token)
    if device == "cuda":
        quant_cfg = BitsAndBytesConfig(load_in_4bit=True, bnb_4bit_compute_dtype=torch.bfloat16, bnb_4bit_use_double_quant=True, bnb_4bit_quant_type="nf4")
        model = AutoModelForCausalLM.from_pretrained(model_id, quantization_config=quant_cfg, device_map="auto", token=hf_token)
    else:
        model = AutoModelForCausalLM.from_pretrained(model_id, torch_dtype=torch.float32, device_map="auto", low_cpu_mem_usage=True, token=hf_token)
    return tokenizer, model


def load_finetuned_model(base_model, base_tokenizer, adapter_id, hf_token):
    """Wrap base model with LoRA adapter via PEFT."""
    from peft import PeftModel
    model_ft = PeftModel.from_pretrained(base_model, adapter_id, token=hf_token)
    return base_tokenizer, model_ft


def _build_chat_input(tokenizer, prompt):
    """Format prompt using tokenizer chat template, fallback to Gemma manual format."""
    messages = [{"role": "user", "content": prompt}]
    try:
        return tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
    except Exception:
        return f"<start_of_turn>user\n{prompt}<end_of_turn>\n<start_of_turn>model\n"


def generate_report(tokenizer, model, features):
    """Generate a report; return (report_text, generation_time_seconds)."""
    import torch
    prompt = generate_medgemma_prompt(features)
    formatted = _build_chat_input(tokenizer, prompt)
    device = next(model.parameters()).device
    inputs = tokenizer(formatted, return_tensors="pt").to(device)
    input_len = inputs["input_ids"].shape[1]
    t0 = time.time()
    with torch.no_grad():
        output_ids = model.generate(**inputs, max_new_tokens=350, temperature=0.7, top_p=0.9, do_sample=True, pad_token_id=tokenizer.eos_token_id)
    elapsed = time.time() - t0
    new_tokens = output_ids[0][input_len:]
    report_text = tokenizer.decode(new_tokens, skip_special_tokens=True).strip()
    return report_text, elapsed


def evaluate_model(label, tokenizer, model, samples):
    """Run all samples through a model; return per_sample list + aggregate dict."""
    per_sample = []
    for i, features in enumerate(samples, 1):
        pattern = features.get("pattern", f"sample_{i}")
        print(f"    Sample {i:2d}/{len(samples)} (pattern: {pattern}) ...", end=" ", flush=True)
        try:
            report, gen_time = generate_report(tokenizer, model, features)
            scores = score_report(report, features)
            entry = {"sample_index": i, "pattern": pattern, "report_preview": report[:200], "scores": scores, "generation_time_s": round(gen_time, 2)}
            print(f"{scores['metrics_passed']}/5 passed  ({gen_time:.1f}s)")
        except Exception as e:
            entry = {"sample_index": i, "pattern": pattern, "error": str(e), "scores": None, "generation_time_s": 0}
            print(f"ERROR: {e}")
        per_sample.append(entry)
    valid = [s for s in per_sample if s["scores"] is not None]
    n = max(len(valid), 1)
    def pass_rate(metric):
        return round(sum(1 for s in valid if s["scores"][metric]["pass"]) / n, 3)
    aggregate = {
        "moa_score_pass_rate": pass_rate("moa_score"),
        "research_score_pass_rate": pass_rate("research_score"),
        "mitigation_score_pass_rate": pass_rate("mitigation_score"),
        "anti_parroting_pass_rate": pass_rate("anti_parroting"),
        "word_count_pass_rate": pass_rate("word_count"),
        "mean_generation_time_s": round(sum(s["generation_time_s"] for s in per_sample) / len(per_sample), 2),
    }
    aggregate["overall_pass_rate"] = round(sum([aggregate["moa_score_pass_rate"], aggregate["research_score_pass_rate"], aggregate["mitigation_score_pass_rate"], aggregate["anti_parroting_pass_rate"], aggregate["word_count_pass_rate"]]) / 5, 3)
    return {"per_sample": per_sample, "aggregate": aggregate}


def compute_go_decision(base_agg, ft_agg):
    """Return (go_bool, reason_str, winner_per_metric_dict)."""
    metric_map = {"moa_score": "moa_score_pass_rate", "research_score": "research_score_pass_rate", "mitigation_score": "mitigation_score_pass_rate", "anti_parroting": "anti_parroting_pass_rate", "word_count": "word_count_pass_rate"}
    winners = {}
    ft_wins = []
    for label, key in metric_map.items():
        b, f = base_agg[key], ft_agg[key]
        if f > b:
            winners[label] = "finetuned"; ft_wins.append(label)
        elif b > f:
            winners[label] = "base"
        else:
            winners[label] = "tie"
    go = len(ft_wins) >= 3
    reason = (f"Fine-tuned wins on {len(ft_wins)}/5 metrics ({', '.join(ft_wins)})" if go else f"Fine-tuned only wins {len(ft_wins)}/5 metrics")
    return go, reason, winners


METRIC_LABELS = [
    ("moa_score", "moa_score_pass_rate", "MoA presence          "),
    ("research_score", "research_score_pass_rate", "Research linkage      "),
    ("mitigation_score", "mitigation_score_pass_rate", "Mitigation pathway    "),
    ("anti_parroting", "anti_parroting_pass_rate", "Anti-parroting        "),
    ("word_count", "word_count_pass_rate", "Word count (250-400)  "),
]


def print_summary(base_agg, ft_agg, go, reason, winners):
    """Print a formatted comparison table to stdout."""
    print("\nBENCHMARK RESULTS: Base vs Fine-tuned MedGemma")
    print("=" * 58)
    print(f"{'Metric':<24} {'Base':>7}  {'Fine-tuned':>10}  {'Winner'}")
    print("-" * 58)
    for m_key, rate_key, display in METRIC_LABELS:
        b_pct = f"{base_agg[rate_key] * 100:.0f}%"
        f_pct = f"{ft_agg[rate_key] * 100:.0f}%"
        w = winners.get(m_key, "tie")
        winner_str = "Fine-tuned" if w == "finetuned" else ("Base" if w == "base" else "Tie")
        print(f"{display} {b_pct:>7}  {f_pct:>10}  {winner_str}")
    print("-" * 58)
    print(f"{'OVERALL':<24} {base_agg['overall_pass_rate'] * 100:.0f}%  {ft_agg['overall_pass_rate'] * 100:.0f}%")
    print(f"\nDECISION: {'GO' if go else 'NO-GO'}  {reason}\n")


def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description="Benchmark base MedGemma 4B-it vs LoRA fine-tuned adapter.")
    p.add_argument("--base-model", default="google/medgemma-4b-it")
    p.add_argument("--adapter-id", default="harshameghadri/medgemma-spatial-pathology-adapter")
    p.add_argument("--hf-token", default=os.environ.get("HF_TOKEN", ""), help="HuggingFace token")
    p.add_argument("--eval-data", default="data/lora_training/eval.jsonl")
    p.add_argument("--n-samples", type=int, default=10)
    p.add_argument("--output", default="outputs/lora_benchmark_results.json")
    p.add_argument("--skip-base", action="store_true")
    return p.parse_args()


def main():
    """Entry point: benchmark base vs fine-tuned and write JSON + console summary."""
    args = parse_args()
    samples = SYNTHETIC_SAMPLES[:min(args.n_samples, len(SYNTHETIC_SAMPLES))]
    print(f"Using {len(samples)} synthetic feature samples.")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results = {"timestamp": datetime.utcnow().isoformat() + "Z", "base_model": args.base_model, "adapter_id": args.adapter_id, "n_samples": len(samples), "results": {}, "go_decision": None, "go_reason": "", "winner_per_metric": {}}
    base_tokenizer = None
    base_model_obj = None

    if not args.skip_base:
        print(f"\nLoading base model: {args.base_model}")
        try:
            base_tokenizer, base_model_obj = load_base_model(args.base_model, args.hf_token)
            print("Evaluating base model ...")
            base_results = evaluate_model("base", base_tokenizer, base_model_obj, samples)
            results["results"]["base"] = base_results
            print(f"  Base overall pass rate: {base_results['aggregate']['overall_pass_rate'] * 100:.0f}%")
        except Exception as e:
            print(f"Base model loading failed: {e}")
            results["results"]["base"] = {"error": str(e)}
    else:
        results["results"]["base"] = {"skipped": True}

    print(f"\nLoading fine-tuned adapter: {args.adapter_id}")
    try:
        if base_model_obj is None:
            base_tokenizer, base_model_obj = load_base_model(args.base_model, args.hf_token)
        ft_tokenizer, ft_model = load_finetuned_model(base_model_obj, base_tokenizer, args.adapter_id, args.hf_token)
        print("Evaluating fine-tuned model ...")
        ft_results = evaluate_model("finetuned", ft_tokenizer, ft_model, samples)
        results["results"]["finetuned"] = ft_results
        print(f"  Fine-tuned overall pass rate: {ft_results['aggregate']['overall_pass_rate'] * 100:.0f}%")
    except Exception as e:
        print(f"\nAdapter not available: {e}")
        print("Run train_lora.ipynb on Kaggle first, then re-run benchmark.")
        results["results"]["finetuned"] = {"error": str(e), "adapter_not_yet_trained": True}

    base_ok = "aggregate" in results["results"].get("base", {})
    ft_ok = "aggregate" in results["results"].get("finetuned", {})
    if base_ok and ft_ok:
        go, reason, winners = compute_go_decision(results["results"]["base"]["aggregate"], results["results"]["finetuned"]["aggregate"])
        results["go_decision"] = go
        results["go_reason"] = reason
        results["winner_per_metric"] = winners
        print_summary(results["results"]["base"]["aggregate"], results["results"]["finetuned"]["aggregate"], go, reason, winners)
    else:
        results["go_decision"] = False
        results["go_reason"] = "Adapter not yet trained or insufficient data."
        print("\nPartial results — run training first for full benchmark.")

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()

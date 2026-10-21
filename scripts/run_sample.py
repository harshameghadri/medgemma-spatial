#!/usr/bin/env python3
"""Run one Visium sample through the same functions the Streamlit app calls.

Two stages so analysis (scanpy/squidpy/CellTypist) and generation (GPU) can use
different environments:

    python scripts/run_sample.py analyse --visium-dir data/real_breast --out runs/real_sample
    python scripts/run_sample.py report --out runs/real_sample --mode multimodal [--adapter runs/adapter_v2]

``analyse`` writes ``features.json`` and ``analysed.h5ad`` (keeps the H&E image).
``report`` writes ``report_<mode>_<base|adapter>.json`` with the raw report, the
grounding review (PubMed-checked citations) and the cleaned text the app shows.
Run one ``report`` per process: each loads its own model.
"""

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def analyse(visium_dir: Path, out: Path, resolution: float) -> dict:
    """Read a Space Ranger folder and run annotation, heterogeneity and feature extraction."""
    import squidpy as sq

    from src.streamlit_adapter import annotate_spatial_regions, build_features, calculate_spatial_heterogeneity

    adata = sq.read.visium(visium_dir)
    adata.var_names_make_unique()
    print(f"Loaded {adata.n_obs:,} spots, {adata.n_vars:,} genes")
    start = time.time()
    adata, annot = annotate_spatial_regions(adata, resolution=resolution, use_markers=True, tissue="Unknown")
    spatial = calculate_spatial_heterogeneity(adata)
    features = build_features(adata, annot, spatial)
    features["analysis_seconds"] = round(time.time() - start, 1)
    out.mkdir(parents=True, exist_ok=True)
    (out / "features.json").write_text(json.dumps(features, indent=2, default=str))
    adata.write_h5ad(out / "analysed.h5ad")
    print(json.dumps(features["annotation"]["cell_type_counts"], indent=2))
    for msg in features["errors"]:
        print("WARN", msg)
    return features


def report(out: Path, mode: str, adapter: str = None, verify: bool = True) -> dict:
    """Generate one report as the app would and attach the grounding review."""
    import torch

    from src.report_generation.medgemma_multimodal import generate_multimodal_report, generate_textonly_fallback
    from src.report_generation.pubmed_check import PubMedChecker, review_report

    features = json.loads((out / "features.json").read_text())
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    start = time.time()
    if mode == "multimodal":
        import anndata

        adata = anndata.read_h5ad(out / "analysed.h5ad")
        text, meta = generate_multimodal_report(adata, features, adapter_id=adapter)
    else:
        text, meta = generate_textonly_fallback(features, adapter_id=adapter)
    meta["elapsed_s"] = round(time.time() - start, 1)
    if torch.cuda.is_available():
        meta["peak_vram_gib"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
    checker = PubMedChecker(cache_path=str(out / "pubmed_cache.json"))
    review = review_report(text, meta["prompt"], verify_citations=verify, checker=checker)
    removed = [s for s in _sentences(text) if s not in review["clean_text"]]
    rec = {"mode": mode, "adapter": adapter, "raw": text, "clean": review["clean_text"],
           "removed_sentences": removed, "review": {k: v for k, v in review.items() if k != "clean_text"},
           "metadata": meta}
    name = f"report_{mode}_{'adapter' if adapter else 'base'}.json"
    (out / name).write_text(json.dumps(rec, indent=2, default=str))
    print(f"== {name}: {review['words']} words, {meta['elapsed_s']}s, "
          f"peak {meta.get('peak_vram_gib')} GiB, removed {len(removed)} sentences")
    print(review["clean_text"])
    return rec


def _sentences(text: str) -> list:
    from src.report_generation.faithfulness import split_sentences

    return [s for para in text.split("\n") if para.strip() for s in split_sentences(para)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="stage", required=True)
    a = sub.add_parser("analyse")
    a.add_argument("--visium-dir", type=Path, required=True)
    a.add_argument("--out", type=Path, required=True)
    a.add_argument("--resolution", type=float, default=0.5)
    r = sub.add_parser("report")
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--mode", choices=["text_only", "multimodal"], default="text_only")
    r.add_argument("--adapter", default=None)
    r.add_argument("--no-pubmed", action="store_true")
    args = ap.parse_args()
    if args.stage == "analyse":
        analyse(args.visium_dir, args.out, args.resolution)
    else:
        report(args.out, args.mode, args.adapter, not args.no_pubmed)


if __name__ == "__main__":
    main()

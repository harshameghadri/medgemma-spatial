#!/usr/bin/env python3
"""Write NOTICE.md into an adapter folder trained on Llama 3.1 output (docs/REVIEW.md L15)."""

import argparse
from pathlib import Path

LICENCES = {
    "llama": "Llama 3.1 Community License Agreement (https://www.llama.com/llama3_1/license/)",
    "medgemma": "Health AI Developer Foundations Terms of Use (https://developers.google.com/health-ai-developer-foundations/terms)",
}


def notice(name: str, base: str, teacher: str, licences: list) -> str:
    """NOTICE.md text for one adapter."""
    lines = [f"# {name}", "", "Built with Llama.", "",
             f"LoRA adapter on `{base}`, fine-tuned on pathology report drafts written by `{teacher}`",
             "and picked by a domain expert (docs/REVIEW.md L15). Research use only; not for clinical decisions.",
             "", "Licences that apply:"]
    lines += [f"- {LICENCES[k]}" for k in licences]
    lines += ["", "If this model is distributed, its name must begin with \"Llama\" (Llama 3.1 Community License, section 1.b.i)."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adapter-dir", type=Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--teacher", default="meta-llama/Llama-3.1-8B-Instruct")
    ap.add_argument("--licences", nargs="+", choices=sorted(LICENCES), default=["llama"])
    args = ap.parse_args(argv)
    if not args.name.startswith("Llama"):
        raise SystemExit("models trained on Llama output must be named starting with 'Llama'")
    args.adapter_dir.mkdir(parents=True, exist_ok=True)
    (args.adapter_dir / "NOTICE.md").write_text(notice(args.name, args.base, args.teacher, args.licences))


if __name__ == "__main__":
    main()

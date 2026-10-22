#!/usr/bin/env python3
"""Check any model's report against the prompt it was given; no GPU and no model needed.

Flags what the input never told the model: numbers, cell types, literature
citations (each looked up in PubMed unless --offline), and 5+ word runs copied
verbatim from the prompt. Prints the cleaned report the app would show.

    python scripts/check_report.py --report report.md --prompt prompt.txt
    python scripts/check_report.py --report report.md --features features.json   # rebuild the prompt
    python scripts/check_report.py --report report.md --prompt prompt.txt --json --offline
"""

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.report_generation.pubmed_check import review_report  # noqa: E402


def check(report: str, prompt: str, verify_citations: bool = True) -> dict:
    """Grounding review of one report, as a JSON-ready dict."""
    review = review_report(report, prompt, verify_citations=verify_citations)
    kept = review.pop("clean_text")
    return {**review, "clean_text": kept, "changed": kept != report.strip()}


def render(result: dict) -> str:
    """Human-readable summary of a check() result."""
    lines = [f"words: {result['words']}   sections: {', '.join(result['sections']) or 'none'}"]
    for key, label in (("unsupported_numbers", "numbers not in the prompt"),
                       ("unsupported_cell_types", "cell types not in the prompt"),
                       ("copied_spans", "copied verbatim from the prompt")):
        items = result[key]
        lines.append(f"{label}: {len(items)}" + "".join(f"\n  - {x}" for x in items))
    lines.append(f"citations: {len(result['pubmed'])}")
    for c in result["pubmed"]:
        where = f" -> PMID {c['pmid']}" if c.get("pmid") else ""
        lines.append(f"  - [{c['status']}] {c['text']}{where}")
    lines.append("\n--- report after removing ungrounded sentences ---" if result["changed"] else "\n--- report ---")
    lines.append(result["clean_text"])
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", type=Path, required=True, help="generated report (text or markdown)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--prompt", type=Path, help="the exact prompt the model received")
    src.add_argument("--features", type=Path, help="features.json; the text-only prompt is rebuilt from it")
    ap.add_argument("--offline", action="store_true", help="do not query PubMed")
    ap.add_argument("--json", action="store_true", help="print JSON instead of a summary")
    args = ap.parse_args(argv)

    if args.prompt:
        prompt = args.prompt.read_text()
    else:
        from src.report_generation.prompt_builder import generate_medgemma_prompt
        prompt = generate_medgemma_prompt(json.loads(args.features.read_text()))
    result = check(args.report.read_text(), prompt, verify_citations=not args.offline)
    print(json.dumps(result, indent=2, default=str) if args.json else render(result))


if __name__ == "__main__":
    main()

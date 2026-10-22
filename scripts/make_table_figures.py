#!/usr/bin/env python3
"""Render the article's three tables as HTML pages (docs/figures/table1-3.html) for PNG export."""

import html
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "docs" / "figures"

TABLES = {
    "table1": {
        "title": "Distillation passed the automatic metrics but not the blinded read",
        "head": ["Model", "Mechanism score, base → student", "Pre-registered threshold", "Automatic verdict",
                 "Blinded read: errors, base → student (7 cases)", "Overall"],
        "rows": [["Llama 3.1 8B", "7.03 → 6.52", "≥ 6.53", "Fail (by 0.01)", "8 → 5", "Fail"],
                 ["MedGemma 4B", "3.13 → 7.00", "≥ 4.13", "Pass", "3 → 4", "Fail"]],
        "num": [1, 2, 4],
        "verdict": [3, 5],
        "widths": ["15%", "21%", "16%", "15%", "21%", "12%"],
    },
    "table2": {
        "title": "Recurring kinds of biology error in the first blinded read",
        "head": ["Error kind", "Example from a report", "Why it is wrong"],
        "rows": [
            ["Inhibitory axis described as activating",
             "“activation of T cells by Stromal CAFs through the PD-1/PD-L1 axis”",
             "PD-1 engagement suppresses T-cell receptor signalling"],
            ["Immunosuppression proposed as anti-tumour therapy",
             "infliximab “enhancing the anti-tumour activity”; “enhancing the activity of regulatory T cells”",
             "Dampening immunity favours tumour escape"],
            ["Receptor placed on the wrong cell",
             "a CD40 agonist “enhancing the activation of Epithelial cells”",
             "CD40 agonism acts on antigen-presenting cells"],
            ["Cell types absent from the data",
             "regulatory T cells, myeloid-derived suppressor cells",
             "Not in the annotation"],
            ["Therapy that ignores the phenotype",
             "checkpoint blockade to increase infiltration of an immune desert; no myeloid-directed agent in a "
             "macrophage-dominant tissue",
             "Checkpoint blockade releases T cells already present"],
        ],
        "num": [],
        "verdict": [],
        "quote": 1,
        "widths": ["24%", "44%", "32%"],
    },
    "table3": {
        "title": "Prompt conditions compared (MedGemma 4B)",
        "head": ["Prompt condition (MedGemma 4B)", "Flagged reports (n = 31)", "Mechanism score, own words",
                 "Words copied from inserted text (median)", "Blinded read: errors (n = 7)",
                 "Blinded read: generic reports"],
        "rows": [["Plain prompt", "6", "3.13", "—", "3", "7 of 7"],
                 ["+ two retrieved approved reports", "1", "1.23", "86%", "not run (failed automatic criteria)", "—"],
                 ["+ curated axis table and rules", "0", "6.87", "0%", "2", "0 of 7"]],
        "num": [1, 2, 3, 4, 5],
        "verdict": [],
        "highlight": 2,
        "bad": 1,
    },
}

PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{title}</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
<style>
 * {{ box-sizing: border-box; }}
 html, body {{ margin: 0; background: #fff; }}
 body {{ font-family: Inter, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: #17212b; }}
 .canvas {{ width: 1200px; padding: 32px 40px 30px; background: #fff; }}
 h1 {{ font-size: 21px; margin: 0 0 16px; letter-spacing: -0.005em; }}
 table {{ width: 100%; border-collapse: collapse; font-size: 16px; line-height: 1.45; font-variant-numeric: tabular-nums; }}
 th {{ text-align: left; vertical-align: bottom; font-weight: 600; font-size: 13px; letter-spacing: 0.03em;
       text-transform: uppercase; color: #5b6875; padding: 0 14px 10px; border-bottom: 2px solid #17212b; }}
 td {{ padding: 13px 14px; border-bottom: 1px solid #e3e8ed; vertical-align: top; }}
 td.num, th.num {{ text-align: center; }}
 td.first {{ font-weight: 600; }}
 td.num {{ white-space: nowrap; }}
 td.quote {{ color: #2d3a46; }}
 tr.hl td {{ background: #e6f4f1; font-weight: 600; }}
 tr.bad td {{ color: #7a2e0e; }}
 .fail {{ color: #b42318; font-weight: 700; }}
 .pass {{ color: #067647; font-weight: 700; }}
 .foot {{ margin-top: 12px; font-size: 12.5px; color: #5b6875; text-align: right; }}
</style></head><body><div class="canvas">
<table>{colgroup}<thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>
<div class="foot">github.com/harshameghadri/medgemma-spatial</div>
</div></body></html>
"""


def cell(t, value, ci):
    cls = []
    if ci == 0:
        cls.append("first")
    if ci in t["num"]:
        cls.append("num")
    if ci == t.get("quote"):
        cls.append("quote")
    v = html.escape(value)
    if ci in t["verdict"]:
        v = f'<span class="{"pass" if value.startswith("Pass") else "fail"}">{v}</span>'
    return f'<td class="{" ".join(cls)}">{v}</td>'


def render(t):
    head = "".join(f'<th class="{"num" if i in t["num"] else ""}">{html.escape(h)}</th>' for i, h in enumerate(t["head"]))
    rows = []
    for ri, r in enumerate(t["rows"]):
        cls = "hl" if ri == t.get("highlight") else "bad" if ri == t.get("bad") else ""
        rows.append(f'<tr class="{cls}">' + "".join(cell(t, v, ci) for ci, v in enumerate(r)) + "</tr>")
    colgroup = ("<colgroup>" + "".join(f'<col style="width:{w}">' for w in t["widths"]) + "</colgroup>"
                if "widths" in t else "")
    return PAGE.format(title=html.escape(t["title"]), head=head, rows="".join(rows), colgroup=colgroup)


def main():
    for name, t in TABLES.items():
        (OUT / f"{name}.html").write_text(render(t))
        print(f"wrote {OUT / name}.html")


if __name__ == "__main__":
    main()

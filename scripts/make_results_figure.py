#!/usr/bin/env python3
"""Render docs/figures/results.html: six SVG panels drawn from the run outputs (docs/REVIEW.md L17-L23).

    python scripts/make_results_figure.py --data figdata.json
"""

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
INK, MUTED, GRID = "#17212b", "#5b6875", "#e3e8ed"
PLAIN, FULL, EVID, BAD = "#9aa5b1", "#d97706", "#0f766e", "#b42318"
W, H = 470, 300
L, R, T, B = 56, 16, 18, 52


def text(x, y, s, size=13, color=INK, anchor="middle", weight=400):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" text-anchor="{anchor}" '
            f'font-weight="{weight}">{s}</text>')


def axis(ymax, step, label):
    out = []
    ph = H - T - B
    v = 0
    while v <= ymax + 1e-9:
        y = T + ph - ph * v / ymax
        out.append(f'<line x1="{L}" x2="{W - R}" y1="{y:.1f}" y2="{y:.1f}" stroke="{GRID}"/>')
        out.append(text(L - 8, y + 4, f"{v:g}", 12, MUTED, "end"))
        v += step
    out.append(f'<text transform="translate(16,{T + ph / 2}) rotate(-90)" font-size="12" fill="{MUTED}" '
               f'text-anchor="middle">{label}</text>')
    return out


def bars(groups, series, ymax, step, ylabel, fmt="{:g}"):
    """groups: [(label, [values per series])]; series: [(name, colour)]."""
    out = axis(ymax, step, ylabel)
    ph, pw = H - T - B, W - L - R
    gw = pw / len(groups)
    bw = min(46, gw * 0.7 / len(series))
    for gi, (glabel, vals) in enumerate(groups):
        cx = L + gw * (gi + 0.5)
        x0 = cx - bw * len(series) / 2
        for si, v in enumerate(vals):
            if v is None:
                out.append(text(x0 + bw * (si + 0.5), T + ph - 6, "n/a", 11, MUTED))
                continue
            h = ph * v / ymax
            x = x0 + bw * si
            out.append(f'<rect x="{x + 2:.1f}" y="{T + ph - h:.1f}" width="{bw - 4:.1f}" height="{h:.1f}" '
                       f'rx="3" fill="{series[si][1]}"/>')
            out.append(text(x + bw / 2, T + ph - h - 6, fmt.format(v), 12.5, INK, weight=600))
        for li, line in enumerate(glabel.split("|")):
            out.append(text(cx, T + ph + 20 + 15 * li, line, 12.5, INK))
    return out


def legend(series):
    items = []
    for name, colour in series:
        items.append(f'<span class="key"><i style="background:{colour}"></i>{name}</span>')
    return "".join(items)


def strip(rows, ylabel):
    """rows: [(label, values 0-1, colour)] as a dot strip per condition with a median line."""
    out = axis(100, 25, ylabel)
    ph, pw = H - T - B, W - L - R
    gw = pw / len(rows)
    for gi, (label, vals, colour) in enumerate(rows):
        cx = L + gw * (gi + 0.5)
        svals = sorted(vals)
        for i, v in enumerate(svals):
            jitter = ((i * 37) % 11 - 5) * 4.2
            y = T + ph - ph * v
            out.append(f'<circle cx="{cx + jitter:.1f}" cy="{y:.1f}" r="4.2" fill="{colour}" fill-opacity="0.75"/>')
        med = svals[len(svals) // 2]
        y = T + ph - ph * med
        out.append(f'<line x1="{cx - 34:.1f}" x2="{cx + 34:.1f}" y1="{y:.1f}" y2="{y:.1f}" stroke="{INK}" stroke-width="2.5"/>')
        out.append(text(cx + 40, y + 4, f"{med * 100:.0f}%", 12.5, INK, "start", 700))
        for li, line in enumerate(label.split("|")):
            out.append(text(cx, T + ph + 20 + 15 * li, line, 12.5, INK))
    return out


def panel(tag, title, note, svg, keys=""):
    return (f'<figure class="panel"><div class="tag">{tag}</div><h2>{title}</h2><p>{note}</p>'
            f'<svg viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Inter, -apple-system, Helvetica, Arial, sans-serif">'
            f'{"".join(svg)}</svg><div class="keys">{keys}</div></figure>')


def build(d):
    s_base = [("Base model", PLAIN), ("Distilled student", FULL)]
    a = panel("A", "Passing the metrics did not mean passing the read",
              "Round 2, blinded read of 7 cases: biology errors before and after distillation.",
              bars([("Llama 3.1 8B|automatic: fail by 0.01", [d["read1_errors"]["llama_base"], d["read1_errors"]["llama_student"]]),
                    ("MedGemma 4B|automatic: pass", [d["read1_errors"]["medgemma_base"], d["read1_errors"]["medgemma_student"]])],
                   s_base, 10, 2, "Biology errors"), legend(s_base))
    kinds = d["kinds_reports"]
    order = [("therapy", "Therapy ignores|phenotype"), ("cellbio", "Wrong cell|biology"), ("direction", "Signalling|direction"),
             ("unsupported", "Not supported|by the data"), ("drug", "Wrong or|invented drug")]
    b = panel("B", "The errors had a shape",
              "First blinded read (28 reports): reports tagged with each error kind; a report can carry several.",
              bars([(lab, [kinds.get(k, 0)]) for k, lab in order], [("Reports", BAD)], 8, 2, "Reports"))
    c = panel("C", "Fewer reports flagged by the biology check",
              "MedGemma 4B, 31 held-out prompts: reports with at least one flag.",
              bars([("Plain", [d["mg"]["base"]["biology_flagged_reports"]]),
                    ("+ retrieved|reports", [d["mg"]["adapter"]["biology_flagged_reports"]]),
                    ("+ curated|axes", [d["mge"]["adapter"]["biology_flagged_reports"]])],
                   [("Flagged reports", INK)], 8, 2, "Flagged reports (of 31)"))
    c = c.replace(f'fill="{INK}"/>', "fill=\"X\"/>")
    for colour in (PLAIN, FULL, EVID):
        c = c.replace('fill="X"/>', f'fill="{colour}"/>', 1)
    s_moa = [("All mechanism terms", "#c3cbd3"), ("Own words only", INK)]
    dm = panel("D", "Retrieval inflated the raw score",
               "Mechanism terms per report; 'own words' excludes runs of 5+ words copied from the inserted text.",
               bars([("Plain", [d["mge"]["base"]["moa_mean"], d["mge"]["base"]["moa_novel_mean"]]),
                     ("+ retrieved|reports", [d["mg"]["adapter"]["moa_mean"], d["mg"]["adapter"]["moa_novel_mean"]]),
                     ("+ curated|axes", [d["mge"]["adapter"]["moa_mean"], d["mge"]["adapter"]["moa_novel_mean"]])],
                    s_moa, 8, 2, "Mechanism terms per report", "{:.1f}"), legend(s_moa))
    e = panel("E", "Retrieved reports were copied, curated facts were not",
              "Share of each report's words in 5+ word runs copied from the inserted text; one dot per report, line = median.",
              strip([("MedGemma|+ retrieved reports", d["copy_full"], FULL),
                     ("Llama|+ retrieved reports", d["copy_full_llama"], FULL),
                     ("MedGemma|+ curated axes", d["copy_evidence"], EVID)], "Report copied (%)"))
    s_read = [("Plain prompt", PLAIN), ("+ curated axes (served)", EVID)]
    f = panel("F", "The blinded read agreed",
              "Second blinded read, same 7 cases: biology errors and reports marked generic.",
              bars([("Biology errors", [3, 2]), ("Generic reports", [7, 0])], s_read, 8, 2, "Reports (of 7)"),
              legend(s_read))
    return [a, b, c, dm, e, f]


PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>MedGemma Spatial results</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&family=JetBrains+Mono:wght@500&display=swap" rel="stylesheet">
<style>
 * {{ box-sizing: border-box; }}
 html, body {{ margin: 0; background: #ffffff; }}
 body {{ font-family: Inter, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: #17212b; }}
 .canvas {{ width: 1600px; padding: 44px 52px 36px; background: #fff; }}
 h1 {{ font-size: 30px; margin: 0 0 6px; letter-spacing: -0.01em; }}
 .sub {{ font-size: 17px; color: #5b6875; margin: 0 0 26px; }}
 .grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 26px 30px; }}
 .panel {{ margin: 0; position: relative; }}
 .tag {{ position: absolute; left: 0; top: 0; font: 700 15px "JetBrains Mono", monospace; color: #0f766e; }}
 .panel h2 {{ font-size: 17px; margin: 0 0 4px 24px; line-height: 1.25; }}
 .panel p {{ font-size: 13.5px; color: #5b6875; margin: 0 0 8px 24px; line-height: 1.4; min-height: 38px; }}
 .keys {{ display: flex; gap: 16px; margin-left: 56px; font-size: 12.5px; color: #2d3a46; min-height: 18px; }}
 .key i {{ display: inline-block; width: 12px; height: 12px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; }}
 .foot {{ margin-top: 20px; font-size: 12.5px; color: #5b6875; display: flex; justify-content: space-between; }}
</style></head><body><div class="canvas">
<h1>What moved the biology error rate</h1>
<p class="sub">MedGemma 4B reports from spatial transcriptomics features: distillation, retrieval of approved reports, and curated ligand&#8211;receptor facts</p>
<div class="grid">{panels}</div>
<div class="foot"><span>31 held-out prompts (synthetic feature summaries) and blinded reads of 7 cases by one domain reviewer. Criteria fixed before generation; details in docs/REVIEW.md.</span><span>github.com/harshameghadri/medgemma-spatial</span></div>
</div></body></html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=str(REPO / "docs/figures/results.html"))
    args = ap.parse_args(argv)
    d = json.loads(Path(args.data).read_text())
    Path(args.out).write_text(PAGE.format(panels="".join(build(d))))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

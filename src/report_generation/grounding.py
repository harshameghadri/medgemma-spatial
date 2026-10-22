"""Inference-time grounding: curated ligand-receptor evidence plus retrieved reviewer-approved reports.

Two blocks are inserted into the serving prompt before its instructions:

* **Evidence**: axes from ``data/knowledge/axes.tsv`` whose sending and receiving
  cell types are both in the sample, ranked by the sample's context (myeloid,
  stromal, desert, immune, TLS). Every pair exists in the OmniPath snapshot
  (tested); the sender, receiver and effect on the receiving cell are curated.
  Three therapy rules from the blinded read (docs/REVIEW.md L18) follow it.
* **Examples**: the two most similar reports from ``approved_reports.jsonl``,
  teacher drafts the domain reviewer picked (distill_v1) minus those the biology
  check flags. Similarity uses only the data-summary lines of the prompt.

``ground_prompt`` returns the new prompt and the inserted text, so copying from
the examples can be measured separately from copying from the original prompt.
"""

import csv
import functools
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .faithfulness import CELL_TYPE_PATTERNS, cell_types_in

KNOWLEDGE = Path(__file__).resolve().parents[2] / "data" / "knowledge"
AXES = KNOWLEDGE / "axes.tsv"
APPROVED = KNOWLEDGE / "approved_reports.jsonl"
MAX_AXES = 6
N_EXAMPLES = 2

SUMMARY = re.compile(r"SPATIAL TRANSCRIPTOMICS DATA SUMMARY:\n(.*?)\n\n", re.S)
FIELD = re.compile(r"^- (?P<key>[^:(]+?)(?: \([^)]*\))?: (?P<value>.*)$", re.M)
NUMBER = re.compile(r"\((-?\d+(?:\.\d+)?)\)")
INSERT_BEFORE = re.compile(r"^Generate a ", re.M)

MYELOID = {"macrophage", "neutrophil"}
CONTEXT_WEIGHT = {"myeloid": 2, "stromal": 2, "desert": 2, "tls": 2, "immune": 1}

RULES = (
    "- If you name a ligand-receptor axis, choose one listed above and keep its direction and effect.",
    "- Checkpoint blockade releases T cells already present; it does not bring T cells into an immune desert.",
    "- Do not propose immunosuppression (anti-TNF, corticosteroids, dampening T cells) as anti-tumour therapy.",
)
MYELOID_RULE = "- Macrophages or neutrophils dominate: include a myeloid-directed agent in the therapy."


@functools.lru_cache(maxsize=1)
def load_axes(path: str = str(AXES)) -> Tuple[Dict[str, str], ...]:
    """Curated axes with sender/receiver/context split into lists."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for key in ("sender", "receiver", "context"):
            r[key] = tuple(x for x in r[key].split(",") if x)
    return tuple(rows)


@functools.lru_cache(maxsize=1)
def load_approved(path: str = str(APPROVED)) -> Tuple[Dict, ...]:
    """Reviewer-approved reports with their parsed data summaries."""
    if not Path(path).exists():
        return ()
    return tuple(json.loads(line) for line in open(path))


def parse_summary(prompt: str) -> Dict:
    """Phenotype, ordered cell types and spatial metrics from a serving prompt's data summary."""
    m = SUMMARY.search(prompt)
    fields = {f.group("key").strip(): f.group("value").strip() for f in FIELD.finditer(m.group(1) if m else "")}
    cells = [c.strip() for c in fields.get("Cell types detected", "").split(",") if c.strip()]

    def num(key):
        n = NUMBER.search(fields.get(key, ""))
        return float(n.group(1)) if n else None

    return {"phenotype": fields.get("Tissue phenotype", ""), "cell_types": cells,
            "morans_i": num("Spatial autocorrelation"), "entropy": num("Neighbourhood entropy")}


def canonical(label: str) -> Optional[str]:
    """Canonical cell-type key ('Stromal CAF' -> 'fibroblast'), None if unknown."""
    hits = cell_types_in(label)
    return next((k for k in CELL_TYPE_PATTERNS if k in hits), None)


def contexts(summary: Dict) -> set:
    """Context tags that rank the axes for this sample."""
    keys = [canonical(c) for c in summary["cell_types"]]
    top2, top3 = set(keys[:2]), set(keys[:3])
    pheno = summary["phenotype"].lower()
    tags = set()
    if top2 & MYELOID:
        tags.add("myeloid")
    if "stromal" in pheno or "fibrotic" in pheno or "fibroblast" in top2:
        tags.add("stromal")
    if "desert" in pheno:
        tags.add("desert")
    if "hot" in pheno or "immune-rich" in pheno:
        tags.add("immune")
    if "b_cell" in top3 and "t_cell" in top3:
        tags.add("tls")
    return tags


def select_axes(summary: Dict, k: int = MAX_AXES) -> List[Dict[str, str]]:
    """Axes whose sender and receiver are both present, context matches first, table order otherwise."""
    present = {canonical(c) for c in summary["cell_types"]} - {None}
    tags = contexts(summary)
    eligible = [(i, a) for i, a in enumerate(load_axes())
                if present & set(a["sender"]) and present & set(a["receiver"])]
    eligible.sort(key=lambda ia: (-sum(CONTEXT_WEIGHT[t] for t in tags & set(ia[1]["context"])), ia[0]))
    return [a for _, a in eligible[:k]]


def _names(keys, summary) -> str:
    by_key = {canonical(c): c for c in summary["cell_types"]}
    return ", ".join(by_key[k] for k in keys if k in by_key)


def evidence_block(summary: Dict) -> str:
    """Candidate axes and therapy rules for the prompt ('' if no axis fits)."""
    axes = select_axes(summary)
    if not axes:
        return ""
    lines = ["CANDIDATE LIGAND-RECEPTOR AXES (pairs in OmniPath; sender -> receiver; effect on the receiver):"]
    for a in axes:
        lines.append(f"- {a['axis']}: {_names(a['sender'], summary)} -> {_names(a['receiver'], summary)}; "
                     f"{a['effect']}; agent: {a['therapy']}")
    lines += ["", "THERAPY RULES:", *RULES]
    if "myeloid" in contexts(summary):
        lines.append(MYELOID_RULE)
    return "\n".join(lines)


def similarity(a: Dict, b: Dict) -> float:
    """Data-summary similarity: phenotype, leading cell types and spatial metrics."""
    score = 2.0 * (a["phenotype"] == b["phenotype"])
    ka = [canonical(c) for c in a["cell_types"]]
    kb = [canonical(c) for c in b["cell_types"]]
    score += len(set(ka[:3]) & set(kb[:3])) / 3 + (ka[:1] == kb[:1])
    for key, scale in (("morans_i", 0.5), ("entropy", 1.0)):
        if a[key] is not None and b[key] is not None:
            score -= abs(a[key] - b[key]) / scale
    return score


def retrieve(summary: Dict, k: int = N_EXAMPLES) -> List[Dict]:
    """The k most similar approved reports (ties broken by id)."""
    pool = sorted(load_approved(), key=lambda r: (-similarity(summary, r["summary"]), r["id"]))
    return pool[:k]


def examples_block(summary: Dict) -> str:
    """Retrieved reports, framed as structure and reasoning to follow, not text to reuse."""
    hits = retrieve(summary)
    if not hits:
        return ""
    parts = ["EXAMPLE REPORTS FOR OTHER SAMPLES (approved by a pathologist reviewer). Follow their structure "
             "and level of specificity; do not reuse their wording, and discuss only this sample's cell types."]
    for i, r in enumerate(hits, 1):
        parts.append(f"--- Example {i} ({r['summary']['phenotype']}) ---\n{r['text'].strip()}")
    return "\n\n".join(parts)


def ground_prompt(prompt: str, evidence: bool = True, examples: bool = True) -> Tuple[str, str]:
    """(grounded prompt, inserted text); the prompt is unchanged if it has no data summary."""
    summary = parse_summary(prompt)
    if not summary["cell_types"]:
        return prompt, ""
    blocks = [b for b in (evidence and evidence_block(summary), examples and examples_block(summary)) if b]
    if not blocks:
        return prompt, ""
    added = "\n\n".join(blocks) + "\n\n"
    m = INSERT_BEFORE.search(prompt)
    at = m.start() if m else len(prompt)
    return prompt[:at] + added + prompt[at:], added

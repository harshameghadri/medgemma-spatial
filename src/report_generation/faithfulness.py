"""Grounding checks for generated reports: citations, numbers and cell types.

A claim is *unsupported* when the report states something the prompt never
gave the model: a number that is not in the prompt, a cell type the prompt did
not list, or a literature citation (the prompt contains none). The same checks
score generated reports (``compare_adapter.py``), clean LoRA targets
(``build_lora_dataset.py``) and flag reports in the app.
"""

import re
from typing import Dict, List, Tuple

ABBREVIATIONS = ("et al.", "e.g.", "i.e.", "vs.", "approx.", "cf.", "Fig.", "ca.")
_DOT = "․"

CITATION = re.compile(
    r"(?P<surname>[A-Z][A-Za-zÀ-ſ'\-]+)(?:\s+and\s+[A-Z][A-Za-z'\-]+)?\s+et\s+al\.?,?"
    r"\s*\(?(?P<middle>[^()]{0,80}?)\(?(?P<year>(?:19|20)\d\d)[a-z]?\)?"
)
CITATION_MARKER = re.compile(r"\bet al\b|\b(?:19|20)\d\d\b")

NUMBER = re.compile(r"(?<![\w\-/.–])\d+(?:\.\d+)?%?(?![\w])")
TRIVIAL_NUMBER = re.compile(r"[1-5]|(?:19|20)\d\d")

_I = re.IGNORECASE
CELL_TYPE_PATTERNS: Dict[str, re.Pattern] = {
    "epithelial": re.compile(r"\bepitheli\w*", _I),
    "t_cell": re.compile(r"\bT[\s-](?:(?i:cells?|lymphocytes?))\b|(?i:\bcytotoxic lymphocytes?\b)"),
    "b_cell": re.compile(r"\bB[\s-](?:(?i:cells?|lymphocytes?))\b"),
    "nk_cell": re.compile(r"\bNK\b|\bnatural killer\b", _I),
    "macrophage": re.compile(r"(?i:\bmacrophages?\b)|\bTAMs?\b"),
    "dendritic_cell": re.compile(r"(?i:\bdendritic\b)|\bc?DCs?\b"),
    "plasma_cell": re.compile(r"\bplasma(?:\s+cells?|blasts?)\b", _I),
    "neutrophil": re.compile(r"\bneutrophils?\b|\bTANs?\b", _I),
    "endothelial": re.compile(r"\bendotheli\w*", _I),
    "fibroblast": re.compile(r"\b(?:myo)?fibroblasts?\b|\bCAFs?\b", _I),
    "mast_cell": re.compile(r"\bmast cells?\b", _I),
    "monocyte": re.compile(r"\bmonocytes?\b", _I),
    "treg": re.compile(r"\bTregs?\b|(?i:\bregulatory\s+)T[\s-](?i:cells?)\b|\bFOXP3\b"),
    "mdsc": re.compile(r"\bMDSCs?\b|\bmyeloid[\s-]derived suppressor\b", _I),
}

HEADING = re.compile(r"^\s*(?:#{1,6}\s*)?\*\*(?P<title>[^*]+?)\*\*:?\s*(?P<rest>.*)$|^\s*#{1,6}\s*(?P<md>.+?)\s*$")
SECTIONS = ("TISSUE MICROENVIRONMENT", "MECHANISM OF ACTION", "THERAPEUTIC IMPLICATIONS", "VALIDATION CAVEATS")
_SECTION_RULES = (
    (re.compile(r"PATHOLOGY REPORT|SPECIMEN ANALYSIS|VISIUM HD|TISSUE ANALYSIS$", _I), None),
    (re.compile(r"MECHANIS", _I), "MECHANISM OF ACTION"),
    (re.compile(r"THERAPEUTIC|TREATMENT", _I), "THERAPEUTIC IMPLICATIONS"),
    (re.compile(r"QUALITY|VALIDATION|UNCERTAINTY|TECHNICAL|LIMITATION|CAVEAT", _I), "VALIDATION CAVEATS"),
    (re.compile(r"MICROENVIRONMENT|ARCHITECTURE|CONTEXTURE|DIAGNOSTIC|SPECIMEN CHARACTERISTICS|COMPOSITION", _I),
     "TISSUE MICROENVIRONMENT"),
)


def split_sentences(text: str) -> List[str]:
    """Split prose into sentences without breaking 'et al.', 'e.g.' and similar."""
    protected = text
    for abbr in ABBREVIATIONS:
        protected = protected.replace(abbr, abbr.replace(".", _DOT))
    parts = re.split(r"(?<=[.!?])\s+", protected.strip())
    return [p.replace(_DOT, ".") for p in parts if p.strip()]


def cell_types_in(text: str) -> set:
    """Canonical cell types mentioned in a text."""
    return {name for name, pat in CELL_TYPE_PATTERNS.items() if pat.search(text)}


def unsupported_numbers(text: str, prompt: str) -> List[str]:
    """Standalone numbers in the report that are absent from the prompt (gene names like CXCL9 are ignored)."""
    allowed = set(NUMBER.findall(prompt))
    return [n for n in NUMBER.findall(text) if n not in allowed and not TRIVIAL_NUMBER.fullmatch(n)]


def unsupported_cell_types(text: str, prompt: str) -> List[str]:
    """Cell types named in the report that the prompt did not list."""
    return sorted(cell_types_in(text) - cell_types_in(prompt))


def extract_citations(text: str) -> List[Dict]:
    """Author-year citations, e.g. 'Helmink et al. (Nature Medicine, 2020)'."""
    out = []
    for m in CITATION.finditer(text):
        journal = re.sub(r"[\s,;:.]+$", "", m.group("middle")).strip(" ,(")
        out.append({"surname": m.group("surname"), "year": int(m.group("year")),
                    "journal": journal or None, "text": m.group(0)})
    return out


def sentence_issues(sentence: str, prompt: str) -> List[str]:
    """Reasons a sentence is not grounded in the prompt (empty when it is)."""
    issues = []
    if CITATION_MARKER.search(sentence):
        issues.append("citation")
    if unsupported_numbers(sentence, prompt):
        issues.append("number")
    if unsupported_cell_types(sentence, prompt):
        issues.append("cell_type")
    return issues


def score_report(text: str, prompt: str) -> Dict:
    """Grounding summary for one generated report."""
    return {
        "words": len(text.split()),
        "citations": [c["text"] for c in extract_citations(text)],
        "citation_sentences": sum("citation" in sentence_issues(s, prompt) for s in split_sentences(text)),
        "unsupported_numbers": unsupported_numbers(text, prompt),
        "unsupported_cell_types": unsupported_cell_types(text, prompt),
        "sections": [s for s in SECTIONS if s in text.upper()],
    }


def canonical_section(title: str):
    """Map a free-form heading to one of SECTIONS; None for report titles, '' when unknown."""
    for pattern, section in _SECTION_RULES:
        if pattern.search(title):
            return section
    return ""


def _parse_line(line: str) -> Tuple[object, str]:
    m = HEADING.match(line)
    if not m:
        return False, line
    title = (m.group("title") or m.group("md") or "").strip()
    rest = (m.group("rest") or "").strip()
    return canonical_section(title), rest


def scrub_report(text: str, prompt: str) -> Tuple[str, List[Dict]]:
    """Drop ungrounded sentences and normalise headings to SECTIONS.

    Returns the cleaned report and the removed sentences with reasons. Unknown
    headings fold their text into the current section; sections left empty are
    dropped.
    """
    sections: List[Tuple[str, List[str]]] = []
    removed: List[Dict] = []
    current = None

    def add(paragraph: str):
        nonlocal current
        if current is None:
            current = ("", [])
            sections.append(current)
        kept = []
        for sentence in split_sentences(paragraph):
            issues = sentence_issues(sentence, prompt)
            if issues:
                removed.append({"sentence": sentence, "issues": issues})
            else:
                kept.append(sentence)
        if kept:
            current[1].append(" ".join(kept))

    paragraph: List[str] = []
    for line in text.splitlines() + [""]:
        section, rest = _parse_line(line)
        if section is False and line.strip():
            paragraph.append(line.strip())
            continue
        if paragraph:
            add(" ".join(paragraph))
            paragraph = []
        if section:
            current = (section, [])
            sections.append(current)
        if section and rest:
            paragraph.append(rest)

    merged: Dict[str, List[str]] = {}
    order: List[str] = []
    for name, paras in sections:
        if name not in merged:
            merged[name], _ = [], order.append(name)
        merged[name].extend(paras)
    blocks = []
    for name in order:
        if not merged[name]:
            continue
        body = "\n\n".join(merged[name])
        blocks.append(f"**{name}**\n\n{body}" if name else body)
    return "\n\n".join(blocks).strip(), removed

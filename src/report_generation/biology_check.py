"""Rule and database checks for the biology errors found in the blinded read (docs/REVIEW.md L18).

Each check returns findings ``{"check", "sentence", "detail"}``; a report with no
findings is not thereby correct, only free of these specific errors.

* ``unannotated_cell_type`` - a cell type the prompt never listed (Tregs, MDSCs...).
* ``inhibitory_axis_as_activating`` - an inhibitory axis (PD-1/PD-L1, CTLA-4,
  CD47-SIRPa, TGF-b, IL-10, LAG-3, TIM-3, TIGIT) described as activating or
  recruiting immune cells, outside a blockade sentence. OmniPath signs describe
  the receptor (PD-L1 *stimulates* PD-1), not the receiving cell (the T cell is
  inhibited), so the cell-level effect comes from ``INHIBITORY_GENES``.
* ``immunosuppressive_therapy`` - immunosuppressants (WHO ATC L04, systemic
  corticosteroids H02) or an explicitly immunosuppressive aim in a tumour.
* ``checkpoint_to_recruit_in_desert`` - checkpoint blockade offered to bring
  immune cells into a cold immune desert.
* ``no_myeloid_target`` - a myeloid-dominant tissue whose therapy names no
  myeloid target.
* ``pair_not_in_omnipath`` - a named ligand-receptor pair ("X-Y axis", "X/Y")
  with no interaction in the OmniPath snapshot (``scripts/fetch_omnipath.py``).

Not covered: receptors placed on the wrong cell (CD40 on epithelium; database
expression annotations are too permissive to flag it) and pathways not grounded
in the sample's data, which need per-sample ligand-receptor evidence in the prompt.
"""

import csv
import functools
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .faithfulness import split_sentences, unsupported_cell_types

KNOWLEDGE = Path(__file__).resolve().parents[2] / "data" / "knowledge"
SNAPSHOT = KNOWLEDGE / "omnipath_lr.tsv.gz"
ROLES = KNOWLEDGE / "omnipath_roles.tsv.gz"

_I = re.IGNORECASE

ALIASES: Dict[str, Tuple[str, ...]] = {
    "PD1": ("PDCD1",), "PDL1": ("CD274",), "PDL2": ("PDCD1LG2",), "B7H1": ("CD274",),
    "CTLA4": ("CTLA4",), "SIRPA": ("SIRPA",), "CD40L": ("CD40LG",), "CD154": ("CD40LG",),
    "CD40LIGAND": ("CD40LG",), "MCSF": ("CSF1",), "CSF1": ("CSF1",), "MCSFR": ("CSF1R",),
    "MCSFRECEPTOR": ("CSF1R",), "CSF1RECEPTOR": ("CSF1R",), "NKG2D": ("KLRK1",),
    "MICAB": ("MICA", "MICB"), "TIM3": ("HAVCR2",), "TGFB": ("TGFB1", "TGFB2", "TGFB3"),
    "TGFBR": ("TGFBR1", "TGFBR2"), "TGFBRECEPTOR": ("TGFBR1", "TGFBR2"),
    "PDGF": ("PDGFA", "PDGFB", "PDGFC", "PDGFD"), "PDGFR": ("PDGFRA", "PDGFRB"),
    "VEGF": ("VEGFA",), "VEGFR": ("KDR", "FLT1"), "TNFA": ("TNF",), "IFNG": ("IFNG",),
    "IL1B": ("IL1B",), "NOTCH": ("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"),
    "CXCL12": ("CXCL12",), "CXCR4": ("CXCR4",),
}
NOT_GENES = {"CAF", "CAFS", "TME", "ECM", "NK", "DC", "DCS", "TCR", "BCR", "ICI", "ICIS", "TAM", "TAMS",
             "TLS", "MOA", "HE", "RNA", "DNA", "AND", "THE", "MDSC", "MDSCS", "CD", "IL", "II", "III"}

GENE_TOKEN = re.compile(
    r"\b(?:M-CSF|CSF-?1|TGF-?β|CD40)\s+(?:ligand|receptor)\b"
    r"|\bMICA/B\b"
    r"|(?<![\w-])[A-Z][A-Za-z0-9αβγ]*(?:[-–][A-Za-z0-9αβγ]+)*(?![\w-])"
)
HYPHEN = re.compile(r"[-–]")
PAIR_JOIN = re.compile(r"^\s*[-–—/]\s*$")

INHIBITORY_GENES = {"PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD47", "SIRPA", "TGFB1", "TGFB2", "TGFB3",
                    "TGFBR1", "TGFBR2", "IL10", "IL10RA", "LAG3", "HAVCR2", "TIGIT"}
IMMUNE_GAIN = re.compile(
    r"\b(?:activat\w*|recruit\w*|stimulat\w*|infiltrat\w*)\s+of\s+(?:the\s+)?"
    r"(?:immune|T[\s-]cells?|NK|effector|lymphocytes?|cytotoxic)"
    r"|\b(?:immune|T[\s-]cell|NK[\s-]cell|lymphocyte|effector)\s+(?:cell\s+)?"
    r"(?:activation|recruitment|infiltration|stimulation|priming)\b")
VALIDATION = re.compile(r"\bvalidat\w*|\bIHC\b|immunohistochem\w*|flow cytometry|\bassays?\b", _I)
BLOCKADE = re.compile(r"\banti-|inhibit\w*|block\w*|antagon\w*|\btarget\w*|disrupt\w*|neutrali\w*"
                      r"|suppress\w*|modulat\w*", _I)

IMMUNOSUPPRESSANTS = re.compile(
    r"\b(?:infliximab|adalimumab|etanercept|golimumab|certolizumab|abatacept|belatacept"
    r"|tacrolimus|ciclosporin|cyclosporine|sirolimus|everolimus|mycophenolate|azathioprine"
    r"|corticosteroids?|glucocorticoids?|prednisone|prednisolone|dexamethasone|methylprednisolone)\b"
    r"|\banti-TNF\b|\bTNF(?:-?α)?\s+(?:inhibitor|blocker|antagonist)s?\b", _I)
SUPPRESSIVE_AIM = re.compile(
    r"\bimmunosuppressive\s+(?:agents?|drugs?|therap\w*)"
    r"|\binhibit\w*\s+T[\s-]cell\s+activation"
    r"|\bpromot\w*\s+T[\s-]cell\s+exhaustion"
    r"|\benhanc\w*\s+(?:the\s+)?(?:activity\s+of\s+)?regulatory\s+T"
    r"|\breduc\w*\s+immune\s+(?:cell\s+)?infiltration", _I)

CHECKPOINT = re.compile(r"checkpoint|\bPD-?L?1\b|CTLA-?4", _I)
RECRUIT_AIM = re.compile(r"(?:enhanc|increas|promot|improv|driv)\w*\s+(?:immune|T[\s-]cell)\s+(?:cell\s+)?"
                         r"(?:infiltration|recruitment|trafficking)", _I)
MYELOID_TARGET = re.compile(r"CSF-?1R?|M-CSF|CCR2|CCL2|CD40|CD47|SIRP|TREM2|CD11b|ITGAM|PI3Kγ|"
                            r"macrophage[\s-](?:targeting|depleting|reprogramming)|\bTAM\b|repolari", _I)

THERAPY = re.compile(r"THERAPEUTIC IMPLICATIONS\**:?(.*?)(?=\n\s*\*\*[A-Z ]+\*\*|\Z)", re.S)


@functools.lru_cache(maxsize=1)
def omnipath_pairs(path: str = str(SNAPSHOT)) -> Tuple[frozenset, frozenset]:
    """(undirected ligand-receptor pairs, gene symbols) from the OmniPath snapshot."""
    pairs, genes = set(), set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pairs.add(frozenset((row["a"], row["b"])))
            genes.update((row["a"], row["b"]))
    return frozenset(pairs), frozenset(genes)


@functools.lru_cache(maxsize=1)
def omnipath_roles(path: str = str(ROLES)) -> Dict[str, frozenset]:
    """Gene symbol -> {'ligand', 'receptor'} from the OmniPath intercell snapshot."""
    with gzip.open(path, "rt") as fh:
        return {r["gene"]: frozenset(k for k in ("ligand", "receptor") if r[k] == "1")
                for r in csv.DictReader(fh, delimiter="\t")}


def _one_sided(a: Tuple[str, ...], b: Tuple[str, ...], roles) -> bool:
    """A list, not a pair: both sides only ligands or only receptors ('CCL2/CCL5'), or co-ligands of one
    receptor ('CXCL9/CXCL10', both bind CXCR3)."""
    sides = [set().union(*(roles.get(g, frozenset()) for g in side)) for side in (a, b)]
    if any(sides[0] == sides[1] == {r} for r in ("ligand", "receptor")):
        return True
    partners = omnipath_partners()
    return any(partners.get(x, set()) & partners.get(y, set()) for x in a for y in b)


@functools.lru_cache(maxsize=1)
def omnipath_partners() -> Dict[str, Set[str]]:
    """Gene symbol -> its interaction partners in the snapshot."""
    out: Dict[str, Set[str]] = {}
    for pair in omnipath_pairs()[0]:
        x, y = tuple(pair)
        out.setdefault(x, set()).add(y)
        out.setdefault(y, set()).add(x)
    return out


def _norm(token: str) -> str:
    return re.sub(r"[\s\-/]", "", token.upper().replace("Α", "A").replace("Β", "B").replace("Γ", "G"))


def _symbols(token: str, genes) -> Optional[Tuple[str, ...]]:
    key = _norm(token)
    if key in NOT_GENES:
        return None
    return ALIASES.get(key) or ((key,) if key in genes else None)


def _split_mentions(token: str, start: int, genes) -> List[Tuple[int, int, Tuple[str, ...]]]:
    """Resolve a hyphenated token as one name ('PD-L1') or as joined names ('CD47-SIRPα'), longest first."""
    whole = _symbols(token, genes)
    if whole:
        return [(start, start + len(token), whole)]
    cuts = [m.start() for m in HYPHEN.finditer(token)]
    for cut in reversed(cuts):
        left = _symbols(token[:cut], genes)
        if left:
            return [(start, start + cut, left)] + _split_mentions(token[cut + 1:], start + cut + 1, genes)
    return []


def gene_mentions(sentence: str, genes: Iterable[str]) -> List[Tuple[int, int, Tuple[str, ...]]]:
    """(start, end, symbols) for each gene or protein named in a sentence."""
    out = []
    for m in GENE_TOKEN.finditer(sentence):
        out.extend(_split_mentions(m.group(0), m.start(), genes))
    return out


def named_pairs(sentence: str, genes: Iterable[str]) -> List[Tuple[Tuple[str, ...], Tuple[str, ...], str]]:
    """Adjacent mentions joined by '-', '–' or '/' ("CD47-SIRPα", "PD-1/PD-L1")."""
    ms = gene_mentions(sentence, genes)
    return [(a[2], b[2], sentence[a[0]:b[1]]) for a, b in zip(ms, ms[1:])
            if PAIR_JOIN.match(sentence[a[1]:b[0]]) and set(a[2]) != set(b[2])]


def therapy_section(text: str) -> str:
    """Body of the THERAPEUTIC IMPLICATIONS section ('' if absent)."""
    m = THERAPY.search(text)
    return m.group(1) if m else ""


def _finding(check: str, sentence: str, detail: str) -> Dict[str, str]:
    return {"check": check, "sentence": sentence.strip(), "detail": detail}


def check_sentence(sentence: str, prompt: Optional[str], genes: Set[str], pairs,
                   in_therapy: bool = False) -> List[Dict[str, str]]:
    """Sentence-level findings."""
    out = []
    if prompt is not None:
        extra = unsupported_cell_types(sentence, prompt)
        if extra:
            out.append(_finding("unannotated_cell_type", sentence, ", ".join(extra)))
    mentioned = {s for _, _, syms in gene_mentions(sentence, genes) for s in syms}
    gain = IMMUNE_GAIN.search(sentence)
    if mentioned & INHIBITORY_GENES and gain and not BLOCKADE.search(sentence) and not VALIDATION.search(sentence):
        out.append(_finding("inhibitory_axis_as_activating", sentence,
                            f"{'/'.join(sorted(mentioned & INHIBITORY_GENES))} with '{gain.group(0)}'"))
    m = IMMUNOSUPPRESSANTS.search(sentence) or (in_therapy and SUPPRESSIVE_AIM.search(sentence))
    if m:
        out.append(_finding("immunosuppressive_therapy", sentence, m.group(0)))
    roles = omnipath_roles()
    for a, b, span in named_pairs(sentence, genes):
        if not _one_sided(a, b, roles) and not any(frozenset((x, y)) in pairs for x in a for y in b):
            out.append(_finding("pair_not_in_omnipath", sentence, span))
    return out


def check_report(text: str, prompt: Optional[str] = None, pattern: Optional[str] = None) -> List[Dict[str, str]]:
    """All findings for one report; ``pattern`` is the dataset's spatial pattern label when known."""
    pairs, genes = omnipath_pairs()
    therapy = therapy_section(text)
    therapy_sentences = split_sentences(therapy)
    out = [f for s in split_sentences(text)
           for f in check_sentence(s, prompt, genes, pairs, in_therapy=any(t in s for t in therapy_sentences))]
    if pattern == "cold_immune_desert" and CHECKPOINT.search(therapy):
        for s in split_sentences(therapy):
            aim = RECRUIT_AIM.search(s)
            if aim:
                out.append(_finding("checkpoint_to_recruit_in_desert", s, aim.group(0)))
    if pattern == "myeloid_dominant" and therapy and not MYELOID_TARGET.search(therapy):
        out.append(_finding("no_myeloid_target", therapy.strip()[:200], "therapy names no myeloid target"))
    return out

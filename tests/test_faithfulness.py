"""Grounding checks, target scrubbing and PubMed citation verification (docs/REVIEW.md L9-L10)."""

import json
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from src.report_generation.faithfulness import (  # noqa: E402
    SECTIONS,
    extract_citations,
    scrub_report,
    split_sentences,
    unsupported_cell_types,
    unsupported_numbers,
)
from src.report_generation.prompt_builder import generate_medgemma_prompt  # noqa: E402
from src.report_generation.pubmed_check import (  # noqa: E402
    PubMedChecker,
    review_report,
    strip_unverified,
)

FEATURES = {
    "annotation": {"cell_type_counts": {"Epithelial": 4438, "Stromal_CAF": 2118, "T_cells": 2805,
                                        "Macrophages": 997, "Unknown": 428},
                   "mean_confidence": 0.864},
    "spatial_heterogeneity": {"morans_i_mean": 0.5639, "entropy_mean": 1.7033, "n_enriched_pairs": 141},
    "uncertainty": {"mean_prediction_entropy": 0.8128},
}
PROMPT = generate_medgemma_prompt(FEATURES)

TARGET = """**SPATIAL TRANSCRIPTOMICS PATHOLOGY REPORT**

**SPECIMEN ANALYSIS:** Visium HD Spatial Profiling

**TISSUE MICROENVIRONMENT CHARACTERIZATION:**

The tissue shows strong spatial autocorrelation (Moran's I: 0.56), indicating organised compartments.

**MECHANISM OF ACTION:**

Macrophages and T cells line the epithelial-stromal interface via CXCL9/CXCL10 signalling. This pattern is consistent with organized anti-tumor immunity, as described by Helmink et al. (Nature Medicine, 2020), who linked B cells to response. Plasma cells support humoral immunity.

**THERAPEUTIC RECOMMENDATION:**

Consider anti-PD-1 with TGF-β inhibition to relieve CAF-mediated IL-10 suppression.

**QUALITY ASSURANCE:**

Mean deconvolution confidence (86.4%) is acceptable."""


def test_split_keeps_et_al_citation_in_one_sentence():
    sents = split_sentences("Pattern seen before, as described by Helmink et al. (Nature, 2020), in melanoma. Next one.")
    assert len(sents) == 2 and "Helmink et al. (Nature, 2020)" in sents[0]


def test_gene_and_marker_names_are_not_numbers():
    text = "CXCL9, IL-10, PD-1, CD8+ and Th17 signalling with 4-1BB agonists in 10x data."
    assert unsupported_numbers(text, PROMPT) == []
    assert unsupported_numbers("About 62% of spots and 86.4% confidence.", PROMPT) == ["62%", "86.4%"]
    assert unsupported_numbers("Moran's I of 0.56 and 141 pairs.", PROMPT) == []


def test_cell_types_must_be_listed_in_prompt():
    assert "Macrophages" in PROMPT and "Stromal CAF" in PROMPT and "Unknown" not in PROMPT
    assert unsupported_cell_types("Macrophages, CAFs and T cells", PROMPT) == []
    assert unsupported_cell_types("B cells, Tregs and MDSCs", PROMPT) == ["b_cell", "mdsc", "treg"]


def test_extract_citation_fields():
    (c,) = extract_citations("as described by Helmink et al. (Nature Medicine, 2020), who")
    assert (c["surname"], c["year"], c["journal"]) == ("Helmink", 2020, "Nature Medicine")
    (c,) = extract_citations("(Kieffer et al., 2020)")
    assert (c["surname"], c["year"], c["journal"]) == ("Kieffer", 2020, None)


def test_scrub_removes_ungrounded_sentences_and_keeps_headings():
    clean, removed = scrub_report(TARGET, PROMPT)
    reasons = {r["sentence"][:20]: r["issues"] for r in removed}
    assert any("citation" in v for v in reasons.values())
    assert ["number"] in reasons.values() and ["cell_type"] in reasons.values()
    assert "Helmink" not in clean and "86.4%" not in clean and "Plasma" not in clean
    assert "CXCL9/CXCL10" in clean and "IL-10" in clean and "0.56" in clean
    assert clean.startswith("**TISSUE MICROENVIRONMENT**")
    assert [s for s in SECTIONS if f"**{s}**" in clean] == list(SECTIONS[:3])
    assert "PATHOLOGY REPORT" not in clean and "Visium HD" not in clean


def test_prompt_and_scrubbed_target_agree_on_sections():
    for section in SECTIONS:
        assert f"**{section}**" in PROMPT
    assert "Do not cite specific studies" in PROMPT


class FakeEutils:
    """Minimal esearch/esummary stand-in keyed by query term."""

    def __init__(self, hits):
        self.hits, self.calls = hits, []

    def __call__(self, url, timeout):
        q = parse_qs(urlparse(url).query)
        self.calls.append(q.get("term", q.get("id"))[0])
        if "esummary" in url:
            pmid = q["id"][0]
            return {"result": {pmid: {"title": "B cells and TLS", "fulljournalname": "Nature",
                                      "pubdate": "2020 Jan"}}}
        return {"esearchresult": {"idlist": self.hits.get(q["term"][0], [])}}


HITS = {
    'Helmink[1au] AND 2020[pdat] AND "Nature"[jour]': ["31942075"],
    "Helmink[1au] AND 2020[pdat]": ["31942075"],
    "Kieffer[1au] AND 2020[pdat]": ["32434947", "33301138"],
}


def test_pubmed_statuses():
    checker = PubMedChecker(fetch=FakeEutils(HITS))
    checker.min_interval = 0
    assert checker.verify("Helmink", 2020, "Nature")["status"] == "verified"
    assert checker.verify("Helmink", 2020, "Nature")["pmid"] == "31942075"
    assert checker.verify("Helmink", 2020, "Nature Medicine")["status"] == "mismatch"
    assert checker.verify("Kieffer", 2020, None)["status"] == "incomplete"
    assert checker.verify("Nobody", 2019, "Cell")["status"] == "not_found"


def test_pubmed_offline_degrades_to_unchecked(tmp_path):
    def boom(url, timeout):
        raise OSError("network down")

    checker = PubMedChecker(cache_path=str(tmp_path / "c.json"), fetch=boom)
    checker.min_interval = 0
    assert checker.verify("Helmink", 2020, "Nature")["status"] == "unchecked"


def test_pubmed_cache_avoids_repeat_requests(tmp_path):
    fake = FakeEutils(HITS)
    first = PubMedChecker(cache_path=str(tmp_path / "c.json"), fetch=fake)
    first.min_interval = 0
    first.verify("Helmink", 2020, "Nature")
    first.save()
    n = len(fake.calls)
    second = PubMedChecker(cache_path=str(tmp_path / "c.json"), fetch=fake)
    second.verify("Helmink", 2020, "Nature")
    assert len(fake.calls) == n


def test_review_strips_only_unconfirmed_citations():
    text = ("**MECHANISM OF ACTION**\n\nT cells infiltrate. Helmink et al. (Nature, 2020) linked TLS to "
            "response. Helmink et al. (Nature Medicine, 2020) said the same. Macrophages persist.")
    checker = PubMedChecker(fetch=FakeEutils(HITS))
    checker.min_interval = 0
    review = review_report(text, PROMPT, checker=checker)
    statuses = [c["status"] for c in review["pubmed"]]
    assert statuses == ["verified", "mismatch"]
    assert "(Nature, 2020)" in review["clean_text"]
    assert "Nature Medicine" not in review["clean_text"]
    assert "Macrophages persist." in review["clean_text"]
    assert review_report(text, PROMPT, verify_citations=False)["clean_text"] == text


def test_review_removes_ungrounded_cell_types_and_numbers():
    text = ("**MECHANISM OF ACTION**\n\nCAFs exclude T cells via TGF-β. MDSCs suppress immunity. "
            "Roughly 62% of spots are stromal. Moran's I of 0.56 supports compartments.")
    review = review_report(text, PROMPT)
    assert review["unsupported_cell_types"] == ["mdsc"] and review["unsupported_numbers"] == ["62%"]
    assert review["clean_text"] == ("**MECHANISM OF ACTION**\n\nCAFs exclude T cells via TGF-β. "
                                    "Moran's I of 0.56 supports compartments.")


def test_strip_unverified_noop_without_bad_citations():
    assert strip_unverified("A. B.", [{"status": "verified", "text": "x"}]) == "A. B."


@pytest.mark.network
def test_live_pubmed_distinguishes_real_and_miscited_papers():
    checker = PubMedChecker()
    cases = {("Helmink", 2020, "Nature"): "verified", ("Helmink", 2020, "Nature Medicine"): "mismatch",
             ("Kieffer", 2020, "Cancer Discovery"): "verified", ("Kieffer", 2020, "Cancer Cell"): "mismatch"}
    got = {case: checker.verify(*case) for case in cases}
    if any(r["status"] == "unchecked" for r in got.values()):
        pytest.skip("PubMed unreachable or rate limited")
    assert {case: r["status"] for case, r in got.items()} == cases
    assert got[("Helmink", 2020, "Nature")]["pmid"] == "31942075"


def test_pubmed_retries_rate_limit(monkeypatch):
    import urllib.error

    fake, state = FakeEutils(HITS), {"n": 0}

    def flaky(url, timeout):
        state["n"] += 1
        if state["n"] == 1:
            raise urllib.error.HTTPError(url, 429, "Too Many Requests", None, None)
        return fake(url, timeout)

    monkeypatch.setattr("time.sleep", lambda s: None)
    checker = PubMedChecker(fetch=flaky)
    assert checker.verify("Helmink", 2020, "Nature")["status"] == "verified"


def test_multimodal_prompt_grounds_listed_cell_types():
    from src.report_generation.medgemma_multimodal import create_multimodal_prompt

    prompt = create_multimodal_prompt(FEATURES)
    text = "**TISSUE MICROENVIRONMENT**\n\nT cells and macrophages border CAF-rich stroma around epithelial nests."
    review = review_report(text, prompt, verify_citations=False)
    assert review["unsupported_cell_types"] == [] and review["clean_text"] == text


def test_section_emptied_by_guard_loses_its_heading():
    text = "**MECHANISM OF ACTION**\n\nMDSCs suppress immunity.\n\n**VALIDATION CAVEATS**\n\nConfirm with IHC."
    assert review_report(text, PROMPT)["clean_text"] == "**VALIDATION CAVEATS**\n\nConfirm with IHC."


def test_atlas_labels_ground_their_lineage():
    from src.report_generation.faithfulness import label_lineage

    assert label_lineage("LummHR-SCGB") == "epithelial cells" and label_lineage("Fibro-major") == "fibroblasts"
    assert label_lineage("CD4-Treg") == "regulatory T cells" and label_lineage("Lymph-major") == "endothelial cells"
    assert label_lineage("T_cells") == "" and label_lineage("Lymphocytes") == ""
    features = {"annotation": {"cell_type_counts": {"LummHR-SCGB": 3000, "plasma_IgG": 400, "Macro-m2": 100}}}
    prompt = generate_medgemma_prompt(features)
    assert "LummHR-SCGB (epithelial cells)" in prompt and "plasma IgG (plasma cells)" in prompt
    text = "Luminal epithelial cells, plasma cells and M2 macrophages dominate; T cells are sparse."
    assert unsupported_cell_types(text, prompt) == ["t_cell"]


def test_multimodal_prompt_matches_training_prompt_plus_image_note():
    from src.report_generation.medgemma_multimodal import create_multimodal_prompt

    prompt = create_multimodal_prompt(FEATURES)
    assert prompt.replace(prompt[prompt.index("H&E IMAGE"):prompt.index("SPATIAL TRANSCRIPTOMICS")], "") == PROMPT
    assert all(f"**{s}**" in prompt for s in SECTIONS)


def test_copied_spans_flag_verbatim_prompt_reuse():
    from src.report_generation.faithfulness import copied_spans

    prompt = "Name the mechanism (e.g., TGF-β-mediated exclusion, IFN-γ/CXCL9 chemokine gradient)."
    text = "**MECHANISM OF ACTION**\n\nAn IFN-γ/CXCL9 chemokine gradient recruits T cells. CAFs restrict access."
    assert copied_spans(text, prompt, n=3) == ["ifn-γ/cxcl9 chemokine gradient"]
    assert copied_spans("**MECHANISM OF ACTION**\n\nT cells cluster near CAFs.", prompt, n=3) == []
    assert "e.g." not in PROMPT and "TGF-β-mediated exclusion" not in PROMPT


def test_moa_novel_ignores_copied_terms():
    sys.path.insert(0, str(ROOT / "scripts"))
    from compare_adapter import score

    prompt = "Explain the TGF-beta signalling pathway mediated by fibroblasts in this sample please."
    copied = "**MECHANISM OF ACTION**\n\nWe see the TGF-beta signalling pathway mediated by fibroblasts here."
    rec = score(copied, prompt)
    assert rec["copied_spans"] and rec["moa_novel"] < rec["moa_hits"]


def test_check_report_cli_offline(tmp_path, capsys):
    from check_report import main

    report = tmp_path / "r.md"
    report.write_text("**MECHANISM OF ACTION**\n\nMDSCs suppress immunity. Helmink et al. (Nature, 2020) agree. "
                      "T cells border CAF-rich stroma.")
    prompt = tmp_path / "p.txt"
    prompt.write_text(PROMPT)
    main(["--report", str(report), "--prompt", str(prompt), "--offline", "--json"])
    out = json.loads(capsys.readouterr().out)
    assert out["unsupported_cell_types"] == ["mdsc"] and out["pubmed"][0]["status"] == "unchecked"
    assert out["clean_text"].endswith("T cells border CAF-rich stroma.") and out["changed"]


def test_strip_thinking_keeps_only_the_answer():
    from src.report_generation.faithfulness import strip_thinking

    assert strip_thinking("<unused94>thought\nPlan: copy the prompt.<unused95>**TISSUE MICROENVIRONMENT**\nX.") \
        == "**TISSUE MICROENVIRONMENT**\nX."
    assert strip_thinking("<unused94>thought\nstill planning when tokens ran out") == ""
    assert strip_thinking("**VALIDATION CAVEATS**\nConfirm.") == "**VALIDATION CAVEATS**\nConfirm."


def test_normalise_headings_fixes_misspelled_sections_only():
    from src.report_generation.faithfulness import normalise_headings

    text = "**VALIDATION CAVES**\nConfirm by IHC.\n\n**MECHANISM OF ACTION**\nX.\n\n**Notes**\nY."
    assert normalise_headings(text) == "**VALIDATION CAVEATS**\nConfirm by IHC.\n\n**MECHANISM OF ACTION**\nX.\n\n**Notes**\nY."
    assert normalise_headings("**VALIDATION CAVESATS**: check it") == "**VALIDATION CAVEATS** check it"

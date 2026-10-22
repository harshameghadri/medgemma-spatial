"""Biology checks against the blinded manual read (docs/REVIEW.md L18)."""

import json
from pathlib import Path

import pytest

from src.report_generation.biology_check import check_report, gene_mentions, named_pairs, omnipath_pairs

READS = [json.loads(line) for line in (Path(__file__).parent / "fixtures/biology_read/reads.jsonl").open()]


def checks(text, prompt=None, pattern=None):
    return {f["check"] for f in check_report(text, prompt, pattern)}


def test_snapshot_has_the_axes_the_reviewer_named():
    pairs, _ = omnipath_pairs()
    for a, b in [("CD47", "SIRPA"), ("CD274", "PDCD1"), ("CD40LG", "CD40"), ("CSF1", "CSF1R"), ("MICA", "KLRK1")]:
        assert frozenset((a, b)) in pairs


def test_hyphenated_names_resolve_to_one_gene_or_a_pair():
    _, genes = omnipath_pairs()
    assert [m[2] for m in gene_mentions("PD-L1 and TGF-β", genes)] == [("CD274",), ("TGFB1", "TGFB2", "TGFB3")]
    assert [(a, b) for a, b, _ in named_pairs("the CD47-SIRPα axis", genes)] == [(("CD47",), ("SIRPA",))]


@pytest.mark.parametrize("sentence, flagged", [
    ("This involves the activation of T cells by CAFs through the PD-1/PD-L1 axis.", True),
    ("CAFs may express PD-L1, which interacts with PD-1 on T cells, suppressing their activity.", False),
    ("Anti-PD-1 therapy restores the activation of T cells.", False),
    ("The TGF-β pathway promotes fibrosis and immune cell recruitment.", True),
    ("Validate PD-L1 by IHC and assess the activation of T cells by flow cytometry.", False),
])
def test_inhibitory_axis_described_as_activating(sentence, flagged):
    assert ("inhibitory_axis_as_activating" in checks(sentence)) is flagged


def test_ligand_lists_are_not_pairs_but_wrong_pairs_are():
    assert "pair_not_in_omnipath" not in checks("A role for the CCL2/CCL5 chemokine axis.")
    assert "pair_not_in_omnipath" in checks("Mediated by the CCL2/CX3CR1 axis.")
    assert "pair_not_in_omnipath" not in checks("Mediated by the CCL2/CCR2 axis.")
    assert "pair_not_in_omnipath" not in checks("Effector T cells are recruited via CXCL9/CXCL10 and CXCR3.")
    assert "pair_not_in_omnipath" not in checks("Dendritic cells present CD80/CD86 to CTLA-4.")


def test_suppressive_aim_counts_only_in_the_therapy_section():
    mech = "**MECHANISM OF ACTION** CAFs upregulate PD-L1 to inhibit T cell activation."
    ther = "**THERAPEUTIC IMPLICATIONS** Inhibiting T cell activation might be a viable approach."
    assert "immunosuppressive_therapy" not in checks(mech)
    assert "immunosuppressive_therapy" in checks(ther)
    assert "immunosuppressive_therapy" in checks("A TNF inhibitor (e.g., infliximab) may help.")


def test_phenotype_rules_need_the_pattern():
    ther = "**THERAPEUTIC IMPLICATIONS** Anti-PD-1 and anti-CTLA-4 may enhance immune cell infiltration."
    assert "checkpoint_to_recruit_in_desert" in checks(ther, pattern="cold_immune_desert")
    assert "no_myeloid_target" in checks(ther, pattern="myeloid_dominant")
    assert not checks(ther, pattern="hot_immune")
    assert "no_myeloid_target" not in checks(ther + " A CSF1R inhibitor targets TAMs.", pattern="myeloid_dominant")


def test_no_flags_on_reports_the_reviewer_passed():
    clean = [r for r in READS if r["errors"] == 0]
    assert len(clean) == 13
    assert {r["case"] + r["model"]: checks(r["text"], r["prompt"], r["pattern"]) for r in clean} == \
        {r["case"] + r["model"]: set() for r in clean}


def test_flags_most_reports_the_reviewer_failed():
    bad = [r for r in READS if r["errors"] > 0]
    caught = [r for r in bad if checks(r["text"], r["prompt"], r["pattern"])]
    assert len(bad) == 15 and len(caught) >= 11

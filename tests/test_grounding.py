"""Inference-time grounding: curated axes, retrieval corpus and prompt assembly."""

import json
from pathlib import Path

from src.report_generation.biology_check import check_report, omnipath_pairs
from src.report_generation.faithfulness import CELL_TYPE_PATTERNS
from src.report_generation.grounding import (contexts, ground_prompt, load_approved, load_axes, parse_summary,
                                             retrieve, select_axes)

READS = [json.loads(line) for line in (Path(__file__).parent / "fixtures/biology_read/reads.jsonl").open()]
EVAL_PROMPTS = {r["case"]: r["prompt"] for r in READS}


def test_every_curated_axis_is_an_omnipath_pair_in_the_stated_direction():
    pairs, _ = omnipath_pairs()
    for a in load_axes():
        assert frozenset((a["ligand"], a["receptor"])) in pairs, a["axis"]
        assert set(a["sender"]) | set(a["receiver"]) <= set(CELL_TYPE_PATTERNS), a["axis"]


def test_cd40_is_received_by_antigen_presenting_cells_not_epithelium():
    cd40 = next(a for a in load_axes() if a["receptor"] == "CD40")
    assert "epithelial" not in cd40["receiver"] and "dendritic_cell" in cd40["receiver"]


def test_axes_need_both_cell_types_present():
    summary = {"phenotype": "immune-rich (high immune density)", "cell_types": ["Epithelial", "T cells"],
               "morans_i": 0.3, "entropy": 1.0}
    axes = {a["axis"] for a in select_axes(summary)}
    assert "PD-L1–PD-1" in axes
    assert "CSF1–CSF1R" not in axes


def test_myeloid_samples_rank_myeloid_axes_first():
    s = parse_summary(EVAL_PROMPTS["p7"])
    assert "myeloid" in contexts(s)
    top = [a["axis"] for a in select_axes(s)][:3]
    assert {"CD47–SIRPα", "CSF1–CSF1R"} & set(top)


def test_ground_prompt_inserts_before_the_instructions_and_keeps_the_data():
    source = EVAL_PROMPTS["p1"]
    grounded, added = ground_prompt(source)
    assert added and added in grounded
    assert grounded.index("CANDIDATE LIGAND-RECEPTOR AXES") < grounded.index("Generate a ")
    assert grounded.replace(added, "") == source
    assert ground_prompt("no summary here") == ("no summary here", "")


def test_corpus_is_clean_and_never_contains_an_eval_prompt():
    corpus = load_approved()
    assert len(corpus) == 131
    eval_summaries = [parse_summary(p) for p in EVAL_PROMPTS.values()]
    assert not [r["id"] for r in corpus if r["summary"] in eval_summaries]
    assert all(r["id"][0] == "r" and ".v" in r["id"] for r in corpus)
    assert not [r["id"] for r in corpus if check_report(r["text"], None, r["pattern"])]


def test_retrieval_prefers_the_same_phenotype():
    s = parse_summary(EVAL_PROMPTS["p6"])
    assert all(r["summary"]["phenotype"] == s["phenotype"] for r in retrieve(s))


def test_serving_respects_the_switch_and_infers_the_pattern(monkeypatch):
    from src.report_generation.grounding import infer_pattern, serve_prompt

    source = EVAL_PROMPTS["p7"]
    monkeypatch.setenv("MEDGEMMA_GROUNDING", "0")
    assert serve_prompt(source) == (source, {"grounding": False})
    monkeypatch.setenv("MEDGEMMA_GROUNDING", "1")
    served, meta = serve_prompt(source)
    assert meta["grounding"] and served != source
    assert infer_pattern(source) == "myeloid_dominant"


def test_multimodal_prompt_is_grounded_too():
    from src.report_generation.grounding import ground_prompt
    from src.report_generation.medgemma_multimodal import create_multimodal_prompt
    from tests.test_faithfulness import FEATURES

    _, added = ground_prompt(create_multimodal_prompt(FEATURES))
    assert "CANDIDATE LIGAND-RECEPTOR AXES" in added

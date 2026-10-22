"""LoRA data + training guards (docs/REVIEW.md L1-L7)."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_lora_dataset import rebuild, write_jsonl  # noqa: E402
from generate_training_data import (  # noqa: E402
    build_example_list,
    build_user_prompt,
    compute_split,
    format_chatml,
)

from src.report_generation.prompt_builder import generate_medgemma_prompt  # noqa: E402

TINY_GEMMA3 = "trl-internal-testing/tiny-Gemma3ForConditionalGeneration"
COMPLETION = "TISSUE MICROENVIRONMENT: synthetic completion used only for tests. " * 12


def legacy_records(n):
    """Records in the original single-'text' format, built from the seeded features."""
    return [format_chatml(build_user_prompt(f), COMPLETION) for _, f in build_example_list(n, seed=42)]


def test_split_matches_published_dataset():
    assert compute_split(240) == 206


def test_rebuild_uses_serving_prompt():
    rows = rebuild(legacy_records(7))
    features = build_example_list(7, seed=42)
    for row, (pattern, feats) in zip(rows, features):
        assert row["pattern"] == pattern
        assert row["prompt"][-1]["content"] == generate_medgemma_prompt(feats)
        assert row["completion"] == [{"role": "assistant", "content": COMPLETION.strip()}]


def test_rebuild_tolerates_last_ulp_float_differences():
    records = legacy_records(3)
    _, feats = build_example_list(3, seed=42)[2]
    value = feats["spatial_heterogeneity"]["entropy_mean"]
    nudged = float(np.nextafter(value, 0))
    records[2]["text"] = records[2]["text"].replace(repr(value), repr(nudged))
    rows = rebuild(records)
    assert f"({nudged:.2f})" in rows[2]["prompt"][-1]["content"]


def test_rebuild_rejects_records_not_from_seed():
    records = legacy_records(3)
    records[1]["text"] = records[1]["text"].replace('"n_spots"', '"n_spotz"')
    with pytest.raises(ValueError, match="Record 1"):
        rebuild(records)


@pytest.fixture(scope="module")
def tiny_dataset(tmp_path_factory):
    rows = rebuild(legacy_records(6))
    d = tmp_path_factory.mktemp("lora")
    write_jsonl(d / "train.jsonl", rows[:4])
    write_jsonl(d / "eval.jsonl", rows[4:])
    return d


def test_legacy_text_format_is_refused(tmp_path):
    pytest.importorskip("datasets")
    import train_lora

    write_jsonl(tmp_path / "train.jsonl", legacy_records(2))
    with pytest.raises(ValueError, match="prompt"):
        train_lora.load_data(str(tmp_path))


@pytest.mark.network
def test_length_guard_rejects_truncating_limits(tiny_dataset):
    pytest.importorskip("trl")
    import train_lora
    from transformers import AutoTokenizer

    tok = AutoTokenizer.from_pretrained(TINY_GEMMA3)
    ds = train_lora.load_data(str(tiny_dataset))
    with pytest.raises(ValueError, match="truncated"):
        train_lora.check_lengths(tok, ds["train"], max_length=256)
    assert train_lora.check_lengths(tok, ds["train"], max_length=1024) <= 1024


@pytest.mark.network
@pytest.mark.slow
def test_training_masks_prompt_and_saves_loadable_adapter(tiny_dataset, tmp_path):
    pytest.importorskip("trl")
    import train_lora
    from src.report_generation import medgemma_multimodal as mm

    out = tmp_path / "adapter"
    trainer = train_lora.main([
        "--model-id", TINY_GEMMA3, "--data-dir", str(tiny_dataset), "--output-dir", str(out),
        "--max-steps", "2", "--grad-accum", "1", "--quantize", "none", "--precision", "fp32",
    ])

    batch = next(iter(trainer.get_train_dataloader()))
    labels = batch["labels"][0]
    n_trained = int((labels != -100).sum())
    n_masked = int((labels == -100).sum())
    assert n_trained >= 100, "completion tokens must carry the loss"
    assert n_masked > n_trained, "prompt tokens must be masked"
    assert all(h["loss"] == h["loss"] for h in trainer.state.log_history if "loss" in h)

    assert json.loads((out / "adapter_config.json").read_text())["r"] == 16
    report, meta = mm.generate_textonly_fallback(
        {}, model_id=TINY_GEMMA3, device="cpu", max_new_tokens=4, adapter_id=str(out))
    assert meta["adapter_id"] == str(out)

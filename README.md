# MedGemma Spatial

Spatial transcriptomics in, a draft pathology report out, and an honest account of how often the biology in that
report is wrong.

[![CI](https://github.com/harshameghadri/medgemma-spatial/actions/workflows/ci.yml/badge.svg)](https://github.com/harshameghadri/medgemma-spatial/actions/workflows/ci.yml)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![License: Unlicense](https://img.shields.io/badge/license-Unlicense-blue.svg)](LICENSE)

> **Experimental, research use only.** Not a medical device. Reports are drafts for a pathologist to check, not clinical
> advice. Provided as is, with no warranty of any kind.

## What it does

```
10x Visium / Visium HD (Space Ranger output or .h5ad)
  → Scanpy + Squidpy: QC, cell-type annotation, Moran's I, neighbourhood entropy, co-enrichment
  → feature summary (cell types, phenotype, spatial metrics)
  → evidence prompt: curated ligand-receptor axes for the cell types present
                     (OmniPath-confirmed pairs, sender → receiver, effect) + therapy rules
  → MedGemma 4B (text, or text + H&E image)
  → checks: unsupported numbers / cell types, PubMed-verified citations, prompt copying,
            and a biology check for the error kinds the reviewer found
  → Streamlit app with the report, the warnings and a download
```

## What I learned building it

The interesting part of this project is not the pipeline; it is measuring whether a small model's biology can be
trusted, and what actually moves that number. Every step is logged with its evidence in
[`docs/REVIEW.md`](docs/REVIEW.md) (findings L1-L20).

1. **The original fine-tune trained on the prompt, not the report.** Examples were truncated at 256–512 tokens, so the
   targets never reached the model, and the targets themselves carried wrong citations (Helmink et al. 2020 cited
   as *Nature Medicine* 34 times; it is *Nature*). Fixed: completion-only loss, length checks, citation scrubbing.
2. **Automatic metrics were not enough.** I distilled a local Llama 3.1 8B teacher into two students, with
   gates, a 189-card human review and criteria fixed before training. The MedGemma student passed every
   automatic criterion (mechanism terms 3.1 → 7.0 per report, prompt copying 4 → 0) and then **failed the blinded
   biology read**: 4 errors against the base model's 3. The Llama student failed the automatic criteria by 0.01.
   Neither adapter ships.
3. **The errors had a shape.** The reviewer's 20 errors fell into a few kinds: inhibitory axes described as
   activating ("PD-1/PD-L1 promotes T cell activation"), immunosuppressants proposed as anti-tumour therapy,
   receptors placed on the wrong cell (CD40 on epithelium), cell types the data never showed (Tregs, MDSCs), and
   therapy that ignores the phenotype. Those labels became a checker (`biology_check.py`) and then the design of
   round 2.
4. **Retrieving prose made the model copy; retrieving facts made it specific.** Round 2 grounds the prompt instead
   of retraining. Adding the two most similar reviewer-approved reports made MedGemma retell them: 86% of each
   report's words were copied, and the raw mechanism count doubled, which would have read as a win if copied text
   had not been subtracted. Llama 3.1 8B did the same. A 23-row axis table instead (every pair confirmed in
   OmniPath, sender → receiver and effect on the receiving cell, signed off by the reviewer) worked.

### Round 2 results (MedGemma 4B, 31 held-out prompts; blinded read on 7)

| Prompt | Biology-check flags | Mechanism terms (own words) | Copied from inserted text | Blinded read: errors | Generic |
|---|---|---|---|---|---|
| Plain | 6 reports | 3.13 | — | 3 | 7 / 7 |
| + 2 retrieved approved reports | 1 | 1.23 | 86% (median) | not run (failed) | — |
| **+ curated axes and rules (served)** | **0** | **6.87** | **0%** | **2** | **0 / 7** |

The evidence prompt is the app's default. Caveats, all in [`REVIEW.md`](docs/REVIEW.md) L20-L23: the pass needed a
checker fix made after seeing the outputs (co-ligand lists were read as wrong pairs); the reviewer could tell the
conditions apart from style, so the read was not truly blind; seven cases is a small sample; and the remaining two
errors (inhibiting a favourable TLS axis; checkpoint monotherapy in stroma-trapped tissue) are the next fixes.

## Quick start

```bash
git clone https://github.com/harshameghadri/medgemma-spatial.git
cd medgemma-spatial
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pytest -q                                   # ~2 min, CPU only, synthetic tissue

streamlit run app/streamlit_app.py          # upload a .h5ad; MedGemma needs a GPU or Apple MPS

# or from the command line (Space Ranger folder with filtered_feature_bc_matrix.h5 + spatial/)
python scripts/run_sample.py analyse --visium-dir data/my_sample --out runs/my_sample
python scripts/run_sample.py report --out runs/my_sample --mode multimodal
```

MedGemma is gated on Hugging Face: accept its terms and set `HF_TOKEN`. `MEDGEMMA_GROUNDING` selects the prompt:
`evidence` (default: curated axes and therapy rules), `off` (plain) or `full` (plus retrieved reports; not
recommended, see results).

## Check any model's report

`scripts/check_report.py` needs no GPU and no model; it works on any LLM's output given the exact prompt.

```bash
python scripts/check_report.py --report report.md --prompt prompt.txt
```

It flags numbers and cell types the prompt never gave, citations that PubMed cannot confirm as written, and 5+ word
runs copied from the prompt. `biology_check.check_report` adds the six biology checks. They catch the specific
errors they were built from; they do not replace reading the report.

To compare models on the same 31 held-out prompts:

```bash
python scripts/compare_adapter.py --model-id google/medgemma-4b-it --eval data/lora_v4/eval.jsonl --out runs/x.json
python scripts/build_grounded_eval.py --eval data/lora_v4/eval.jsonl --out data/grounded_v1/eval.jsonl
```

## Repository map

| Path | What it is |
|---|---|
| `src/spatial_analysis/` | annotation, spatial statistics, uncertainty |
| `src/report_generation/prompt_builder.py` | feature summary → prompt |
| `src/report_generation/grounding.py` | axis evidence + retrieval of approved reports |
| `src/report_generation/faithfulness.py`, `pubmed_check.py` | grounding checks, PubMed citation lookup |
| `src/report_generation/biology_check.py` | the six biology checks, OmniPath pair lookup |
| `data/knowledge/` | OmniPath snapshot, curated axes, approved reports (with sources and licences) |
| `scripts/train_lora.py`, `draft_targets.py`, `build_distill_set.py` | QLoRA training and the distillation pipeline |
| `scripts/compare_adapter.py` | pre-registered model comparisons |
| `app/` | Streamlit app and Dockerfile |
| `tests/` | 120+ tests on synthetic tissue, including the blinded-read labels |
| `docs/REVIEW.md` | every finding, decision and result, with evidence |

## Data and licences

- Code: [The Unlicense](LICENSE) (public domain, provided as is). This covers the code only, not the data below.
- OmniPath snapshot: commercial-use resources only; citation and terms in `data/knowledge/README.md`.
- `data/knowledge/approved_reports.jsonl`: Llama 3.1 8B outputs picked by the reviewer. Built with Llama; the Llama
  3.1 Community License applies to these outputs.
- MedGemma: Health AI Developer Foundations terms. No fine-tuned weights are published.

## Acknowledgements

Scanpy (Wolf et al., 2018), Squidpy (Palla et al., 2022), OmniPath (Türei et al., 2016, 2021), MedGemma (Google),
Llama 3.1 (Meta), 10x Genomics public Visium datasets.

Built by Sri Harsha Meghadri.

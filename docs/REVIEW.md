# Code Review — October 2026

An objective review of the repository after the Kaggle MedGemma challenge
closed (Feb 24, 2026). The goal now is portfolio hardening: make the pipeline
correct, make CI able to fail, and make LoRA training reproducible.

Every finding below was reproduced, not inferred. Severity: **P0** = output is
wrong or the component cannot work; **P1** = will fail under realistic use;
**P2** = hygiene/compliance.

## 1. LoRA fine-tuning (why Kaggle "kept failing")

| # | Sev | Finding | Evidence |
|---|-----|---------|----------|
| L1 | P0 | **No training example contains any target text after truncation.** Prompt alone is 453–473 tokens; full examples are 877–970 tokens. Notebook uses `max_seq_length=256`; earlier runs used 512. | MedGemma tokenizer over `data/lora_training/train.jsonl`: 0/206 examples have a completion token inside 256; 0/206 fit in 512. Every past run trained the model to reproduce the input JSON. |
| L2 | P0 | **Last Kaggle run never reached training.** Cell 1 calls `kernel.do_shutdown(restart=True)`; under Kaggle's papermill runner that kills the run. | Kernel log (`kaggle kernels output`): `Restarting kernel...` → `nbclient.exceptions.DeadKernelError: Kernel died`. Status `ERROR`. |
| L3 | P1 | Unpinned `pip install -U transformers accelerate peft …` next to a pinned Unsloth release. | Same log: pip upgraded torch to 2.10 and reported `fastai requires torch<2.10` / `cuda-bindings` conflicts. |
| L4 | P1 | Loss is computed over prompt + completion (no completion-only masking), so ~50% of the gradient is spent on a templated JSON dump. | `SFTTrainer(dataset_text_field="text")`, no assistant-only mask. |
| L5 | P1 | Train/serve skew: training prompts are a raw `json.dumps(features, indent=2)`; serving prompts come from `prompt_builder.generate_medgemma_prompt`. The adapter is also never loaded by `src/` or `app/`. | `grep -rn PeftModel src app` → nothing. |
| L6 | P1 | `scripts/train_lora.py` uses removed TRL arguments (`tokenizer=`, `max_seq_length=` in `SFTConfig`) and `transformers==4.45.1` (cannot load Gemma 3, see I1). | — |
| L8 | P2 | Synthetic features are out of distribution for the real pipeline: the generator draws `n_enriched_pairs` up to 200, but with ~11 cell types there are at most 66 unordered pairs, and the cell-type vocabulary is fixed to 11 labels while the pipeline emits up to 20 panels plus CellTypist labels. | `generate_training_data.py` |
| L9 | P0 | **Adapter v1 (trained 2026-10-20, RTX 3080) is less faithful than base MedGemma; do not deploy it.** The targets were written by an LLM that saw the full feature JSON (counts, exact metrics), while the serving prompt deliberately hides those numbers, so the adapter learns to state specifics it was never given. Training itself was healthy (eval loss 2.55 → 1.445, train ≈ eval). | `scripts/compare_adapter.py` on all 34 held-out prompts. Base: 0/34 with citations, 6/34 with numbers absent from the prompt, median 216 words. Adapter: **15/34 with citations, 26/34 with unsupported numbers (44 total)**, median 301 words (18/34 within the 280–350 target), fewer mechanism terms (3.9 vs 4.8). |
| L10 | — | **Adapter v2 (2026-10-20, RTX 3080, no API calls): far more grounded than v1 but still fails one pre-set criterion (unsupported cell types). Default to base MedGemma with the new prompt; adapter v2 is opt-in via `MEDGEMMA_ADAPTER`.** Targets lose every sentence with a citation, a number absent from the prompt, or a cell type the prompt does not list (645 sentences; 233/240 targets kept, median 204 words). The prompt now lists detected cell types and annotation confidence, uses the targets' four sections, forbids citations and asks for 170-230 words. Re-scoring v1 with the corrected checks (the old splitter cut `et al.` citations in half and the old number check counted gene names such as CXCL9) gives 16/34 reports with citations, 35 unsupported numbers and 165 unsupported cell types. The new prompt alone makes base MedGemma grounded (v1's base had 22 unsupported cell types; its earlier "6 unsupported numbers" were gene names). The remaining adapter failure is MDSCs (pretraining prior, absent from the cleaned targets); the app now removes ungrounded sentences and sentences citing papers PubMed cannot confirm. Grounding checks cannot catch wrong biology, e.g. MDSCs described as "macrophages, neutrophils". | `compare_adapter.py --pubmed` on 31 held-out prompts, criteria fixed before training. Base / adapter v2: citations 0 / 0; unsupported numbers 0 / 0; unsupported cell types 0 / 6 (all MDSC); 170-250 words 26 / 27; all four sections 29 / 31; MoA terms 3.48 / 3.35. Eval loss 1.32, 325 s. Verdict: fails only on cell types. |
| L11 | P1 | **Real breast-atlas labels broke the cell-type guard (fixed 2026-10-21).** The annotated breast sample uses CellTypist atlas labels (`LummHR-SCGB`, `Fibro-major`, `Macro-m2`, `plasma_IgG`, `Lymph-major`, ...); 34 of its 38 labels named no lineage the guard recognises, so on real data correct sentences about epithelium, fibroblasts or macrophages would have been stripped. Fix: `faithfulness.label_lineage` maps atlas labels to a plain lineage and the prompt lists them as `LummHR-SCGB (epithelial cells)`. Plain labels (`T_cells`, `Stromal_CAF`) get no suffix, so the training prompts and adapter v2 are unchanged. Remaining: atlas labels never appear in the synthetic training data. | `src/report_generation/faithfulness.py`, `tests/test_faithfulness.py` |
| L7 | P2 | `max_steps=150 × grad_accum 8` on 206 examples ≈ 5.8 epochs of a tiny synthetic set — overfitting risk. Gemma 3 in fp16 on T4 (no bf16) is prone to NaN loss. | — |

## 2. Spatial pipeline → prompt data contract

All of the following were reproduced by running `annotate_spatial_regions` +
`calculate_spatial_heterogeneity` on a synthetic 400-spot tissue whose left
half is T-cell-rich and right half epithelial (see `tests/conftest.py`).

| # | Sev | Finding | Evidence |
|---|-----|---------|----------|
| S1 | P0 | **Key mismatch zeroes every spatial metric in the prompt.** Adapter emits nested `morans_i.mean` / `neighborhood_enrichment.n_enriched_pairs`; prompt builder reads flat `morans_i_mean` / `n_enriched_pairs`. | Pipeline computed Moran's I = 0.087; prompt said `Moran's I: weak (0.00)` and `co-enriched pairs: 0`. |
| S2 | P0 | **Immune cells not recognised.** Prompt builders match `'t cells'`, but labels are `T_cells`. | 48% T-cell tissue rendered as `Tissue phenotype: cold immune desert`. |
| S3 | P0 | **Fabricated metrics.** `mean_confidence` reads a `leiden_confidence` column that is never created → silently defaults to `0.8`. `mean_prediction_entropy` reads a never-created column → hard-coded `1.0`. | Both values appear identically for every input. |
| S4 | P1 | `compute_multiscale_neighborhood_enrichment` never returns `n_significant_pairs`/`max_enrichment`, which the adapter reads → always `0`/`1.0`. It also overwrites `obsp['spatial_connectivities']` with the radius-3 graph, which spatial entropy then uses. | — |
| S5 | P1 | Moran's I builds a dense `n×n` `np.outer` per gene per permutation (50 genes × 100 perms). Fine at 400 spots; infeasible for Visium HD bins. | `uncertainty_spatial_analysis.py` |
| S6 | P1 | "Already preprocessed" heuristic `n_vars < 3000` sends raw counts to PCA/CellTypist for targeted panels. | CellTypist: `Invalid expression matrix … expect log1p normalized`. |
| S6b | P1 | A typical processed upload (X = 2,000 scaled HVGs, full log-normalised matrix in `.raw`) was scored on the scaled HVGs only. Now restored from `.raw`. | Real breast Visium (`outputs/processed_visium.h5ad`, 4,895 spots): 2,000 → 21,351 genes used. |
| S8 | P1 | **Open, not fixed: pan-tissue marker panels mislabel spots.** Per-gene z-score + argmax always picks a label, and sparse genes get large z-scores. On the breast sample, 230 spots are called `Hepatocyte` (ALB/HP/APOA1 are genuinely spatially enriched there) and 238 `Neuronal`; a 1.0 z threshold still leaves 151. This is a method limitation, so the fix is a design choice: tissue-aware panels, a reference-based deconvolution (cell2location / RCTD), or CellTypist on all spots with a tissue model. | Same sample. |
| S7 | P1 | Every stage is wrapped in `except Exception: print(...); return zeros`. Failures are indistinguishable from biology. | — |

## 3. Inference / deployment

| # | Sev | Finding | Evidence |
|---|-----|---------|----------|
| I1 | P0 | **MedGemma cannot load with the pinned dependencies.** `transformers==4.45.1` predates Gemma 3 (added in 4.50); the `image-text-to-text` pipeline does not exist in 4.45.1. | `AutoConfig.from_pretrained("google/medgemma-4b-it")` → `model type gemma3 … not recognize this architecture`. |
| I2 | P1 | 4-bit bitsandbytes is requested unconditionally, but the HF Space is CPU-only and M1 has no CUDA. | `generate_textonly_fallback` |
| I3 | P1 | No chat template on an instruction-tuned model; output sliced by string length instead of generated token ids. | — |
| I4 | P1 | Model reloaded on every Streamlit rerun (no `st.cache_resource`). | — |
| I5 | P1 | Deploy workflow pushed to the Space on every push to `main`, independent of CI, using an inline Dockerfile that had drifted from `app/Dockerfile`. | — |

## 4. CI/CD and compliance

| # | Sev | Finding |
|---|-----|---------|
| C1 | P0 | CI could not fail: every step ended in `|| true` / `|| echo`; `tests/` did not exist. |
| C2 | P1 | Three conflicting dependency specs (`requirements.txt`, `app/requirements.txt`, `requirements_full.txt`). |
| C2b | P1 | The scverse stack on Python 3.10 is stuck on squidpy 1.6 + legacy dask; with transformers < 5 the lock carried 8 known CVEs (pip-audit). Moving the lock to Python 3.12 (Kaggle's current runtime) gives squidpy 1.8, transformers 5.x and a clean audit. |
| C2c | P1 | On macOS/arm64 + Python 3.12, `import torch` before `sc.pp.neighbors` segfaults (numba OpenMP vs torch libomp). Fixed by `numba.config.THREADING_LAYER = "workqueue"` in `src/__init__.py`. |
| C3 | P1 | `.gitleaks.toml` allow-listed `.env` (the file most likely to hold a real secret) and its Kaggle rule used an unsupported `paths` key, producing 140+ false positives on MD5 sample IDs. |
| C4 | P1 | `scripts/push_lora_kernel.sh` substitutes the real `HF_TOKEN` into notebook source before pushing a **public** Kaggle kernel (`is_private: false`). The current notebook uses Kaggle Secrets so nothing is substituted today, but any earlier version pushed with the placeholder would have published the token. **Recommendation: rotate the HF token.** |
| C5 | P1 | **Training targets contain wrong citations.** Spot-checked against PubMed: Helmink et al. 2020 is *Nature* ([10.1038/s41586-019-1922-8](https://doi.org/10.1038/s41586-019-1922-8)) but is cited as *Nature Medicine* 34×; Kieffer et al. 2020 is *Cancer Discovery* ([10.1158/2159-8290.CD-19-1384](https://doi.org/10.1158/2159-8290.CD-19-1384)) but is cited as *Nature* 32× and *Cancer Cell* 7×; Salmon et al. 2012 is *J Clin Invest* ([10.1172/JCI45817](https://doi.org/10.1172/JCI45817)) but is cited as *Immunity* 18×. Fine-tuning on these teaches the model to emit confident, wrong references in a clinical-facing report. |
| C6 | P2 | Training data was generated by an external LLM API. Check that provider's terms on using outputs to train other models before publishing the adapter. |
| C7 | P2 | The app recommends therapies; it needs a visible "research use only, not clinical advice" disclaimer. |

## Remediation plan (one PR each, into `devel`)

1. **CI that can fail** — ruff (correctness rules), pytest on synthetic AnnData, notebook validation, gitleaks over full history, pip-audit, Docker build + container health check; single locked `requirements.txt`; deploy gated on CI. *(this PR)*
2. **Pipeline data contract** — S1–S7, with regression tests that pin the evidence above.
3. **Inference** — I1–I4: Gemma 3–capable transformers, chat template, CPU-safe loading, model caching, disclaimer.
4. **LoRA** — training pipeline fixed and run (see L9: data must be regenerated before the adapter is usable); L1–L7: max length ≥ 1024, completion-only loss, train/serve prompt parity, a notebook that runs under "Save & Run All", and a test that every example keeps its completion tokens.

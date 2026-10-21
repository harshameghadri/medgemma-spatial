#!/usr/bin/env python3
"""LoRA / QLoRA fine-tuning of MedGemma (or a text-only causal LM such as Llama 3.1) on spatial pathology reports.

Runs the same way on a Kaggle T4, a CUDA workstation, an M1 Mac or CPU (the
last two only for smoke tests). Expects prompt-completion JSONL produced by
``scripts/build_lora_dataset.py``.

Guards against the failure modes found in docs/REVIEW.md:
  * every example must fit in --max-length (no silent truncation of targets)
  * loss is computed on completion tokens only
  * training stops on a non-finite loss instead of saving a broken adapter

Examples:
    # 2-step CPU smoke test with a tiny checkpoint
    python scripts/train_lora.py --model-id trl-internal-testing/tiny-Gemma3ForConditionalGeneration \
        --data-dir tests/fixtures/lora --max-steps 2 --quantize none

    # Kaggle T4
    python scripts/train_lora.py --data-dir /kaggle/input/medgemma-spatial-data \
        --output-dir /kaggle/working/adapter --push-to-hub user/medgemma-spatial-adapter
"""

import argparse
import math
import os
import sys
from pathlib import Path

DEFAULT_TARGETS = r".*language_model.*\.(q_proj|k_proj|v_proj|o_proj)"
CAUSAL_LM_TARGETS = r".*\.(q_proj|k_proj|v_proj|o_proj)"


def default_targets(config) -> str:
    """Attention projections of the language model: Gemma 3 VLMs nest it, text-only LMs (Llama) do not."""
    return DEFAULT_TARGETS if hasattr(config, "vision_config") else CAUSAL_LM_TARGETS


def parse_args(argv=None):
    """Parse CLI arguments."""
    p = argparse.ArgumentParser(description="MedGemma LoRA fine-tuning")
    p.add_argument("--model-id", default="google/medgemma-4b-it")
    p.add_argument("--data-dir", default="data/lora_training_v2",
                   help="Directory (searched recursively) containing train.jsonl and eval.jsonl")
    p.add_argument("--output-dir", default="outputs/lora_adapter")
    p.add_argument("--max-length", type=int, default=1024,
                   help="Longest rendered example is 944 tokens with the serving prompt")
    p.add_argument("--epochs", type=float, default=2.0)
    p.add_argument("--max-steps", type=int, default=-1, help="-1 = use --epochs")
    p.add_argument("--learning-rate", type=float, default=1e-4)
    p.add_argument("--batch-size", type=int, default=1)
    p.add_argument("--grad-accum", type=int, default=8)
    p.add_argument("--lora-r", type=int, default=16)
    p.add_argument("--lora-alpha", type=int, default=32)
    p.add_argument("--lora-dropout", type=float, default=0.05)
    p.add_argument("--target-modules", default=None,
                   help="Regex over module names; default = language-model attention only")
    p.add_argument("--quantize", choices=["auto", "4bit", "none"], default="auto",
                   help="auto = 4-bit when CUDA + bitsandbytes are available")
    p.add_argument("--precision", choices=["auto", "bf16", "fp16", "fp32"], default="auto",
                   help="auto = bf16 if supported, else fp32 compute (Gemma 3 overflows in fp16)")
    p.add_argument("--push-to-hub", default=None, help="HF repo id for the adapter")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def find_data_files(data_dir: str):
    """Locate train.jsonl / eval.jsonl, searching subdirectories (Kaggle inputs nest)."""
    root = Path(data_dir)
    direct = root / "train.jsonl"
    hits = [direct] if direct.exists() else sorted(root.rglob("train.jsonl"))
    if not hits:
        raise FileNotFoundError(f"train.jsonl not found under {root}")
    train = hits[0]
    evals = train.parent / "eval.jsonl"
    return train, (evals if evals.exists() else None)


def load_data(data_dir: str):
    """Load prompt-completion JSONL into a DatasetDict."""
    from datasets import load_dataset

    train, evals = find_data_files(data_dir)
    files = {"train": str(train)}
    if evals:
        files["validation"] = str(evals)
    ds = load_dataset("json", data_files=files)
    missing = {"prompt", "completion"} - set(ds["train"].column_names)
    if missing:
        raise ValueError(
            f"{train} lacks {sorted(missing)}; build it with scripts/build_lora_dataset.py "
            "(the old single-'text' format trains on the prompt)")
    print(f"Data: {train.parent} — train {len(ds['train'])}"
          + (f", eval {len(ds['validation'])}" if evals else ""))
    return ds


def token_lengths(tokenizer, example):
    """(prompt_tokens, total_tokens) for one prompt-completion example via the chat template."""
    prompt_ids = tokenizer.apply_chat_template(
        example["prompt"], tokenize=True, add_generation_prompt=True)
    full_ids = tokenizer.apply_chat_template(
        example["prompt"] + example["completion"], tokenize=True)
    if hasattr(prompt_ids, "keys"):
        prompt_ids, full_ids = prompt_ids["input_ids"], full_ids["input_ids"]
    return len(prompt_ids), len(full_ids)


def check_lengths(tokenizer, dataset, max_length: int, min_completion_tokens: int = 50):
    """Fail fast if any example would be truncated or has almost no target tokens."""
    worst_total, worst_completion = 0, math.inf
    for i, ex in enumerate(dataset):
        n_prompt, n_total = token_lengths(tokenizer, ex)
        n_completion = n_total - n_prompt
        worst_total = max(worst_total, n_total)
        worst_completion = min(worst_completion, n_completion)
        if n_total > max_length:
            raise ValueError(
                f"Example {i} is {n_total} tokens > --max-length {max_length}; it would be "
                f"truncated and lose {n_total - max_length} target tokens")
        if n_completion < min_completion_tokens:
            raise ValueError(f"Example {i} has only {n_completion} completion tokens")
    print(f"Length check: longest example {worst_total} tokens (limit {max_length}); "
          f"shortest completion {worst_completion} tokens")
    return worst_total


def resolve_precision(choice: str):
    """Return (torch_dtype, bf16_flag, fp16_flag)."""
    import torch

    if choice == "auto":
        if torch.cuda.is_available() and torch.cuda.is_bf16_supported():
            choice = "bf16"
        else:
            choice = "fp32"
    if choice == "bf16":
        return torch.bfloat16, True, False
    if choice == "fp16":
        return torch.float16, False, True
    return torch.float32, False, False


def load_model(model_id: str, quantize: str, dtype):
    """Load MedGemma (Gemma 3 image-text-to-text) or a text-only causal LM, optionally in 4-bit NF4."""
    import torch
    from transformers import AutoConfig, AutoModelForCausalLM, AutoModelForImageTextToText

    token = os.environ.get("HF_TOKEN")
    config = AutoConfig.from_pretrained(model_id, token=token)
    vlm = hasattr(config, "vision_config")
    prequantized = bool(getattr(config, "quantization_config", None))
    use_4bit = prequantized or quantize == "4bit" or (quantize == "auto" and torch.cuda.is_available())
    kwargs = {"dtype": dtype, "token": token}
    if vlm:
        kwargs["attn_implementation"] = "eager"
    if prequantized:
        kwargs["device_map"] = {"": 0}
    elif use_4bit:
        from transformers import BitsAndBytesConfig

        kwargs["quantization_config"] = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=dtype,
        )
        kwargs["device_map"] = {"": 0}
    model = (AutoModelForImageTextToText if vlm else AutoModelForCausalLM).from_pretrained(model_id, **kwargs)
    model.config.use_cache = False
    print(f"Model: {model_id} ({type(model).__name__}, 4-bit={use_4bit}, dtype={dtype})")
    return model, use_4bit


def build_trainer(args, model, tokenizer, dataset, use_4bit, bf16, fp16):
    """Configure LoRA + TRL SFTTrainer with completion-only loss."""
    from peft import LoraConfig, prepare_model_for_kbit_training
    from transformers import TrainerCallback
    from trl import SFTConfig, SFTTrainer

    class StopOnNonFiniteLoss(TrainerCallback):
        def on_log(self, _args, state, control, logs=None, **kwargs):
            loss = (logs or {}).get("loss")
            if loss is not None and not math.isfinite(loss):
                raise RuntimeError(f"Non-finite training loss at step {state.global_step}: {loss}")

    if use_4bit:
        model = prepare_model_for_kbit_training(model, use_gradient_checkpointing=True)

    lora = LoraConfig(
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        target_modules=args.target_modules or default_targets(model.config),
        bias="none",
        task_type="CAUSAL_LM",
    )
    has_eval = "validation" in dataset
    config = SFTConfig(
        output_dir=args.output_dir,
        num_train_epochs=args.epochs,
        max_steps=args.max_steps,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        gradient_accumulation_steps=args.grad_accum,
        gradient_checkpointing=True,
        gradient_checkpointing_kwargs={"use_reentrant": False},
        learning_rate=args.learning_rate,
        lr_scheduler_type="cosine",
        warmup_steps=0.05,
        weight_decay=0.0,
        max_grad_norm=1.0,
        max_length=args.max_length,
        completion_only_loss=True,
        packing=False,
        bf16=bf16,
        fp16=fp16,
        optim="paged_adamw_8bit" if use_4bit else "adamw_torch",
        logging_steps=5 if args.max_steps < 0 else 1,
        eval_strategy="epoch" if has_eval and args.max_steps < 0 else "no",
        save_strategy="epoch" if args.max_steps < 0 else "no",
        save_total_limit=1,
        report_to="none",
        seed=args.seed,
        dataloader_pin_memory=False,
    )
    return SFTTrainer(
        model=model,
        args=config,
        train_dataset=dataset["train"],
        eval_dataset=dataset["validation"] if has_eval else None,
        processing_class=tokenizer,
        peft_config=lora,
        callbacks=[StopOnNonFiniteLoss()],
    )


def main(argv=None):
    """Validate data, train, save (and optionally push) the adapter."""
    args = parse_args(argv)
    from transformers import AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(args.model_id, token=os.environ.get("HF_TOKEN"))
    dataset = load_data(args.data_dir)
    for split in dataset:
        check_lengths(tokenizer, dataset[split], args.max_length)

    dtype, bf16, fp16 = resolve_precision(args.precision)
    model, use_4bit = load_model(args.model_id, args.quantize, dtype)
    trainer = build_trainer(args, model, tokenizer, dataset, use_4bit, bf16, fp16)
    trainer.model.print_trainable_parameters()

    result = trainer.train()
    print(f"Training complete: {result.metrics}")
    if trainer.eval_dataset is not None:
        print(f"Eval: {trainer.evaluate()}")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    trainer.model.save_pretrained(str(out))
    tokenizer.save_pretrained(str(out))
    print(f"Adapter saved to {out.resolve()}")

    if args.push_to_hub:
        trainer.model.push_to_hub(args.push_to_hub, token=os.environ.get("HF_TOKEN"), private=True)
        tokenizer.push_to_hub(args.push_to_hub, token=os.environ.get("HF_TOKEN"), private=True)
        print(f"Pushed (private) to https://huggingface.co/{args.push_to_hub}")
    return trainer


if __name__ == "__main__":
    sys.exit(0 if main() else 1)

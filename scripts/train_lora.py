#!/usr/bin/env python3
"""Local smoke-test script for MedGemma LoRA fine-tuning (CPU / M1 Mac compatible)."""

import argparse
import os
from pathlib import Path


def parse_args():
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="MedGemma LoRA fine-tuning (local test)")
    parser.add_argument("--max-steps", type=int, default=-1,
                        help="Maximum training steps (-1 = full training, 1 = smoke test)")
    parser.add_argument("--no-gpu", action="store_true",
                        help="Skip BitsAndBytes quantization; use float32 on CPU")
    parser.add_argument("--data-dir", type=str, default="data/lora_training/",
                        help="Directory containing train.jsonl and eval.jsonl")
    parser.add_argument("--output-dir", type=str, default="outputs/lora_adapter/",
                        help="Output directory for the trained adapter")
    return parser.parse_args()


def load_data(data_dir: str):
    """Load train/eval JSONL datasets."""
    from datasets import load_dataset

    data_path = Path(data_dir)
    train_file = data_path / "train.jsonl"
    eval_file = data_path / "eval.jsonl"

    if not train_file.exists():
        raise FileNotFoundError(f"Training data not found at {train_file}")
    if not eval_file.exists():
        raise FileNotFoundError(f"Eval data not found at {eval_file}")

    dataset = load_dataset("json", data_files={
        "train": str(train_file),
        "validation": str(eval_file),
    })
    print(f"Data loaded — Train: {len(dataset['train'])} | Eval: {len(dataset['validation'])}")
    return dataset


def load_model_and_tokenizer(no_gpu: bool, hf_token):
    """Load model with optional 4-bit quantization."""
    import torch
    from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig

    model_id = "google/medgemma-4b-it"

    tokenizer = AutoTokenizer.from_pretrained(model_id, token=hf_token)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    if no_gpu:
        print("CPU mode: loading in float32, no quantization")
        model = AutoModelForCausalLM.from_pretrained(
            model_id,
            torch_dtype=torch.float32,
            device_map="cpu",
            token=hf_token,
        )
    else:
        print("GPU mode: loading with 4-bit NF4 quantization")
        bnb_config = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_compute_dtype=torch.float16,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
        )
        model = AutoModelForCausalLM.from_pretrained(
            model_id,
            quantization_config=bnb_config,
            device_map="auto",
            token=hf_token,
        )

    model.config.use_cache = False
    device = next(model.parameters()).device
    print(f"Model loaded on {device}")
    return model, tokenizer


def attach_lora(model, no_gpu: bool):
    """Attach LoRA adapter to model."""
    from peft import LoraConfig, get_peft_model, prepare_model_for_kbit_training

    if not no_gpu:
        model = prepare_model_for_kbit_training(model)

    lora_config = LoraConfig(
        r=16,
        lora_alpha=32,
        lora_dropout=0.05,
        target_modules=["q_proj", "k_proj", "v_proj", "o_proj",
                        "gate_proj", "up_proj", "down_proj"],
        bias="none",
        task_type="CAUSAL_LM",
    )
    model = get_peft_model(model, lora_config)
    model.print_trainable_parameters()
    print("LoRA adapter attached")
    return model


def build_training_args(output_dir: str, max_steps: int, no_gpu: bool):
    """Build SFTConfig training arguments."""
    from trl import SFTConfig

    kwargs = dict(
        output_dir=output_dir,
        num_train_epochs=3 if max_steps == -1 else 1,
        per_device_train_batch_size=1,
        gradient_accumulation_steps=1 if no_gpu else 8,
        learning_rate=2e-4,
        max_seq_length=128 if no_gpu else 768,
        lr_scheduler_type="cosine",
        warmup_ratio=0.05,
        logging_steps=1,
        save_strategy="no" if max_steps == 1 else "epoch",
        eval_strategy="no" if max_steps == 1 else "epoch",
        fp16=False if no_gpu else True,
        optim="adamw_torch" if no_gpu else "paged_adamw_8bit",
        report_to="none",
        dataset_text_field="text",
    )
    if max_steps > 0:
        kwargs["max_steps"] = max_steps

    return SFTConfig(**kwargs)


def main():
    """Run LoRA fine-tuning smoke test or full training."""
    args = parse_args()

    hf_token = os.environ.get("HF_TOKEN")
    if not hf_token:
        print("WARNING: HF_TOKEN not set. Will attempt unauthenticated access.")

    print("=== MedGemma LoRA Training ===")
    print(f"  data-dir   : {args.data_dir}")
    print(f"  output-dir : {args.output_dir}")
    print(f"  max-steps  : {args.max_steps} ({'smoke test' if args.max_steps == 1 else 'full'})")
    print(f"  no-gpu     : {args.no_gpu}")

    dataset = load_data(args.data_dir)
    model, tokenizer = load_model_and_tokenizer(no_gpu=args.no_gpu, hf_token=hf_token)
    model = attach_lora(model, no_gpu=args.no_gpu)

    from trl import SFTTrainer

    training_args = build_training_args(
        output_dir=args.output_dir,
        max_steps=args.max_steps,
        no_gpu=args.no_gpu,
    )

    trainer = SFTTrainer(
        model=model,
        train_dataset=dataset["train"],
        eval_dataset=dataset["validation"] if args.max_steps != 1 else None,
        args=training_args,
        tokenizer=tokenizer,
    )

    print("Training started...")
    trainer.train()
    print("Training complete")

    if args.max_steps != 1:
        metrics = trainer.evaluate()
        print(f"Final eval loss: {metrics.get('eval_loss', 'N/A'):.4f}")

    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    trainer.model.save_pretrained(str(output_path))
    tokenizer.save_pretrained(str(output_path))
    print(f"Adapter saved to {output_path.resolve()}")


if __name__ == "__main__":
    main()

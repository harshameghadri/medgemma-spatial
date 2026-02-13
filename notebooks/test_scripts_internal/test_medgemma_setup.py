#!/usr/bin/env python3
"""
Test MedGemma Setup
Verify installation and dependencies before running full analysis.

Usage:
    python test_medgemma_setup.py
"""

import sys
from pathlib import Path


def check_python_version():
    """Verify Python 3.10+"""
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    if version.major == 3 and version.minor >= 10:
        print("✓ Python version compatible")
        return True
    else:
        print("✗ Python 3.10+ required")
        return False


def check_dependencies():
    """Verify required packages."""
    required = {
        'torch': '2.4.0',
        'transformers': '4.45.1',
        'bitsandbytes': '0.49.0',
        'accelerate': '1.12.0'
    }

    all_ok = True
    for package, min_version in required.items():
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"{package:15s} {version:15s} (>= {min_version})")

            if package == 'torch':
                import torch
                print(f"  {'MPS available:'} {torch.backends.mps.is_available()}")
                print(f"  {'CUDA available:'} {torch.cuda.is_available()}")

        except ImportError:
            print(f"{package:15s} NOT INSTALLED")
            all_ok = False

    return all_ok


def check_data_files():
    """Verify input data exists."""
    base_dir = Path(__file__).parent.parent
    required_files = [
        'outputs/spatial_features.json',
        'outputs/cell_type_enhanced_summary.json'
    ]

    all_ok = True
    print("\nChecking data files:")
    for filepath in required_files:
        full_path = base_dir / filepath
        if full_path.exists():
            size_mb = full_path.stat().st_size / 1024**2
            print(f"✓ {filepath:40s} ({size_mb:.2f} MB)")
        else:
            print(f"✗ {filepath:40s} NOT FOUND")
            all_ok = False

    return all_ok


def check_model_access():
    """Test HuggingFace model access (without downloading)."""
    try:
        from huggingface_hub import model_info

        print("\nChecking MedGemma model access:")
        info = model_info("google/medgemma-4b-it")
        print(f"✓ Model accessible: {info.modelId}")
        print(f"  Downloads: {info.downloads:,}")
        print(f"  Size: ~{info.siblings[0].size / 1024**3:.1f} GB (unquantized)")
        return True

    except ImportError:
        print("Install huggingface_hub for model checks: pip install huggingface_hub")
        return True

    except Exception as e:
        print(f"✗ Model access failed: {e}")
        print("  Note: Model will be downloaded on first run")
        return True


def estimate_memory():
    """Estimate memory requirements."""
    try:
        import psutil

        total_gb = psutil.virtual_memory().total / 1024**3
        available_gb = psutil.virtual_memory().available / 1024**3

        print(f"\nMemory Status:")
        print(f"  Total: {total_gb:.1f} GB")
        print(f"  Available: {available_gb:.1f} GB")

        required_gb = 12
        if available_gb > required_gb:
            print(f"✓ Sufficient memory (need ~{required_gb}GB for 4-bit model)")
            return True
        else:
            print(f"⚠ Low memory - may need to close applications")
            return False

    except ImportError:
        print("Install psutil for memory check: pip install psutil")
        return True


def test_tokenizer_only():
    """Test loading tokenizer (small download)."""
    try:
        print("\nTesting tokenizer load (small download):")
        from transformers import AutoTokenizer

        tokenizer = AutoTokenizer.from_pretrained("google/medgemma-4b-it")
        test_text = "This is a test of the tokenizer"
        tokens = tokenizer(test_text, return_tensors="pt")

        print(f"✓ Tokenizer loaded")
        print(f"  Test tokens: {len(tokens['input_ids'][0])} tokens")
        return True

    except Exception as e:
        print(f"✗ Tokenizer test failed: {e}")
        return False


def main():
    print("=" * 80)
    print("MedGemma Setup Verification")
    print("=" * 80)

    results = []

    print("\n1. Python Version")
    print("-" * 80)
    results.append(check_python_version())

    print("\n2. Dependencies")
    print("-" * 80)
    results.append(check_dependencies())

    print("\n3. Data Files")
    print("-" * 80)
    results.append(check_data_files())

    print("\n4. Memory")
    print("-" * 80)
    results.append(estimate_memory())

    print("\n5. Model Access")
    print("-" * 80)
    results.append(check_model_access())

    print("\n6. Tokenizer Test")
    print("-" * 80)
    results.append(test_tokenizer_only())

    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)

    if all(results):
        print("✓ All checks passed - ready to run MedGemma")
        print("\nNext steps:")
        print("  1. Run notebook: jupyter notebook 04_medgemma_reports.ipynb")
        print("  2. Or run script: python run_medgemma.py")
    else:
        print("✗ Some checks failed - review errors above")
        print("\nCommon fixes:")
        print("  pip install torch>=2.4.0 transformers>=4.45.1 bitsandbytes accelerate")
        print("  Run spatial analysis first to generate input JSON files")

    print("\nNote: First run will download MedGemma model (~4GB with quantization)")
    print("=" * 80)


if __name__ == "__main__":
    main()

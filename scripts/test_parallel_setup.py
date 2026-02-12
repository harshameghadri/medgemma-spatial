#!/usr/bin/env python3
"""
Quick validation test for parallel pipeline setup.
Verifies all components without running full analysis.
"""

import sys
import json
from pathlib import Path

def test_obfuscation():
    """Test obfuscation system."""
    print("Testing obfuscation system...")

    mapping_path = Path("data/file_mapping.json")

    if not mapping_path.exists():
        print("  ✗ file_mapping.json not found")
        return False

    with open(mapping_path, 'r') as f:
        mapping = json.load(f)

    total_files = mapping['summary']['total_files']
    expected_files = 12  # 3 files × 4 tissues

    if total_files != expected_files:
        print(f"  ✗ Expected {expected_files} files, found {total_files}")
        return False

    print(f"  ✓ Found {total_files} obfuscated files")

    # Check symlinks exist
    obf_dir = Path("data/obfuscated")
    if not obf_dir.exists():
        print("  ✗ data/obfuscated directory not found")
        return False

    symlink_count = sum(1 for f in obf_dir.iterdir() if f.is_symlink())
    print(f"  ✓ {symlink_count} symlinks verified")

    return True


def test_scripts():
    """Test required scripts exist."""
    print("\nTesting pipeline scripts...")

    required_scripts = [
        "scripts/obfuscate_data.py",
        "scripts/parallel_pipeline.py",
        "notebooks/uncertainty_spatial_analysis.py",
        "notebooks/medgemma_v2_pipeline.py"
    ]

    missing = []
    for script in required_scripts:
        if not Path(script).exists():
            missing.append(script)
            print(f"  ✗ Missing: {script}")
        else:
            print(f"  ✓ Found: {script}")

    return len(missing) == 0


def test_data_integrity():
    """Verify raw data files."""
    print("\nTesting raw data integrity...")

    mapping_path = Path("data/file_mapping.json")
    with open(mapping_path, 'r') as f:
        mapping = json.load(f)

    total_size_gb = mapping['summary']['total_size_gb']
    expected_min_gb = 35.0  # Should be ~39 GB

    if total_size_gb < expected_min_gb:
        print(f"  ✗ Total size {total_size_gb} GB < expected {expected_min_gb} GB")
        return False

    print(f"  ✓ Total data size: {total_size_gb} GB")

    # Check each tissue has binned_outputs
    tissues = set()
    for info in mapping['files'].values():
        if 'binned_outputs' in info['original_name']:
            tissues.add(info['tissue'])

    expected_tissues = {'colon', 'lung', 'pancreas', 'prostate'}
    if tissues != expected_tissues:
        print(f"  ✗ Found tissues {tissues}, expected {expected_tissues}")
        return False

    print(f"  ✓ All tissues present: {', '.join(sorted(tissues))}")

    return True


def test_imports():
    """Test Python dependencies."""
    print("\nTesting Python dependencies...")

    required_packages = {
        'scanpy': 'scanpy',
        'squidpy': 'squidpy',
        'torch': 'torch',
        'transformers': 'transformers',
        'anndata': 'anndata'
    }

    missing = []
    for package, import_name in required_packages.items():
        try:
            __import__(import_name)
            print(f"  ✓ {package}")
        except ImportError:
            print(f"  ✗ {package} not installed")
            missing.append(package)

    if missing:
        print(f"\nInstall missing packages:")
        print(f"  pip install {' '.join(missing)}")
        return False

    return True


def main():
    """Run all validation tests."""
    print("=" * 80)
    print("PARALLEL PIPELINE VALIDATION TEST")
    print("=" * 80)

    tests = [
        ("Obfuscation System", test_obfuscation),
        ("Pipeline Scripts", test_scripts),
        ("Data Integrity", test_data_integrity),
        ("Python Dependencies", test_imports)
    ]

    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\n  ERROR: {e}")
            results.append((test_name, False))

    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)

    all_passed = all(passed for _, passed in results)

    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status:8} {test_name}")

    print("=" * 80)

    if all_passed:
        print("\n✓ ALL TESTS PASSED - Ready to run parallel pipeline!")
        print("\nNext steps:")
        print("  1. Test single sample:")
        print("     python scripts/parallel_pipeline.py --tissues colon --skip-medgemma")
        print("\n  2. Run full pipeline:")
        print("     python scripts/parallel_pipeline.py --workers 4")
    else:
        print("\n✗ SOME TESTS FAILED - Fix issues before running pipeline")
        sys.exit(1)


if __name__ == "__main__":
    main()

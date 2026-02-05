#!/usr/bin/env python3
"""
Parallel processing pipeline for MedGemma spatial analysis.
Handles multiple obfuscated samples with optimal resource utilization.
"""

import os
import sys
import json
import torch
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple
import multiprocessing as mp
from datetime import datetime

sys.path.append(str(Path(__file__).parent.parent / 'notebooks'))


def get_available_cores():
    """
    Determine optimal number of parallel workers.

    Rules:
    - M1 Mac: Use 4-6 cores max (MedGemma is memory-intensive)
    - CPU-only: Use n_cores - 2
    - GPU: Single process (GPU not parallelizable for single model)
    """
    total_cores = mp.cpu_count()

    if torch.cuda.is_available():
        print("GPU detected - using single process (GPU exclusive)")
        return 1

    # M1 Mac or CPU
    optimal_cores = min(6, max(1, total_cores - 2))
    print(f"CPU mode: Using {optimal_cores} of {total_cores} cores")
    return optimal_cores


def load_file_mapping(mapping_path: str) -> Dict:
    """Load obfuscated file mapping."""
    with open(mapping_path, 'r') as f:
        return json.load(f)


def find_binned_outputs(obfuscated_dir: Path, mapping: Dict) -> List[Tuple[str, str, str]]:
    """
    Find all binned_outputs.tar.gz files.

    Returns:
        List of (hex_name, tissue_type, original_name) tuples
    """
    samples = []

    for hex_name, info in mapping['files'].items():
        if 'binned_outputs.tar.gz' in info['original_name']:
            tissue = info['tissue']
            original = info['original_name']
            samples.append((hex_name, tissue, original))

    return sorted(samples)


def detect_available_bins(binned_tar_path: str) -> list:
    """Detect all available bin sizes in archive."""
    import tarfile
    import re

    with tarfile.open(binned_tar_path, 'r:gz') as tar:
        members = tar.getnames()

        # Find all square_XXXum patterns
        bin_pattern = re.compile(r'square_(\d{3}um)')
        bins = set()
        for m in members:
            match = bin_pattern.search(m)
            if match:
                bins.add(match.group(1))

        # Sort by numeric value (064um > 032um > 016um > 008um)
        return sorted(bins, key=lambda x: int(x[:3]), reverse=True)


def extract_binned_h5ad(binned_tar_path: str, output_dir: Path, bin_size: str = None) -> str:
    """
    Extract and convert h5 to h5ad from binned_outputs.tar.gz.

    Args:
        binned_tar_path: Path to binned_outputs.tar.gz
        output_dir: Where to extract
        bin_size: Resolution (008um, 016um, etc.) - auto-detects if None

    Returns:
        Path to h5ad file
    """
    import tarfile
    import scanpy as sc

    output_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect if not specified
    if bin_size is None:
        available = detect_available_bins(binned_tar_path)
        if not available:
            raise FileNotFoundError("No square_XXXum bins found in archive")
        bin_size = available[0]  # Highest bin = lowest memory
        print(f"  Auto-detected bin sizes: {available}")
        print(f"  Using: {bin_size} (highest = least spots = lowest memory)")

    with tarfile.open(binned_tar_path, 'r:gz') as tar:
        # Try .h5ad first (newer format)
        h5ad_members = [
            m for m in tar.getmembers()
            if f"square_{bin_size}" in m.name and m.name.endswith('.h5ad')
        ]

        if h5ad_members:
            member = h5ad_members[0]
            tar.extract(member, path=output_dir)
            return str(output_dir / member.name)

        # Fall back to .h5 (older 10x format)
        h5_members = [
            m for m in tar.getmembers()
            if f"square_{bin_size}" in m.name
            and 'filtered_feature_bc_matrix.h5' in m.name
        ]

        if not h5_members:
            raise FileNotFoundError(
                f"No .h5ad or .h5 found for bin {bin_size}. "
                f"Available bins: {detect_available_bins(binned_tar_path)}"
            )

        # Extract and convert .h5 → .h5ad
        member = h5_members[0]
        tar.extract(member, path=output_dir)
        h5_path = output_dir / member.name

        print(f"  Converting {h5_path.name} to h5ad...")
        adata = sc.read_10x_h5(h5_path)

        # Extract spatial coordinates
        spatial_members = [
            m for m in tar.getmembers()
            if f"square_{bin_size}" in m.name
            and 'tissue_positions.parquet' in m.name
        ]

        if spatial_members:
            import pandas as pd
            spatial_member = spatial_members[0]
            tar.extract(spatial_member, path=output_dir)
            spatial_path = output_dir / spatial_member.name

            print(f"  Loading spatial coordinates from {spatial_path.name}...")
            spatial_df = pd.read_parquet(spatial_path)

            # Add spatial coordinates to adata.obsm['spatial']
            # Match barcodes in adata to spatial_df
            if 'barcode' in spatial_df.columns:
                spatial_df = spatial_df.set_index('barcode')

            # Get x, y coordinates (column names vary by 10x version)
            coord_cols = [c for c in spatial_df.columns if 'pxl' in c.lower() or 'imagerow' in c.lower() or 'imagecol' in c.lower()]

            if len(coord_cols) >= 2:
                # Use first two coordinate columns
                coords = spatial_df.loc[adata.obs_names, coord_cols[:2]].values
                adata.obsm['spatial'] = coords
                print(f"  Added spatial coordinates: {coords.shape}")
            else:
                print(f"  Warning: Could not find coordinate columns in {list(spatial_df.columns)}")
        else:
            print(f"  Warning: No tissue_positions.parquet found for {bin_size}")

        h5ad_path = output_dir / f"square_{bin_size}.h5ad"
        adata.write_h5ad(h5ad_path)
        print(f"  Saved: {h5ad_path}")

        return str(h5ad_path)


def process_single_sample(
    hex_name: str,
    tissue_type: str,
    obfuscated_dir: Path,
    output_base: Path,
    bin_size: str = "008um",
    skip_medgemma: bool = False
) -> Dict:
    """
    Process single sample through full pipeline.

    Steps:
        1. Extract h5ad from binned_outputs.tar.gz
        2. Run uncertainty_spatial_analysis.py
        3. Run medgemma_v2_pipeline.py (optional)

    Returns:
        Results dictionary with paths and status
    """
    sample_id = hex_name.replace('.tar.gz', '')
    sample_dir = output_base / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    import hashlib
    sample_hash = hashlib.md5(sample_id.encode()).hexdigest()[:8]

    results = {
        'sample_hash': sample_hash,
        'status': 'STARTED',
        'timestamp': datetime.now().isoformat()
    }

    try:
        # Step 1: Extract h5ad
        print(f"[{sample_id}] Extracting h5ad from {hex_name}...")
        binned_tar = obfuscated_dir / hex_name

        h5ad_path = extract_binned_h5ad(
            str(binned_tar),
            sample_dir / "extracted",
            bin_size=bin_size
        )

        results['h5ad_path'] = h5ad_path
        print(f"[{sample_id}] Extracted: {h5ad_path}")

        # Step 2: Uncertainty spatial analysis
        print(f"[{sample_id}] Running uncertainty spatial analysis...")
        features_output = sample_dir / "uncertainty_features.json"

        from uncertainty_spatial_analysis import run_uncertainty_aware_spatial_analysis

        uncertainty_results = run_uncertainty_aware_spatial_analysis(
            adata_path=h5ad_path,
            output_path=str(features_output)
        )

        results['features_path'] = str(features_output)
        results['uncertainty_status'] = 'SUCCESS'
        print(f"[{sample_id}] Features saved: {features_output}")

        # Step 3: MedGemma report generation (optional)
        if not skip_medgemma:
            print(f"[{sample_id}] Running MedGemma report generation...")
            report_output = sample_dir / "medgemma_report.json"

            from medgemma_v2_pipeline import run_medgemma_v2_pipeline

            report = run_medgemma_v2_pipeline(
                uncertainty_features_path=str(features_output),
                output_path=str(report_output),
                tissue_type=map_tissue_type(tissue_type)
            )

            results['report_path'] = str(report_output)
            results['medgemma_status'] = 'SUCCESS'
            print(f"[{sample_id}] Report saved: {report_output}")

        results['status'] = 'COMPLETED'

    except Exception as e:
        results['status'] = 'FAILED'
        results['error'] = str(e)
        print(f"[{sample_id}] ERROR: {e}")

    # Save results
    results_path = sample_dir / "processing_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)

    return results


def map_tissue_type(raw_tissue: str) -> str:
    """Map obfuscated tissue names to reference cohort keys."""
    mapping = {
        'colon': 'colon_cancer',
        'lung': 'lung_cancer',
        'pancreas': 'pan_tissue',  # No pancreas-specific cohort
        'prostate': 'pan_tissue',
        'breast': 'breast_cancer',
        'brain': 'brain_glioma'
    }

    for key, value in mapping.items():
        if key in raw_tissue.lower():
            return value

    return 'pan_tissue'  # Default fallback


def run_parallel_pipeline(
    obfuscated_dir: str,
    mapping_file: str,
    output_dir: str,
    max_workers: int = None,
    bin_size: str = "008um",
    skip_medgemma: bool = False,
    tissue_filter: List[str] = None
):
    """
    Run pipeline on all obfuscated samples in parallel.

    Args:
        obfuscated_dir: Path to data/obfuscated
        mapping_file: Path to file_mapping.json
        output_dir: Base output directory
        max_workers: Number of parallel processes (auto-detected if None)
        bin_size: Resolution to process (008um, 016um, etc.)
        skip_medgemma: Skip MedGemma generation (faster testing)
        tissue_filter: Only process specific tissues (e.g., ['colon', 'lung'])
    """
    print("=" * 80)
    print("PARALLEL MEDGEMMA PIPELINE")
    print("=" * 80)

    # Load mapping
    print(f"\n[1/5] Loading file mapping from {mapping_file}")
    mapping = load_file_mapping(mapping_file)
    print(f"  Total files: {mapping['summary']['total_files']}")

    # Find samples
    print(f"\n[2/5] Finding binned_outputs.tar.gz files...")
    obf_dir = Path(obfuscated_dir)
    samples = find_binned_outputs(obf_dir, mapping)

    # Apply tissue filter
    if tissue_filter:
        samples = [s for s in samples if s[1] in tissue_filter]
        print(f"  Filtered to tissues: {tissue_filter}")

    print(f"  Found {len(samples)} samples to process")
    for hex_name, tissue, original in samples:
        print(f"    - {hex_name[:16]}... ({tissue})")

    # Determine workers
    if max_workers is None:
        max_workers = get_available_cores()

    print(f"\n[3/5] Starting parallel processing with {max_workers} workers...")
    print(f"  Bin size: {bin_size}")
    print(f"  MedGemma: {'DISABLED' if skip_medgemma else 'ENABLED'}")

    output_base = Path(output_dir)
    output_base.mkdir(parents=True, exist_ok=True)

    # Process samples
    all_results = []

    if max_workers == 1:
        # Serial processing (GPU or single-core)
        for hex_name, tissue, original in samples:
            result = process_single_sample(
                hex_name, tissue, obf_dir, output_base, bin_size, skip_medgemma
            )
            all_results.append(result)
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_single_sample,
                    hex_name, tissue, obf_dir, output_base, bin_size, skip_medgemma
                ): (hex_name, tissue)
                for hex_name, tissue, _ in samples
            }

            for future in as_completed(futures):
                hex_name, tissue = futures[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    print(f"  ✓ {hex_name[:16]}... ({tissue}): {result['status']}")
                except Exception as e:
                    print(f"  ✗ {hex_name[:16]}... ({tissue}): FAILED - {e}")
                    all_results.append({
                        'sample_id': hex_name,
                        'tissue': tissue,
                        'status': 'FAILED',
                        'error': str(e)
                    })

    # Summary
    print(f"\n[4/5] Processing complete!")
    completed = sum(1 for r in all_results if r['status'] == 'COMPLETED')
    failed = sum(1 for r in all_results if r['status'] == 'FAILED')

    print(f"  Completed: {completed}/{len(samples)}")
    print(f"  Failed: {failed}/{len(samples)}")

    # Save summary
    summary_path = output_base / "pipeline_summary.json"
    summary = {
        'timestamp': datetime.now().isoformat(),
        'config': {
            'max_workers': max_workers,
            'bin_size': bin_size,
            'skip_medgemma': skip_medgemma,
            'tissue_filter': tissue_filter
        },
        'results': all_results
    }

    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n[5/5] Summary saved to: {summary_path}")
    print("=" * 80)

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parallel MedGemma pipeline for obfuscated Visium HD data"
    )

    parser.add_argument(
        "--obfuscated-dir",
        default="data/obfuscated",
        help="Directory with obfuscated files"
    )

    parser.add_argument(
        "--mapping",
        default="data/file_mapping.json",
        help="File mapping JSON"
    )

    parser.add_argument(
        "--output",
        default="outputs/parallel_results",
        help="Output directory"
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (auto-detect if not specified)"
    )

    parser.add_argument(
        "--bin-size",
        default=None,
        type=str,
        help="Visium HD bin resolution (auto-detects highest if not specified)"
    )

    parser.add_argument(
        "--skip-medgemma",
        action="store_true",
        help="Skip MedGemma report generation (faster testing)"
    )

    parser.add_argument(
        "--tissues",
        nargs="+",
        default=None,
        help="Filter to specific tissues (e.g., colon lung)"
    )

    args = parser.parse_args()

    run_parallel_pipeline(
        obfuscated_dir=args.obfuscated_dir,
        mapping_file=args.mapping,
        output_dir=args.output,
        max_workers=args.workers,
        bin_size=args.bin_size,
        skip_medgemma=args.skip_medgemma,
        tissue_filter=args.tissues
    )

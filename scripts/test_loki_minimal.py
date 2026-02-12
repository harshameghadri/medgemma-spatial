"""
Minimal Loki Integration Test
Tests gene preprocessing and embedding pipeline WITHOUT pretrained model
"""

import sys
import os
import time
import psutil
import json
sys.path.insert(0, '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/Loki/src')

import numpy as np
import pandas as pd
import scanpy as sc

import loki.preprocess


def get_memory_usage_mb():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


def test_gene_preprocessing():
    """Test gene preprocessing functionality."""
    print("\n=== Loki Gene Preprocessing Test ===")
    start_time = time.time()
    start_mem = get_memory_usage_mb()

    try:
        adata_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/annotated_visium.h5ad'
        house_keeping_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/housekeeping_genes.csv'

        print(f"\nLoading data...")
        adata = sc.read_h5ad(adata_path)
        print(f"Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

        house_keeping_genes = pd.read_csv(house_keeping_path, index_col=0)
        print(f"Loaded {len(house_keeping_genes)} housekeeping genes")

        print(f"\nGenerating top-50 gene representations...")
        top_k_genes_str = loki.preprocess.generate_gene_df(
            adata,
            house_keeping_genes,
            todense=False
        )

        print(f"\nResults:")
        print(f"  Generated gene strings for {len(top_k_genes_str)} spots")
        print(f"\n  Sample spot gene representation (first 150 chars):")
        print(f"  {top_k_genes_str['label'].iloc[0][:150]}...")

        elapsed = time.time() - start_time
        memory_delta = get_memory_usage_mb() - start_mem

        print(f"\nPerformance:")
        print(f"  Time: {elapsed:.2f} seconds")
        print(f"  Memory delta: {memory_delta:.1f} MB")
        print(f"  Current memory: {get_memory_usage_mb():.1f} MB")

        stats = {
            'n_spots': len(top_k_genes_str),
            'avg_genes_per_spot': top_k_genes_str['label'].str.split().str.len().mean(),
            'sample_gene_string': top_k_genes_str['label'].iloc[0][:200]
        }

        result = {
            'status': 'PASS',
            'elapsed_time': elapsed,
            'memory_delta_mb': memory_delta,
            'current_memory_mb': get_memory_usage_mb(),
            'stats': stats
        }

        print(f"\n✓ Gene preprocessing PASSED")
        return result

    except Exception as e:
        import traceback
        print(f"\n✗ Gene preprocessing FAILED")
        print(f"Error: {str(e)}")
        print(f"\nTraceback:")
        print(traceback.format_exc())

        return {
            'status': 'FAIL',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'elapsed_time': time.time() - start_time
        }


def main():
    """Run minimal Loki test."""
    print("=" * 70)
    print("LOKI MINIMAL INTEGRATION TEST")
    print("=" * 70)
    print(f"Test: Gene preprocessing (no model required)")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    result = test_gene_preprocessing()

    output_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/loki_minimal_test_results.json'
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    if result['status'] == 'PASS':
        print("\n" + "=" * 70)
        print("NEXT STEPS:")
        print("=" * 70)
        print("1. Download Loki pretrained weights from HuggingFace:")
        print("   https://huggingface.co/WangGuangyuLab/Loki")
        print("2. Place checkpoint.pt in:")
        print("   /Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/loki_checkpoint.pt")
        print("3. Run full integration test:")
        print("   python scripts/test_loki_integration.py")
        return 0
    else:
        print("\n✗ Test failed - check error details above")
        return 1


if __name__ == "__main__":
    sys.exit(main())

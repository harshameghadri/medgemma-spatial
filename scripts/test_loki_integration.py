"""
Loki Integration Test Script
Tests Loki spatial foundation model on validated Visium sample
Time limit: 2 hours
"""

import sys
import os
import time
import psutil
import json
import traceback
sys.path.insert(0, '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/Loki/src')

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from PIL import Image

# PyTorch 2.6 compatibility: override open_clip's checkpoint loading to use weights_only=False
# The Loki checkpoint was saved with older PyTorch and contains numpy globals
import open_clip.factory as _ocf
_orig_load_state_dict = _ocf.load_state_dict
def _patched_load_state_dict(checkpoint_path, device=None, weights_only=True):
    return _orig_load_state_dict(checkpoint_path, device=device, weights_only=False)
_ocf.load_state_dict = _patched_load_state_dict

import loki.utils
import loki.preprocess


def get_memory_usage_mb():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


def test_loki_model_loading():
    """Test 1: Load Loki OmiCLIP model."""
    print("\n=== TEST 1: Model Loading ===")
    start_time = time.time()

    try:
        model_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/checkpoint.pt'
        device = 'cpu'

        if not os.path.exists(model_path):
            return {
                'status': 'SKIP',
                'reason': 'Model checkpoint not found. Download from HuggingFace.',
                'download_url': 'https://huggingface.co/WangGuangyuLab/Loki',
                'elapsed_time': time.time() - start_time
            }

        model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
        model.eval()

        memory_mb = get_memory_usage_mb()
        elapsed = time.time() - start_time

        return {
            'status': 'PASS',
            'elapsed_time': elapsed,
            'memory_mb': memory_mb,
            'device': device
        }

    except Exception as e:
        return {
            'status': 'FAIL',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'elapsed_time': time.time() - start_time
        }


def test_generate_gene_embeddings(adata):
    """Test 2: Generate text embeddings from gene expression."""
    print("\n=== TEST 2: Gene Text Embeddings ===")
    start_time = time.time()

    try:
        house_keeping_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/housekeeping_genes.csv'

        if not os.path.exists(house_keeping_path):
            return {
                'status': 'SKIP',
                'reason': 'Housekeeping genes file not found',
                'elapsed_time': time.time() - start_time
            }

        house_keeping_genes = pd.read_csv(house_keeping_path, index_col=0)

        import scipy.sparse as sp
        is_sparse = sp.issparse(adata.X)
        top_k_genes_str = loki.preprocess.generate_gene_df(
            adata,
            house_keeping_genes,
            todense=is_sparse
        )

        elapsed = time.time() - start_time
        memory_mb = get_memory_usage_mb()

        return {
            'status': 'PASS',
            'n_spots': len(top_k_genes_str),
            'sample_genes': top_k_genes_str['label'].iloc[0][:100],
            'elapsed_time': elapsed,
            'memory_mb': memory_mb
        }

    except Exception as e:
        return {
            'status': 'FAIL',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'elapsed_time': time.time() - start_time
        }


def test_encode_spatial_spots(adata, model, tokenizer, device='cpu'):
    """Test 3: Encode spatial spots with Loki."""
    print("\n=== TEST 3: Encode Spatial Spots ===")
    start_time = time.time()

    try:
        house_keeping_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/housekeeping_genes.csv'
        house_keeping_genes = pd.read_csv(house_keeping_path, index_col=0)

        import scipy.sparse as sp
        is_sparse = sp.issparse(adata.X)
        top_k_genes_str = loki.preprocess.generate_gene_df(
            adata,
            house_keeping_genes,
            todense=is_sparse
        )

        text_embeddings = loki.utils.encode_text_df(
            model,
            tokenizer,
            top_k_genes_str,
            'label',
            device
        )

        elapsed = time.time() - start_time
        memory_mb = get_memory_usage_mb()

        return {
            'status': 'PASS',
            'embedding_shape': list(text_embeddings.shape),
            'embedding_dim': text_embeddings.shape[1],
            'mean_embedding': text_embeddings.mean(dim=0).cpu().numpy().tolist()[:10],
            'elapsed_time': elapsed,
            'memory_mb': memory_mb
        }

    except Exception as e:
        return {
            'status': 'FAIL',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'elapsed_time': time.time() - start_time
        }


def test_extract_loki_features(adata, embeddings):
    """Test 4: Extract Loki features for pipeline integration."""
    print("\n=== TEST 4: Feature Extraction ===")
    start_time = time.time()

    try:
        embeddings_np = embeddings.cpu().numpy()

        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=5, random_state=42)
        loki_clusters = kmeans.fit_predict(embeddings_np)

        cluster_counts = pd.Series(loki_clusters).value_counts().to_dict()

        from scipy.spatial.distance import pdist, squareform
        distances = squareform(pdist(embeddings_np, metric='euclidean'))
        avg_intra_cluster_dist = {}

        for cluster_id in range(5):
            mask = loki_clusters == cluster_id
            if mask.sum() > 1:
                cluster_distances = distances[np.ix_(mask, mask)]
                avg_intra_cluster_dist[f'cluster_{cluster_id}'] = float(
                    cluster_distances[np.triu_indices_from(cluster_distances, k=1)].mean()
                )

        features = {
            'loki_embeddings': {
                'embedding_dim': int(embeddings.shape[1]),
                'n_spots': int(embeddings.shape[0]),
                'cluster_counts': {f'cluster_{k}': int(v) for k, v in cluster_counts.items()},
                'avg_intra_cluster_distance': avg_intra_cluster_dist,
                'embedding_mean_norm': float(np.linalg.norm(embeddings_np.mean(axis=0)))
            }
        }

        elapsed = time.time() - start_time
        memory_mb = get_memory_usage_mb()

        return {
            'status': 'PASS',
            'features': features,
            'elapsed_time': elapsed,
            'memory_mb': memory_mb
        }

    except Exception as e:
        return {
            'status': 'FAIL',
            'error': str(e),
            'traceback': traceback.format_exc(),
            'elapsed_time': time.time() - start_time
        }


def run_loki_tests():
    """Main test runner."""
    print("=" * 60)
    print("LOKI INTEGRATION TEST")
    print("=" * 60)
    print(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Device: CPU (M1 Mac)")

    overall_start = time.time()
    results = {}

    adata_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/annotated_visium.h5ad'
    print(f"\nLoading sample data: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    results['test_1_model_loading'] = test_loki_model_loading()
    print(f"Status: {results['test_1_model_loading']['status']}")

    if results['test_1_model_loading']['status'] == 'PASS':
        model_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/checkpoint.pt'
        device = 'cpu'
        model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
        model.eval()

        results['test_2_gene_embeddings'] = test_generate_gene_embeddings(adata)
        print(f"Status: {results['test_2_gene_embeddings']['status']}")

        if results['test_2_gene_embeddings']['status'] == 'PASS':
            results['test_3_encode_spots'] = test_encode_spatial_spots(adata, model, tokenizer, device)
            print(f"Status: {results['test_3_encode_spots']['status']}")

            if results['test_3_encode_spots']['status'] == 'PASS':
                house_keeping_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/data/housekeeping_genes.csv'
                house_keeping_genes = pd.read_csv(house_keeping_path, index_col=0)
                import scipy.sparse as sp
                top_k_genes_str = loki.preprocess.generate_gene_df(adata, house_keeping_genes, todense=sp.issparse(adata.X))
                embeddings = loki.utils.encode_text_df(model, tokenizer, top_k_genes_str, 'label', device)

                results['test_4_feature_extraction'] = test_extract_loki_features(adata, embeddings)
                print(f"Status: {results['test_4_feature_extraction']['status']}")
    else:
        results['test_2_gene_embeddings'] = {'status': 'SKIPPED', 'reason': 'Model loading failed'}
        results['test_3_encode_spots'] = {'status': 'SKIPPED', 'reason': 'Model loading failed'}
        results['test_4_feature_extraction'] = {'status': 'SKIPPED', 'reason': 'Model loading failed'}

    total_elapsed = time.time() - overall_start

    results['summary'] = {
        'total_elapsed_time': total_elapsed,
        'total_elapsed_formatted': f"{total_elapsed/60:.2f} minutes",
        'final_memory_mb': get_memory_usage_mb(),
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
    }

    output_path = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/loki_test_results.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 60)
    print("TEST RESULTS SUMMARY")
    print("=" * 60)
    for test_name, result in results.items():
        if test_name != 'summary':
            status = result.get('status', 'UNKNOWN')
            print(f"{test_name}: {status}")
    print(f"\nTotal Time: {total_elapsed/60:.2f} minutes")
    print(f"Final Memory: {get_memory_usage_mb():.1f} MB")
    print(f"\nResults saved to: {output_path}")

    return results


if __name__ == "__main__":
    results = run_loki_tests()

    all_passed = all(
        r.get('status') == 'PASS'
        for k, r in results.items()
        if k != 'summary' and r.get('status') != 'SKIPPED'
    )

    if all_passed:
        print("\n✓ ALL TESTS PASSED - Loki integration feasible")
        sys.exit(0)
    else:
        print("\n✗ SOME TESTS FAILED - Review results for details")
        sys.exit(1)

#!/usr/bin/env python3
"""
Robustness Test #3: Full spatial analysis pipeline test
Tests annotation, heterogeneity analysis, feature extraction, and anti-parroting
"""

import sys
import time
import json
import traceback
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import scanpy as sc

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from notebooks.marker_annotation import (
    load_panglaodb_markers,
    get_tissue_specific_markers,
    score_clusters_for_markers
)
from notebooks.uncertainty_spatial_analysis import (
    detect_doublets_scrublet,
    annotate_spatial_regions,
    assess_annotation_confidence,
    compute_morans_i_with_uncertainty,
    compute_spatial_entropy_with_bootstrapping
)


def validate_h5ad_structure(adata) -> dict:
    """Validate H5AD file has required fields."""
    checks = {
        'has_X_matrix': hasattr(adata, 'X') and adata.X is not None,
        'has_obs': hasattr(adata, 'obs') and len(adata.obs) > 0,
        'has_var': hasattr(adata, 'var') and len(adata.var) > 0,
        'n_obs': int(adata.n_obs),
        'n_vars': int(adata.n_vars),
        'X_shape': list(adata.X.shape) if hasattr(adata, 'X') else None,
        'obs_columns': list(adata.obs.columns),
        'var_columns': list(adata.var.columns),
    }

    # Check for spatial coordinates
    has_spatial = False
    spatial_key = None
    for key in ['spatial', 'X_spatial']:
        if key in adata.obsm:
            has_spatial = True
            spatial_key = key
            break

    checks['has_spatial_coords'] = has_spatial
    checks['spatial_key'] = spatial_key

    return checks


def test_spatial_coordinates(adata) -> dict:
    """Test spatial coordinate extraction."""
    print("\n=== Testing Spatial Coordinate Extraction ===")
    start_time = time.time()

    try:
        # Extract spatial coordinates from obsm
        has_spatial = False
        spatial_key = None

        for key in ['spatial', 'X_spatial']:
            if key in adata.obsm:
                has_spatial = True
                spatial_key = key
                coords = adata.obsm[key]
                adata.obs['spatial_x'] = coords[:, 0]
                adata.obs['spatial_y'] = coords[:, 1]
                break

        if not has_spatial:
            raise ValueError("No spatial coordinates found in adata.obsm")

        spatial_metrics = {
            'has_spatial': True,
            'spatial_key': spatial_key,
            'n_spots': int(adata.n_obs),
            'coord_range_x': [float(adata.obs['spatial_x'].min()), float(adata.obs['spatial_x'].max())],
            'coord_range_y': [float(adata.obs['spatial_y'].min()), float(adata.obs['spatial_y'].max())]
        }

        result = {
            'status': 'SUCCESS',
            'execution_time': time.time() - start_time,
            'metrics': spatial_metrics,
            'has_coords': True,
            'n_spots_with_coords': int(adata.n_obs)
        }
    except Exception as e:
        result = {
            'status': 'FAILED',
            'execution_time': time.time() - start_time,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

    return result


def test_cell_type_annotation(adata) -> dict:
    """Test cell type annotation with marker genes."""
    print("\n=== Testing Cell Type Annotation ===")
    start_time = time.time()

    try:
        adata, annotation_metrics = annotate_spatial_regions(
            adata,
            resolution=0.5,
            use_markers=True,
            marker_file="data/PanglaoDB_markers_27_Mar_2020.tsv",
            tissue=None  # Let it auto-detect
        )

        # Verify annotation results
        has_cell_types = 'cell_type' in adata.obs
        n_cell_types = adata.obs['cell_type'].nunique() if has_cell_types else 0

        result = {
            'status': 'SUCCESS',
            'execution_time': time.time() - start_time,
            'metrics': annotation_metrics,
            'has_cell_types': has_cell_types,
            'n_cell_types': int(n_cell_types),
            'cell_type_distribution': adata.obs['cell_type'].value_counts().to_dict() if has_cell_types else {}
        }
    except Exception as e:
        result = {
            'status': 'FAILED',
            'execution_time': time.time() - start_time,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

    return result


def test_spatial_heterogeneity(adata) -> dict:
    """Test spatial heterogeneity analysis."""
    print("\n=== Testing Spatial Heterogeneity Analysis ===")
    start_time = time.time()

    try:
        # Compute Moran's I with uncertainty
        morans_results = compute_morans_i_with_uncertainty(adata, n_genes=50, n_permutations=100)

        # Compute spatial entropy
        entropy_results = compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=100)

        heterogeneity_metrics = {
            'morans_i': morans_results,
            'spatial_entropy': entropy_results
        }

        result = {
            'status': 'SUCCESS',
            'execution_time': time.time() - start_time,
            'metrics': heterogeneity_metrics,
            'has_morans_i': True,
            'has_diversity': True
        }
    except Exception as e:
        result = {
            'status': 'FAILED',
            'execution_time': time.time() - start_time,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

    return result


def generate_features_json(adata, spatial_metrics, annotation_metrics, heterogeneity_metrics) -> dict:
    """Generate features JSON for report generation."""
    print("\n=== Generating Features JSON ===")
    start_time = time.time()

    try:
        features = {
            'sample_info': {
                'n_spots': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'has_spatial': bool(spatial_metrics.get('has_spatial', False))
            },
            'spatial_features': spatial_metrics,
            'annotation': {
                'n_cell_types': int(annotation_metrics.get('n_cell_types', 0)),
                'cell_type_counts': annotation_metrics.get('cell_type_counts', {}),
                'annotation_quality': annotation_metrics.get('annotation_quality', 'unknown')
            },
            'heterogeneity': heterogeneity_metrics,
            'timestamp': datetime.now().isoformat()
        }

        result = {
            'status': 'SUCCESS',
            'execution_time': time.time() - start_time,
            'features': features,
            'json_size_bytes': len(json.dumps(features))
        }
    except Exception as e:
        result = {
            'status': 'FAILED',
            'execution_time': time.time() - start_time,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

    return result


def test_anti_parroting_prompt(features_json) -> dict:
    """Test anti-parroting prompt generation."""
    print("\n=== Testing Anti-Parroting Prompt Generation ===")
    start_time = time.time()

    try:
        # Check for tissue type leakage
        features_str = json.dumps(features_json).lower()

        # Common tissue types that should NOT appear in generic features
        tissue_keywords = ['colon', 'intestine', 'brain', 'breast', 'lung', 'liver', 'kidney']
        tissue_leakage = [kw for kw in tissue_keywords if kw in features_str]

        # Generate anti-parroting prompt
        prompt = f"""You are a spatial transcriptomics expert analyzing tissue architecture.

SPATIAL DATA SUMMARY:
- Total spots: {features_json['sample_info']['n_spots']}
- Cell types identified: {features_json['annotation']['n_cell_types']}
- Spatial autocorrelation (Moran's I): {features_json['heterogeneity'].get('morans_i', 'N/A')}
- Shannon diversity: {features_json['heterogeneity'].get('shannon_diversity', 'N/A')}

IMPORTANT INSTRUCTIONS:
1. Describe the spatial organization patterns you observe
2. Discuss cellular diversity and mixing
3. Interpret Moran's I (spatial autocorrelation) values
4. DO NOT mention specific tissue types (colon, intestine, etc.)
5. DO NOT copy exact numbers from the data
6. Focus on biological interpretation and spatial patterns

Generate a 150-200 word pathology-style report describing the spatial architecture."""

        result = {
            'status': 'SUCCESS',
            'execution_time': time.time() - start_time,
            'prompt_length': len(prompt),
            'has_tissue_leakage': len(tissue_leakage) > 0,
            'tissue_keywords_found': tissue_leakage,
            'prompt_preview': prompt[:200] + '...'
        }
    except Exception as e:
        result = {
            'status': 'FAILED',
            'execution_time': time.time() - start_time,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

    return result


def run_robustness_test(h5ad_path: str, output_path: str):
    """Run complete robustness test pipeline."""
    print("="*80)
    print("ROBUSTNESS TEST AGENT #3")
    print("="*80)
    print(f"Test file: {h5ad_path}")
    print(f"Output: {output_path}")
    print(f"Start time: {datetime.now().isoformat()}")

    test_results = {
        'test_info': {
            'test_name': 'Robustness Test #3',
            'h5ad_path': h5ad_path,
            'start_time': datetime.now().isoformat(),
        },
        'steps': {}
    }

    # Step 1: Load H5AD file
    print("\n=== Step 1: Loading H5AD File ===")
    step1_start = time.time()
    try:
        adata = sc.read_h5ad(h5ad_path)
        file_size_mb = Path(h5ad_path).stat().st_size / (1024**2)

        test_results['steps']['1_load_file'] = {
            'status': 'SUCCESS',
            'execution_time': time.time() - step1_start,
            'file_size_mb': round(file_size_mb, 2),
            'validation': validate_h5ad_structure(adata)
        }
        print(f"✓ Loaded successfully: {adata.n_obs} spots × {adata.n_vars} genes ({file_size_mb:.1f} MB)")
    except Exception as e:
        test_results['steps']['1_load_file'] = {
            'status': 'FAILED',
            'execution_time': time.time() - step1_start,
            'error': str(e),
            'traceback': traceback.format_exc()
        }
        print(f"✗ Failed to load: {e}")
        save_results(test_results, output_path)
        return

    # Step 2: Extract spatial coordinates
    test_results['steps']['2_spatial_coords'] = test_spatial_coordinates(adata)
    status = test_results['steps']['2_spatial_coords']['status']
    print(f"{('✓' if status == 'SUCCESS' else '✗')} Spatial coordinates: {status}")

    # Step 3: Cell type annotation
    test_results['steps']['3_cell_type_annotation'] = test_cell_type_annotation(adata)
    status = test_results['steps']['3_cell_type_annotation']['status']
    print(f"{('✓' if status == 'SUCCESS' else '✗')} Cell type annotation: {status}")

    if test_results['steps']['3_cell_type_annotation']['status'] != 'SUCCESS':
        print("⚠ Annotation failed, skipping heterogeneity analysis")
        test_results['test_info']['end_time'] = datetime.now().isoformat()
        save_results(test_results, output_path)
        return

    # Step 4: Spatial heterogeneity
    test_results['steps']['4_spatial_heterogeneity'] = test_spatial_heterogeneity(adata)
    status = test_results['steps']['4_spatial_heterogeneity']['status']
    print(f"{('✓' if status == 'SUCCESS' else '✗')} Spatial heterogeneity: {status}")

    # Step 5: Generate features JSON
    spatial_metrics = test_results['steps']['2_spatial_coords'].get('metrics', {})
    annotation_metrics = test_results['steps']['3_cell_type_annotation'].get('metrics', {})
    heterogeneity_metrics = test_results['steps']['4_spatial_heterogeneity'].get('metrics', {})

    test_results['steps']['5_features_json'] = generate_features_json(
        adata, spatial_metrics, annotation_metrics, heterogeneity_metrics
    )
    status = test_results['steps']['5_features_json']['status']
    print(f"{('✓' if status == 'SUCCESS' else '✗')} Features JSON: {status}")

    # Step 6: Anti-parroting prompt
    if test_results['steps']['5_features_json']['status'] == 'SUCCESS':
        features_json = test_results['steps']['5_features_json']['features']
        test_results['steps']['6_anti_parroting'] = test_anti_parroting_prompt(features_json)
        status = test_results['steps']['6_anti_parroting']['status']
        print(f"{('✓' if status == 'SUCCESS' else '✗')} Anti-parroting prompt: {status}")

    # Summary
    test_results['test_info']['end_time'] = datetime.now().isoformat()
    total_time = sum(
        step.get('execution_time', 0)
        for step in test_results['steps'].values()
    )
    test_results['test_info']['total_execution_time'] = round(total_time, 2)

    # Quality metrics
    success_count = sum(
        1 for step in test_results['steps'].values()
        if step.get('status') == 'SUCCESS'
    )
    total_steps = len(test_results['steps'])

    test_results['summary'] = {
        'total_steps': total_steps,
        'successful_steps': success_count,
        'failed_steps': total_steps - success_count,
        'success_rate': round(success_count / total_steps * 100, 1),
        'total_time_seconds': round(total_time, 2),
        'meets_time_constraint': total_time < 600  # <10 minutes
    }

    # Save results
    save_results(test_results, output_path)

    # Print summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    print(f"Steps completed: {success_count}/{total_steps} ({test_results['summary']['success_rate']}%)")
    print(f"Total execution time: {test_results['summary']['total_time_seconds']}s")
    print(f"Time constraint (<10min): {'✓ PASS' if test_results['summary']['meets_time_constraint'] else '✗ FAIL'}")
    print(f"Results saved to: {output_path}")


def save_results(results: dict, output_path: str):
    """Save test results to JSON file."""
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n✓ Results saved to: {output_path}")


if __name__ == "__main__":
    h5ad_path = "/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/test_full_markers_20260202_063705/e2735493904baddea063a55d5e676d24/extracted/square_008um.h5ad"
    output_path = "/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/robustness_test_sample3.json"

    run_robustness_test(h5ad_path, output_path)

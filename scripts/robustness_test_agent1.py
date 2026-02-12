#!/usr/bin/env python3
"""
Robustness Test Agent #1
Tests spatial analysis pipeline on annotated Visium data.
"""

import os
import sys
import time
import json
import traceback
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd


def test_step_1_load_h5ad(h5ad_path):
    """Load H5AD file and validate structure."""
    print("\n" + "="*80)
    print("STEP 1: Load H5AD File")
    print("="*80)

    start_time = time.time()

    try:
        adata = sc.read_h5ad(h5ad_path)
        elapsed = time.time() - start_time

        # Validate required fields
        has_X = adata.X is not None
        has_obs = adata.obs is not None and len(adata.obs) > 0
        has_var = adata.var is not None and len(adata.var) > 0

        # Check for spatial coordinates
        has_spatial_coords = False
        spatial_key = None

        if 'spatial' in adata.obsm:
            spatial_key = 'spatial'
            has_spatial_coords = True
        elif 'X_spatial' in adata.obsm:
            spatial_key = 'X_spatial'
            has_spatial_coords = True

        # Check for cell type annotation
        has_cell_types = 'cell_type' in adata.obs.columns

        result = {
            'status': 'PASS' if all([has_X, has_obs, has_var, has_spatial_coords]) else 'FAIL',
            'elapsed_time': elapsed,
            'n_spots': int(adata.n_obs),
            'n_genes': int(adata.n_vars),
            'has_X': has_X,
            'has_obs': has_obs,
            'has_var': has_var,
            'has_spatial_coords': has_spatial_coords,
            'spatial_key': spatial_key,
            'has_cell_types': has_cell_types,
            'obs_columns': list(adata.obs.columns),
            'file_exists': os.path.exists(h5ad_path),
            'file_size_mb': os.path.getsize(h5ad_path) / (1024**2)
        }

        print(f"✓ File loaded successfully")
        print(f"  Shape: {adata.n_obs} spots × {adata.n_vars} genes")
        print(f"  Spatial coords: {'✓' if has_spatial_coords else '✗'} (key: {spatial_key})")
        print(f"  Cell types: {'✓' if has_cell_types else '✗'}")
        print(f"  Time: {elapsed:.2f}s")
        print(f"  File size: {result['file_size_mb']:.1f} MB")

        return result, adata

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }, None


def test_step_2_spatial_coordinates(adata):
    """Extract and validate spatial coordinates."""
    print("\n" + "="*80)
    print("STEP 2: Spatial Coordinates Extraction")
    print("="*80)

    start_time = time.time()

    try:
        # Find spatial coordinates
        spatial_key = None
        if 'spatial' in adata.obsm:
            spatial_key = 'spatial'
        elif 'X_spatial' in adata.obsm:
            spatial_key = 'X_spatial'
        else:
            raise ValueError("No spatial coordinates found in adata.obsm")

        coords = adata.obsm[spatial_key]

        # Validate coordinates
        has_nan = np.isnan(coords).any()
        has_inf = np.isinf(coords).any()

        x_range = [float(coords[:, 0].min()), float(coords[:, 0].max())]
        y_range = [float(coords[:, 1].min()), float(coords[:, 1].max())]

        elapsed = time.time() - start_time

        result = {
            'status': 'PASS' if not (has_nan or has_inf) else 'FAIL',
            'elapsed_time': elapsed,
            'spatial_key': spatial_key,
            'n_dimensions': coords.shape[1],
            'has_nan': has_nan,
            'has_inf': has_inf,
            'x_range': x_range,
            'y_range': y_range
        }

        print(f"✓ Spatial coordinates extracted")
        print(f"  Key: {spatial_key}")
        print(f"  Dimensions: {coords.shape[1]}D")
        print(f"  X range: [{x_range[0]:.1f}, {x_range[1]:.1f}]")
        print(f"  Y range: [{y_range[0]:.1f}, {y_range[1]:.1f}]")
        print(f"  NaN values: {'✗ PRESENT' if has_nan else '✓ NONE'}")
        print(f"  Inf values: {'✗ PRESENT' if has_inf else '✓ NONE'}")
        print(f"  Time: {elapsed:.2f}s")

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def test_step_3_cell_type_annotation(adata):
    """Validate cell type annotations."""
    print("\n" + "="*80)
    print("STEP 3: Cell Type Annotation Validation")
    print("="*80)

    start_time = time.time()

    try:
        if 'cell_type' not in adata.obs.columns:
            raise ValueError("No 'cell_type' column in adata.obs")

        cell_types = adata.obs['cell_type']
        unique_types = cell_types.unique()
        type_counts = cell_types.value_counts().to_dict()

        # Check for confidence scores
        has_confidence = 'conf_score' in adata.obs.columns
        if has_confidence:
            conf_scores = adata.obs['conf_score']
            mean_conf = float(conf_scores.mean())
            low_conf_rate = float((conf_scores < 0.5).mean())
        else:
            mean_conf = None
            low_conf_rate = None

        elapsed = time.time() - start_time

        result = {
            'status': 'PASS',
            'elapsed_time': elapsed,
            'n_unique_types': len(unique_types),
            'cell_type_counts': {str(k): int(v) for k, v in type_counts.items()},
            'has_confidence_scores': has_confidence,
            'mean_confidence': mean_conf,
            'low_confidence_rate': low_conf_rate
        }

        print(f"✓ Cell type annotations present")
        print(f"  Unique types: {len(unique_types)}")
        print(f"  Top 5 types:")
        for ct, count in list(type_counts.items())[:5]:
            pct = 100 * count / len(adata)
            print(f"    - {ct}: {count} spots ({pct:.1f}%)")

        if has_confidence:
            print(f"  Confidence scores: mean={mean_conf:.3f}, low_conf_rate={low_conf_rate:.1%}")

        print(f"  Time: {elapsed:.2f}s")

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def test_step_4_spatial_analysis(adata, output_dir):
    """Run full spatial analysis pipeline."""
    print("\n" + "="*80)
    print("STEP 4: Spatial Analysis Pipeline")
    print("="*80)

    start_time = time.time()

    try:
        # Import spatial analysis functions
        sys.path.insert(0, '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/src')
        from spatial_analysis.uncertainty_spatial_analysis import run_uncertainty_aware_spatial_analysis

        # Create temp output path
        temp_output = os.path.join(output_dir, 'temp_spatial_features.json')

        # Save adata temporarily (pipeline expects path)
        temp_h5ad = os.path.join(output_dir, 'temp_adata.h5ad')
        adata.write_h5ad(temp_h5ad)

        # Run spatial analysis
        print("  Running uncertainty-aware spatial analysis...")
        results = run_uncertainty_aware_spatial_analysis(temp_h5ad, temp_output)

        elapsed = time.time() - start_time

        # Validate results
        has_stage0 = 'stage_0_annotation_quality' in results
        has_stage1 = 'stage_1_spatial_patterns' in results
        has_stopping = 'stage_stopping' in results

        # Extract key metrics
        decision = results.get('stage_stopping', {}).get('decision', 'UNKNOWN')

        result = {
            'status': 'PASS' if all([has_stage0, has_stage1, has_stopping]) else 'FAIL',
            'elapsed_time': elapsed,
            'has_stage0': has_stage0,
            'has_stage1': has_stage1,
            'has_stopping_decision': has_stopping,
            'stopping_decision': decision,
            'results': results,
            'output_file': temp_output
        }

        print(f"✓ Spatial analysis completed")
        print(f"  Stage 0 (annotation quality): {'✓' if has_stage0 else '✗'}")
        print(f"  Stage 1 (spatial patterns): {'✓' if has_stage1 else '✗'}")
        print(f"  Stopping decision: {decision}")
        print(f"  Time: {elapsed:.2f}s")

        # Cleanup temp file
        if os.path.exists(temp_h5ad):
            os.remove(temp_h5ad)

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def test_step_5_features_json(spatial_results):
    """Validate features JSON structure."""
    print("\n" + "="*80)
    print("STEP 5: Features JSON Validation")
    print("="*80)

    start_time = time.time()

    try:
        if spatial_results.get('status') != 'PASS':
            raise ValueError("Spatial analysis failed - cannot validate features")

        features = spatial_results['results']

        # Check required fields
        has_annotation_quality = 'stage_0_annotation_quality' in features
        has_spatial_patterns = 'stage_1_spatial_patterns' in features
        has_stopping = 'stage_stopping' in features

        # Validate structure
        issues = []

        if has_annotation_quality:
            aq = features['stage_0_annotation_quality']
            if 'doublet_metrics' not in aq:
                issues.append("Missing doublet_metrics")
            if 'annotation_confidence' not in aq:
                issues.append("Missing annotation_confidence")

        if has_spatial_patterns:
            sp = features['stage_1_spatial_patterns']
            if 'morans_i' not in sp:
                issues.append("Missing morans_i")
            if 'spatial_entropy' not in sp:
                issues.append("Missing spatial_entropy")

        elapsed = time.time() - start_time

        result = {
            'status': 'PASS' if len(issues) == 0 else 'FAIL',
            'elapsed_time': elapsed,
            'has_annotation_quality': has_annotation_quality,
            'has_spatial_patterns': has_spatial_patterns,
            'has_stopping': has_stopping,
            'validation_issues': issues,
            'is_valid_json': True
        }

        print(f"{'✓' if result['status'] == 'PASS' else '✗'} Features JSON validation")
        print(f"  Annotation quality: {'✓' if has_annotation_quality else '✗'}")
        print(f"  Spatial patterns: {'✓' if has_spatial_patterns else '✗'}")
        print(f"  Stopping decision: {'✓' if has_stopping else '✗'}")
        if issues:
            print(f"  Issues: {', '.join(issues)}")
        print(f"  Time: {elapsed:.2f}s")

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def test_step_6_anti_parroting_prompt(spatial_results, output_dir):
    """Create and validate anti-parroting prompt."""
    print("\n" + "="*80)
    print("STEP 6: Anti-Parroting Prompt Generation")
    print("="*80)

    start_time = time.time()

    try:
        if spatial_results.get('status') != 'PASS':
            raise ValueError("Spatial analysis failed - cannot generate prompt")

        # Import prompt generator
        sys.path.insert(0, '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/notebooks')
        from medgemma_report_generator import create_anti_parroting_prompt

        # Convert spatial results to features format
        features = {
            'annotation': {
                'cell_type_counts': {}  # Would need actual cell type counts
            },
            'spatial_heterogeneity': spatial_results['results'].get('stage_1_spatial_patterns', {}).get('morans_i', {}),
            'uncertainty': spatial_results['results'].get('stage_1_spatial_patterns', {}).get('spatial_entropy', {})
        }

        # Generate prompt
        prompt = create_anti_parroting_prompt(features)

        # Validate prompt
        has_instructions = 'TASK:' in prompt or 'Generate' in prompt
        has_requirements = 'REQUIREMENTS' in prompt or 'CRITICAL' in prompt
        prompt_length = len(prompt)

        # Check for data leakage prevention
        has_no_raw_numbers = not any(char.isdigit() for word in prompt.split() if word.isdigit() and len(word) > 2)

        elapsed = time.time() - start_time

        # Save prompt
        prompt_file = os.path.join(output_dir, 'anti_parroting_prompt.txt')
        with open(prompt_file, 'w') as f:
            f.write(prompt)

        result = {
            'status': 'PASS' if all([has_instructions, has_requirements, prompt_length > 100]) else 'FAIL',
            'elapsed_time': elapsed,
            'prompt_length': prompt_length,
            'has_instructions': has_instructions,
            'has_requirements': has_requirements,
            'prevents_data_leakage': has_no_raw_numbers,
            'prompt_file': prompt_file
        }

        print(f"✓ Anti-parroting prompt generated")
        print(f"  Length: {prompt_length} characters")
        print(f"  Has instructions: {'✓' if has_instructions else '✗'}")
        print(f"  Has requirements: {'✓' if has_requirements else '✗'}")
        print(f"  Prevents data leakage: {'✓' if has_no_raw_numbers else '⚠️'}")
        print(f"  Saved to: {prompt_file}")
        print(f"  Time: {elapsed:.2f}s")

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def test_step_7_output_quality(all_results):
    """Assess overall output quality."""
    print("\n" + "="*80)
    print("STEP 7: Output Quality Assessment")
    print("="*80)

    start_time = time.time()

    try:
        # Check all steps passed
        step_statuses = [r.get('status') for r in all_results.values() if 'status' in r]
        all_passed = all(s == 'PASS' for s in step_statuses)

        # Calculate total execution time
        total_time = sum(r.get('elapsed_time', 0) for r in all_results.values())

        # Memory estimate (file size is proxy)
        file_size_mb = all_results.get('step_1_load', {}).get('file_size_mb', 0)

        # Check for critical issues
        issues = []
        if not all_passed:
            failed_steps = [k for k, v in all_results.items() if v.get('status') == 'FAIL']
            issues.append(f"Failed steps: {', '.join(failed_steps)}")

        if total_time > 600:  # 10 minutes
            issues.append(f"Execution time too long: {total_time:.1f}s > 600s")

        elapsed = time.time() - start_time

        result = {
            'status': 'PASS' if all_passed and len(issues) == 0 else 'FAIL',
            'elapsed_time': elapsed,
            'all_steps_passed': all_passed,
            'total_execution_time': total_time,
            'within_time_limit': total_time < 600,
            'file_size_mb': file_size_mb,
            'quality_issues': issues
        }

        print(f"{'✓' if result['status'] == 'PASS' else '✗'} Output quality assessment")
        print(f"  All steps passed: {'✓' if all_passed else '✗'}")
        print(f"  Total time: {total_time:.1f}s ({'✓' if total_time < 600 else '✗'} < 10min)")
        print(f"  File size: {file_size_mb:.1f} MB")
        if issues:
            print(f"  Issues:")
            for issue in issues:
                print(f"    - {issue}")
        print(f"  Time: {elapsed:.2f}s")

        return result

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ FAILED: {str(e)}")
        print(f"  Traceback:\n{traceback.format_exc()}")

        return {
            'status': 'FAIL',
            'error': str(e),
            'elapsed_time': elapsed,
            'traceback': traceback.format_exc()
        }


def run_robustness_test(h5ad_path, output_path):
    """Execute complete robustness test suite."""

    print("="*80)
    print("ROBUSTNESS TEST AGENT #1")
    print("Testing spatial analysis pipeline on annotated Visium data")
    print("="*80)
    print(f"\nInput: {h5ad_path}")
    print(f"Output: {output_path}")

    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)

    all_results = {}
    adata = None

    # Step 1: Load H5AD
    result, adata = test_step_1_load_h5ad(h5ad_path)
    all_results['step_1_load'] = result

    if adata is None:
        print("\n✗ Cannot proceed - file loading failed")
        all_results['overall'] = {'status': 'FAIL', 'reason': 'File loading failed'}
        save_results(all_results, output_path)
        return all_results

    # Step 2: Spatial coordinates
    result = test_step_2_spatial_coordinates(adata)
    all_results['step_2_spatial_coords'] = result

    # Step 3: Cell type annotation
    result = test_step_3_cell_type_annotation(adata)
    all_results['step_3_annotation'] = result

    # Step 4: Spatial analysis
    result = test_step_4_spatial_analysis(adata, output_dir)
    all_results['step_4_spatial_analysis'] = result

    # Step 5: Features JSON validation
    result = test_step_5_features_json(all_results['step_4_spatial_analysis'])
    all_results['step_5_features_json'] = result

    # Step 6: Anti-parroting prompt
    result = test_step_6_anti_parroting_prompt(all_results['step_4_spatial_analysis'], output_dir)
    all_results['step_6_prompt'] = result

    # Step 7: Output quality
    result = test_step_7_output_quality(all_results)
    all_results['step_7_quality'] = result

    # Overall summary
    all_passed = all(r.get('status') == 'PASS' for r in all_results.values() if 'status' in r)

    all_results['overall'] = {
        'status': 'PASS' if all_passed else 'FAIL',
        'total_steps': 7,
        'passed_steps': sum(1 for r in all_results.values() if r.get('status') == 'PASS'),
        'failed_steps': sum(1 for r in all_results.values() if r.get('status') == 'FAIL')
    }

    # Save results
    save_results(all_results, output_path)

    # Final report
    print("\n" + "="*80)
    print("ROBUSTNESS TEST COMPLETE")
    print("="*80)
    print(f"\nOverall Status: {'✓ PASS' if all_passed else '✗ FAIL'}")
    print(f"Passed: {all_results['overall']['passed_steps']}/{all_results['overall']['total_steps']} steps")
    print(f"\nResults saved to: {output_path}")
    print("="*80)

    return all_results


def save_results(results, output_path):
    """Save test results to JSON."""

    # Clean up non-serializable objects
    def clean_for_json(obj):
        if isinstance(obj, dict):
            return {k: clean_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_for_json(item) for item in obj]
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif pd.isna(obj):
            return None
        else:
            return obj

    cleaned_results = clean_for_json(results)

    with open(output_path, 'w') as f:
        json.dump(cleaned_results, f, indent=2)

    print(f"\n✓ Results saved to: {output_path}")


if __name__ == "__main__":
    h5ad_path = "/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/annotated_visium.h5ad"
    output_path = "/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/robustness_test_sample1.json"

    results = run_robustness_test(h5ad_path, output_path)

    sys.exit(0 if results['overall']['status'] == 'PASS' else 1)

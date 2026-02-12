"""
End-to-end pipeline validation (no Streamlit required).
Tests the full pipeline from H5AD → annotation → spatial metrics → report prompt → JSON output.
Run as: python scripts/validate_pipeline.py
"""
import sys
import json
import time
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import scanpy as sc
from src.streamlit_adapter import annotate_spatial_regions, calculate_spatial_heterogeneity

H5AD = Path("outputs/annotated_visium.h5ad")

if not H5AD.exists():
    print(f"ERROR: {H5AD} not found. Run notebooks/01_scanpy_baseline.ipynb first.")
    sys.exit(1)


def safe_json_check(obj, path="root"):
    """Recursively check for non-JSON-serializable values."""
    issues = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            issues.extend(safe_json_check(v, f"{path}.{k}"))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            issues.extend(safe_json_check(v, f"{path}[{i}]"))
    elif isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            issues.append(f"{path}: invalid float value {obj}")
    elif isinstance(obj, np.integer):
        pass  # numpy ints are not JSON-serializable but will be caught below
    elif not isinstance(obj, (int, float, str, bool, type(None))):
        issues.append(f"{path}: non-JSON type {type(obj).__name__}")
    return issues


print("=" * 70)
print("END-TO-END PIPELINE VALIDATION")
print(f"Input: {H5AD}")
print("=" * 70)

results = {}
all_pass = True

# Test 1: Load data
print("\n[1/5] Loading H5AD...")
t0 = time.time()
try:
    adata = sc.read_h5ad(H5AD)
    print(f"  PASS: {adata.n_obs:,} spots x {adata.n_vars:,} genes ({time.time()-t0:.1f}s)")
    results['load'] = {'status': 'PASS', 'n_obs': adata.n_obs, 'n_vars': adata.n_vars}
except Exception as e:
    print(f"  FAIL: {e}")
    results['load'] = {'status': 'FAIL', 'error': str(e)}
    all_pass = False
    sys.exit(1)

# Test 2: Annotation
print("\n[2/5] annotate_spatial_regions()...")
t0 = time.time()
try:
    adata_ann, metrics = annotate_spatial_regions(adata, resolution=0.5, use_markers=True, tissue="Unknown")
    elapsed = time.time() - t0
    assert 'n_clusters' in metrics, "Missing n_clusters"
    assert 'n_cell_types' in metrics, "Missing n_cell_types"
    assert 'cell_type_counts' in metrics, "Missing cell_type_counts"
    assert 'spatial_region' in adata_ann.obs.columns, "Missing spatial_region column"
    assert 'cell_type' in adata_ann.obs.columns, "Missing cell_type column"
    print(f"  PASS: {metrics['n_clusters']} clusters, {metrics['n_cell_types']} cell types ({elapsed:.1f}s)")
    results['annotation'] = {'status': 'PASS', **metrics}
except Exception as e:
    print(f"  FAIL: {e}")
    import traceback; traceback.print_exc()
    results['annotation'] = {'status': 'FAIL', 'error': str(e)}
    all_pass = False

# Test 3: Spatial heterogeneity
print("\n[3/5] calculate_spatial_heterogeneity()...")
t0 = time.time()
try:
    spatial_metrics = calculate_spatial_heterogeneity(adata_ann)
    elapsed = time.time() - t0
    assert 'morans_i' in spatial_metrics, "Missing morans_i"
    assert 'spatial_entropy' in spatial_metrics, "Missing spatial_entropy"
    mi = spatial_metrics['morans_i']['mean']
    se = spatial_metrics['spatial_entropy']['mean']
    print(f"  PASS: Moran's I={mi:.3f}, entropy={se:.3f} ({elapsed:.1f}s)")
    results['spatial'] = {'status': 'PASS', **{k: v for k, v in spatial_metrics.items()}}
except Exception as e:
    print(f"  FAIL: {e}")
    results['spatial'] = {'status': 'FAIL', 'error': str(e)}
    all_pass = False

# Test 4: MedGemma prompt generation (demo mode)
print("\n[4/5] MedGemma prompt generation (demo mode)...")
t0 = time.time()
try:
    sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
    from medgemma_report_generator import create_anti_parroting_prompt

    features = {
        'annotation': {
            'n_spots': metrics.get('n_spots', 0),
            'n_clusters': metrics.get('n_clusters', 0),
            'n_cell_types': metrics.get('n_cell_types', 0),
            'mean_confidence': metrics.get('mean_confidence', 0.8),
            'cell_type_counts': metrics.get('cell_type_counts', {}),
        },
        'spatial': {
            'morans_i': spatial_metrics.get('morans_i', {}).get('mean', 0),
            'spatial_entropy': spatial_metrics.get('spatial_entropy', {}).get('mean', 0),
            'n_enriched_pairs': spatial_metrics.get('neighborhood_enrichment', {}).get('n_enriched_pairs', 0),
        }
    }
    prompt = create_anti_parroting_prompt(features)
    assert len(prompt) > 100, "Prompt too short"
    print(f"  PASS: Prompt generated ({len(prompt)} chars, {time.time()-t0:.1f}s)")
    results['prompt'] = {'status': 'PASS', 'prompt_length': len(prompt)}
except Exception as e:
    print(f"  FAIL: {e}")
    import traceback; traceback.print_exc()
    results['prompt'] = {'status': 'FAIL', 'error': str(e)}
    all_pass = False

# Test 5: JSON serialization (no NaN, no numpy types)
print("\n[5/5] JSON serialization check...")
try:
    output = {
        'annotation': results.get('annotation', {}),
        'spatial': results.get('spatial', {}),
    }
    # Convert numpy types for JSON check
    def convert(obj):
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return obj

    output_str = json.dumps(output, default=convert)
    parsed = json.loads(output_str)
    issues = safe_json_check(parsed)
    if issues:
        print(f"  FAIL: {len(issues)} serialization issues:")
        for issue in issues[:5]:
            print(f"    - {issue}")
        results['serialization'] = {'status': 'FAIL', 'issues': issues}
        all_pass = False
    else:
        print(f"  PASS: {len(output_str)} bytes, no NaN/inf/invalid types")
        results['serialization'] = {'status': 'PASS', 'output_bytes': len(output_str)}
except Exception as e:
    print(f"  FAIL: {e}")
    results['serialization'] = {'status': 'FAIL', 'error': str(e)}
    all_pass = False

print()
print("=" * 70)
if all_pass:
    print("RESULT: ALL CHECKS PASSED")
    print("Pipeline is validated end-to-end without Streamlit.")
else:
    failed = [k for k, v in results.items() if v.get('status') == 'FAIL']
    print(f"RESULT: FAILED checks: {failed}")
print("=" * 70)

sys.exit(0 if all_pass else 1)

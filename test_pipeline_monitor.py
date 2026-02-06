#!/usr/bin/env python3
"""
Pipeline monitoring script - runs spatial analysis with real-time progress.
"""

import sys
import time
from datetime import datetime

print("="*80)
print("PIPELINE EXECUTION MONITOR")
print("="*80)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

# Track execution time
start_time = time.time()

# Import and run the spatial analysis
sys.path.insert(0, 'src/spatial_analysis')
from uncertainty_spatial_analysis import run_uncertainty_aware_spatial_analysis

input_file = "outputs/annotated_visium.h5ad"
output_file = "outputs/pipeline_test_20260205/features.json"

print(f"Input: {input_file}")
print(f"Output: {output_file}")
print()

try:
    results = run_uncertainty_aware_spatial_analysis(input_file, output_file)

    elapsed = time.time() - start_time

    print()
    print("="*80)
    print(f"✓ STEP 1 COMPLETED SUCCESSFULLY")
    print(f"Execution time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Output file: {output_file}")
    print("="*80)

except Exception as e:
    elapsed = time.time() - start_time

    print()
    print("="*80)
    print(f"✗ STEP 1 FAILED")
    print(f"Error: {str(e)}")
    print(f"Execution time before failure: {elapsed:.1f} seconds")
    print("="*80)

    import traceback
    traceback.print_exc()

    sys.exit(1)

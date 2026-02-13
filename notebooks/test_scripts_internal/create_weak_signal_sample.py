#!/usr/bin/env python
"""
Create synthetic weak signal sample for testing stage stopping logic.
"""
import json
import numpy as np
from typing import Dict

def create_weak_signal_features() -> Dict:
    """
    Create synthetic uncertainty features with weak spatial signals.

    Stopping criteria:
    1. No STRONG signals (p < 0.001, |I| > 0.3)
    2. < 5 MODERATE signals (p < 0.05, |I| > 0.1)
    3. Spatial entropy < 0.15 (very homogeneous)
    4. All enrichment p-values > 0.1
    """
    # Stage 0: Annotation quality (normal)
    stage0 = {
        "doublet_metrics": {
            "doublet_rate": 0.02,
            "threshold": 0.3,
            "method": "scrublet"
        },
        "annotation_confidence": {
            "high_confidence_spots": 4500,
            "low_confidence_spots": 200,
            "confidence_threshold": 0.7
        }
    }

    # Stage 1: Weak spatial signals
    n_genes = 30
    genes = [f"GENE{i}" for i in range(n_genes)]

    # All genes have weak Moran's I (< 0.1) and high p-values
    morans_i_values = np.random.uniform(0.01, 0.08, n_genes)
    p_values = np.random.uniform(0.2, 0.8, n_genes)

    # All signals are WEAK or NONE
    signal_strength = ['WEAK'] * n_genes

    stage1 = {
        "morans_i": {
            "genes": genes,
            "morans_i": morans_i_values.tolist(),
            "p_value": p_values.tolist(),
            "ci_lower": (morans_i_values - 0.02).tolist(),
            "ci_upper": (morans_i_values + 0.02).tolist(),
            "signal_strength": signal_strength
        },
        "spatial_entropy": {
            "mean_entropy": 0.12,  # Very low (homogeneous)
            "ci_lower": 0.10,
            "ci_upper": 0.14,
            "n_bootstrap": 1000
        },
        "multi_scale_enrichment": {
            "radii": [50, 100, 200],
            "enrichment_pvalues": {
                "CD8_Tcell": [0.3, 0.4, 0.5],
                "CD4_Tcell": [0.6, 0.7, 0.8],
                "Macrophage": [0.2, 0.3, 0.4]
            },
            "scale_stability": 0.2  # Low stability
        }
    }

    # Stage stopping decision
    stopping = {
        "decision": "STOP_WEAK_SIGNAL",
        "rationale": (
            "WEAK SIGNAL DETECTED - Analysis stopped.\n\n"
            "Reasons:\n"
            f"1. STRONG signals: 0 (threshold: ≥1)\n"
            f"2. MODERATE signals: 0 (threshold: ≥5)\n"
            f"3. Spatial entropy: 0.12 (threshold: ≥0.20)\n"
            f"4. Max enrichment confidence: p=0.2 (threshold: p<0.1)\n\n"
            "CONCLUSION: Insufficient spatial signal quality for mechanistic inference."
        ),
        "conditions_failed": [
            "zero_strong_signals",
            "insufficient_moderate_signals",
            "low_entropy",
            "weak_enrichment"
        ]
    }

    return {
        "stage_0_annotation_quality": stage0,
        "stage_1_spatial_patterns": stage1,
        "stage_stopping": stopping,
        "metadata": {
            "created_by": "create_weak_signal_sample.py",
            "purpose": "Test stopping logic",
            "expected_outcome": "Pipeline should stop and generate stopping report"
        }
    }


def main():
    """Generate and save weak signal features."""
    print("Creating synthetic weak signal sample...")

    features = create_weak_signal_features()

    output_path = "outputs/weak_signal_features.json"
    with open(output_path, 'w') as f:
        json.dump(features, f, indent=2)

    print(f"✅ Weak signal features saved to: {output_path}")
    print()
    print("Summary:")
    print(f"  - STRONG signals: 0")
    print(f"  - MODERATE signals: 0")
    print(f"  - WEAK signals: 30")
    print(f"  - Spatial entropy: 0.12 (95% CI: [0.10, 0.14])")
    print(f"  - Stopping decision: STOP_WEAK_SIGNAL")
    print()
    print("Test command:")
    print(f"  python notebooks/medgemma_v2_pipeline.py \\")
    print(f"    --features {output_path} \\")
    print(f"    --output outputs/stopping_report.json \\")
    print(f"    --tissue breast_cancer")


if __name__ == "__main__":
    main()

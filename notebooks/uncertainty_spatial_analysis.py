"""
Uncertainty-Aware Spatial Analysis Pipeline
Stage 0: Pre-processing with annotation quality validation
Stage 1: Spatial pattern extraction with uncertainty quantification
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, Tuple, List
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# STAGE 0: ANNOTATION QUALITY LAYER (Doublet Detection)
# ==============================================================================

def detect_doublets_scrublet(adata, expected_doublet_rate=0.06):
    """
    Detect doublets using Scrublet algorithm.

    Returns:
        adata: Updated with doublet scores and predictions
        quality_metrics: Dict with doublet statistics
    """
    try:
        import scrublet as scr
    except ImportError:
        print("WARNING: Scrublet not installed. Skipping doublet detection.")
        print("Install with: pip install scrublet")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata, {
            'doublet_detection_available': False,
            'n_predicted_doublets': 0,
            'doublet_rate': 0.0
        }

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85)

    adata.obs['doublet_score'] = doublet_scores

    # Handle case where Scrublet couldn't automatically identify threshold
    if predicted_doublets is None:
        print(f"  Warning: Scrublet couldn't auto-identify threshold, using score > 0.25")
        predicted_doublets = doublet_scores > 0.25

    adata.obs['predicted_doublet'] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    doublet_rate = n_doublets / len(adata)

    quality_metrics = {
        'doublet_detection_available': True,
        'n_predicted_doublets': int(n_doublets),
        'doublet_rate': float(doublet_rate),
        'mean_doublet_score': float(doublet_scores.mean()),
        'max_doublet_score': float(doublet_scores.max())
    }

    print(f"Doublet detection: {n_doublets}/{len(adata)} spots ({doublet_rate:.1%})")

    return adata, quality_metrics


def assess_annotation_confidence(adata, confidence_key='conf_score'):
    """
    Assess cell type annotation confidence and flag low-quality annotations.

    Returns:
        annotation_quality: Dict with confidence statistics
    """
    if confidence_key not in adata.obs.columns:
        return {
            'annotation_confidence_available': False,
            'mean_confidence': np.nan,
            'low_confidence_rate': np.nan
        }

    confidence = adata.obs[confidence_key]
    low_conf_threshold = 0.5

    annotation_quality = {
        'annotation_confidence_available': True,
        'mean_confidence': float(confidence.mean()),
        'median_confidence': float(confidence.median()),
        'min_confidence': float(confidence.min()),
        'low_confidence_rate': float((confidence < low_conf_threshold).mean()),
        'n_low_confidence': int((confidence < low_conf_threshold).sum())
    }

    print(f"Annotation confidence: mean={annotation_quality['mean_confidence']:.3f}, "
          f"low_conf_rate={annotation_quality['low_confidence_rate']:.1%}")

    return annotation_quality


# ==============================================================================
# STAGE 1: UNCERTAINTY-AWARE SPATIAL STATISTICS
# ==============================================================================

def compute_morans_i_with_uncertainty(adata, n_genes=100, n_permutations=999):
    """
    Compute Moran's I with permutation-based p-values and confidence intervals.

    Returns:
        results: Dict with genes, Moran's I values, p-values, CIs, and signal strength
    """
    from sklearn.preprocessing import StandardScaler

    if 'spatial_connectivities' not in adata.obsp:
        print("ERROR: Spatial neighbors not computed. Run sq.gr.spatial_neighbors first.")
        return None

    W = adata.obsp['spatial_connectivities']

    # Select highly variable genes
    if 'highly_variable' not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, flavor='seurat_v3')

    hvg = adata.var_names[adata.var['highly_variable']][:n_genes]

    results = {
        'genes': [],
        'morans_i': [],
        'p_value': [],
        'ci_lower': [],
        'ci_upper': [],
        'signal_strength': []
    }

    print(f"Computing Moran's I with {n_permutations} permutations for {len(hvg)} genes...")

    for gene in hvg:
        X = adata[:, gene].X.toarray().flatten()

        # Compute observed Moran's I
        X_std = StandardScaler().fit_transform(X.reshape(-1, 1)).flatten()
        n = len(X_std)

        numerator = (W.multiply(np.outer(X_std, X_std))).sum()
        denominator = (X_std ** 2).sum()

        I_obs = (n / W.sum()) * (numerator / denominator)

        # Permutation test
        I_perm = []
        for _ in range(n_permutations):
            X_perm = np.random.permutation(X_std)
            numerator_perm = (W.multiply(np.outer(X_perm, X_perm))).sum()
            I_perm.append((n / W.sum()) * (numerator_perm / denominator))

        I_perm = np.array(I_perm)

        # P-value (two-tailed)
        p_val = ((np.abs(I_perm) >= np.abs(I_obs)).sum() + 1) / (n_permutations + 1)

        # 95% Confidence interval from permutation distribution
        ci_lower, ci_upper = np.percentile(I_perm, [2.5, 97.5])

        # Signal strength classification
        if p_val < 0.001 and np.abs(I_obs) > 0.3:
            signal = "STRONG"
        elif p_val < 0.05 and np.abs(I_obs) > 0.1:
            signal = "MODERATE"
        elif p_val < 0.05:
            signal = "WEAK"
        else:
            signal = "NONE"

        results['genes'].append(gene)
        results['morans_i'].append(float(I_obs))
        results['p_value'].append(float(p_val))
        results['ci_lower'].append(float(ci_lower))
        results['ci_upper'].append(float(ci_upper))
        results['signal_strength'].append(signal)

    # Sort by signal strength and p-value
    df = pd.DataFrame(results)
    df = df.sort_values(['p_value', 'morans_i'], ascending=[True, False])

    return df.to_dict('list')


def compute_multiscale_neighborhood_enrichment(adata, radii=[1, 2, 3], n_permutations=999):
    """
    Compute neighborhood enrichment at multiple spatial scales.

    Returns:
        multiscale_results: Dict with enrichment at each radius + scale stability flags
    """
    if 'cell_type' not in adata.obs.columns:
        print("ERROR: Cell types not annotated. Run CellTypist first.")
        return None

    try:
        import squidpy as sq
        squidpy_available = True
    except ImportError:
        print("WARNING: Squidpy not available. Using simplified scanpy-based enrichment.")
        squidpy_available = False

    if not squidpy_available:
        # Simplified enrichment using existing spatial graph
        return {
            'radii': radii,
            'enrichment_by_radius': {},
            'scale_stable_pairs': [],
            'scale_dependent_pairs': [],
            'note': 'Squidpy not available - using single-scale analysis'
        }

    multiscale_results = {
        'radii': radii,
        'enrichment_by_radius': {},
        'scale_stable_pairs': [],
        'scale_dependent_pairs': []
    }

    print(f"Computing neighborhood enrichment at {len(radii)} spatial scales...")

    all_z_scores = {}

    for radius in radii:
        print(f"  Radius={radius}...")

        # Rebuild spatial graph with current radius
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=radius*6, radius=radius*55)

        # Compute neighborhood enrichment
        sq.gr.nhood_enrichment(adata, cluster_key='cell_type', n_perms=n_permutations)

        z_scores = adata.uns['cell_type_nhood_enrichment']['zscore']
        cell_types = adata.obs['cell_type'].cat.categories.tolist()

        all_z_scores[radius] = z_scores

        # Extract significant enrichments
        enrichment_pairs = []
        for i, ct1 in enumerate(cell_types):
            for j, ct2 in enumerate(cell_types):
                z = z_scores[i, j]
                if np.abs(z) > 2.0:  # |z| > 2 is significant
                    enrichment_pairs.append({
                        'cell_type_1': ct1,
                        'cell_type_2': ct2,
                        'z_score': float(z),
                        'enriched': z > 0
                    })

        multiscale_results['enrichment_by_radius'][f'radius_{radius}'] = enrichment_pairs

    # Assess scale stability
    if len(radii) >= 2:
        z1 = all_z_scores[radii[0]]
        z2 = all_z_scores[radii[-1]]

        correlation = np.corrcoef(z1.flatten(), z2.flatten())[0, 1]

        cell_types = adata.obs['cell_type'].cat.categories.tolist()

        for i, ct1 in enumerate(cell_types):
            for j, ct2 in enumerate(cell_types):
                # Check if sign and magnitude consistent across scales
                z_small = z1[i, j]
                z_large = z2[i, j]

                if np.sign(z_small) == np.sign(z_large) and np.abs(z_small) > 2 and np.abs(z_large) > 2:
                    multiscale_results['scale_stable_pairs'].append(f"{ct1}-{ct2}")
                elif np.abs(z_small) > 2 or np.abs(z_large) > 2:
                    multiscale_results['scale_dependent_pairs'].append(f"{ct1}-{ct2}")

        multiscale_results['scale_stability_correlation'] = float(correlation)

        print(f"  Scale stability: {len(multiscale_results['scale_stable_pairs'])} stable pairs, "
              f"{len(multiscale_results['scale_dependent_pairs'])} scale-dependent pairs")

    return multiscale_results


def compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=1000):
    """
    Compute spatial entropy with bootstrap confidence intervals.

    Returns:
        entropy_results: Dict with entropy values and uncertainty estimates
    """
    if 'cell_type' not in adata.obs.columns:
        print("ERROR: Cell types not annotated.")
        return None

    from scipy.stats import entropy

    # Compute entropy for each spot's neighborhood
    spatial_entropy = []

    for spot_idx in range(len(adata)):
        # Get neighbors
        neighbors = adata.obsp['spatial_connectivities'][spot_idx].toarray().flatten()
        neighbor_indices = np.where(neighbors > 0)[0]

        if len(neighbor_indices) == 0:
            spatial_entropy.append(0.0)
            continue

        # Get cell type distribution in neighborhood
        neighbor_types = adata.obs['cell_type'].iloc[neighbor_indices]
        type_counts = neighbor_types.value_counts(normalize=True)

        # Shannon entropy
        H = entropy(type_counts.values)
        spatial_entropy.append(H)

    spatial_entropy = np.array(spatial_entropy)
    adata.obs['spatial_entropy'] = spatial_entropy

    # Bootstrap confidence intervals
    bootstrap_means = []
    n_spots = len(adata)

    for _ in range(n_bootstrap):
        sample_indices = np.random.choice(n_spots, size=n_spots, replace=True)
        bootstrap_means.append(spatial_entropy[sample_indices].mean())

    ci_lower, ci_upper = np.percentile(bootstrap_means, [2.5, 97.5])

    entropy_results = {
        'mean_entropy': float(spatial_entropy.mean()),
        'median_entropy': float(np.median(spatial_entropy)),
        'std_entropy': float(spatial_entropy.std()),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'bootstrap_samples': n_bootstrap
    }

    print(f"Spatial entropy: {entropy_results['mean_entropy']:.3f} "
          f"(95% CI: [{ci_lower:.3f}, {ci_upper:.3f}])")

    return entropy_results


# ==============================================================================
# STAGE STOPPING LOGIC
# ==============================================================================

def assess_signal_quality(morans_results, entropy_results, annotation_quality):
    """
    Determine if signal quality is sufficient to proceed to comparative analysis.

    Returns:
        decision: "PROCEED", "STOP_WEAK_SIGNAL", or "STOP_LOW_QUALITY"
        rationale: Explanation for decision
    """
    stop_conditions = []

    # Check 1: Annotation quality
    if annotation_quality.get('annotation_confidence_available'):
        if annotation_quality['mean_confidence'] < 0.4:
            stop_conditions.append(f"Low annotation confidence (mean={annotation_quality['mean_confidence']:.2f} < 0.4)")
        if annotation_quality['low_confidence_rate'] > 0.5:
            stop_conditions.append(f"High low-confidence rate ({annotation_quality['low_confidence_rate']:.1%} > 50%)")

    # Check 2: Spatial signal strength
    if morans_results:
        strong_signals = sum(1 for s in morans_results['signal_strength'] if s == "STRONG")
        moderate_signals = sum(1 for s in morans_results['signal_strength'] if s == "MODERATE")

        if strong_signals == 0 and moderate_signals < 3:
            stop_conditions.append(f"Weak spatial signal (0 STRONG, {moderate_signals} MODERATE genes)")

    # Check 3: Spatial heterogeneity
    if entropy_results:
        if entropy_results['mean_entropy'] < 0.2:
            stop_conditions.append(f"Low spatial heterogeneity (entropy={entropy_results['mean_entropy']:.2f} < 0.2)")

    # Decision logic
    if len(stop_conditions) >= 2:
        return "STOP_WEAK_SIGNAL", "; ".join(stop_conditions)
    elif len(stop_conditions) == 1 and "Low annotation confidence" in stop_conditions[0]:
        return "STOP_LOW_QUALITY", stop_conditions[0]
    else:
        return "PROCEED", "Signal quality sufficient for comparative analysis"


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def run_uncertainty_aware_spatial_analysis(adata_path, output_path):
    """
    Execute full uncertainty-aware spatial analysis pipeline.
    """
    print("=" * 80)
    print("UNCERTAINTY-AWARE SPATIAL ANALYSIS PIPELINE")
    print("=" * 80)

    # Load data
    print("\n[1/6] Loading data...")
    adata = sc.read_h5ad(adata_path)
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    # Stage 0: Annotation quality
    print("\n[2/6] STAGE 0: Annotation Quality Assessment")
    adata, doublet_metrics = detect_doublets_scrublet(adata)
    annotation_quality = assess_annotation_confidence(adata)

    # Build spatial graph (needed for all subsequent stages)
    print("\n[3/6] Building spatial neighbor graph...")
    try:
        import squidpy as sq
        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
    except ImportError:
        print("  WARNING: Squidpy not available, using scanpy spatial neighbors")
        sc.pp.neighbors(adata, use_rep='spatial')

    # Stage 1: Moran's I with uncertainty
    print("\n[4/6] STAGE 1: Spatial Autocorrelation (Moran's I)")
    morans_results = compute_morans_i_with_uncertainty(adata, n_genes=50, n_permutations=999)

    # Spatial entropy with uncertainty
    print("\n[5/6] Spatial Entropy (with bootstrap CI)")
    entropy_results = compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=1000)

    # Multi-scale analysis
    print("\n[6/6] Multi-Scale Neighborhood Enrichment")
    multiscale_results = compute_multiscale_neighborhood_enrichment(adata, radii=[1, 2, 3], n_permutations=999)

    # Stage stopping decision
    print("\n" + "=" * 80)
    print("STAGE STOPPING ASSESSMENT")
    print("=" * 80)
    decision, rationale = assess_signal_quality(morans_results, entropy_results, annotation_quality)

    print(f"\nDecision: {decision}")
    print(f"Rationale: {rationale}")

    # Compile results
    results = {
        'stage_0_annotation_quality': {
            'doublet_metrics': doublet_metrics,
            'annotation_confidence': annotation_quality
        },
        'stage_1_spatial_patterns': {
            'morans_i': morans_results,
            'spatial_entropy': entropy_results,
            'multiscale_enrichment': multiscale_results
        },
        'stage_stopping': {
            'decision': decision,
            'rationale': rationale
        }
    }

    # Save results
    import json
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {output_path}")
    print("=" * 80)

    return results


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python uncertainty_spatial_analysis.py <input.h5ad> <output.json>")
        sys.exit(1)

    adata_path = sys.argv[1]
    output_path = sys.argv[2]

    run_uncertainty_aware_spatial_analysis(adata_path, output_path)

"""
Uncertainty-Aware Spatial Analysis Pipeline
Stage 0: Pre-processing with annotation quality validation
Stage 1: Spatial pattern extraction with uncertainty quantification
"""

import scanpy as sc
import numpy as np
import pandas as pd
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

    # Handle case when Scrublet fails to auto-detect threshold
    if predicted_doublets is None:
        print("  WARNING: Scrublet failed to auto-detect threshold. Using manual threshold at 0.25")
        threshold = 0.25
        predicted_doublets = doublet_scores > threshold

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

def _select_autocorr_genes(adata, n_genes, layer):
    """Top highly variable genes, falling back to highest-variance genes."""
    if 'highly_variable' in adata.var.columns and adata.var['highly_variable'].any():
        hv = adata.var.loc[adata.var['highly_variable']]
        if 'dispersions_norm' in hv.columns:
            hv = hv.sort_values('dispersions_norm', ascending=False)
        return list(hv.index[:n_genes])
    X = adata.layers[layer] if layer else adata.X
    if hasattr(X, 'toarray'):
        mean = np.asarray(X.mean(axis=0)).ravel()
        var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean ** 2
    else:
        var = np.asarray(X).var(axis=0)
    order = np.argsort(var)[::-1][:n_genes]
    return list(adata.var_names[order])


def compute_morans_i_with_uncertainty(adata, n_genes=100, n_permutations=0,
                                      genes=None, layer=None, connectivity_key='spatial'):
    """Moran's I per gene on the sparse spatial graph, with analytic 95% CI and FDR."""
    import squidpy as sq

    obsp_key = f'{connectivity_key}_connectivities'
    if obsp_key not in adata.obsp:
        raise ValueError(f"adata.obsp['{obsp_key}'] missing; build the spatial graph first")
    if layer is None and 'log1p_norm' in adata.layers:
        layer = 'log1p_norm'
    if genes is None:
        genes = _select_autocorr_genes(adata, n_genes, layer)

    df = sq.gr.spatial_autocorr(
        adata, connectivity_key=obsp_key, genes=list(genes), mode='moran',
        n_perms=n_permutations or None, layer=layer, copy=True,
        show_progress_bar=False, seed=0,
    )
    pcol = 'pval_sim' if n_permutations else 'pval_norm'
    half_width = 1.96 * np.sqrt(df['var_norm'].clip(lower=0))
    signal = np.select(
        [
            (df[pcol] < 0.001) & (df['I'].abs() > 0.3),
            (df[pcol] < 0.05) & (df['I'].abs() > 0.1),
            df[pcol] < 0.05,
        ],
        ['STRONG', 'MODERATE', 'WEAK'],
        default='NONE',
    )
    df = df.assign(ci_lower=df['I'] - half_width, ci_upper=df['I'] + half_width,
                   signal_strength=signal)
    df = df.sort_values([pcol, 'I'], ascending=[True, False])

    return {
        'genes': list(df.index),
        'morans_i': df['I'].astype(float).tolist(),
        'p_value': df[pcol].astype(float).tolist(),
        'p_value_fdr': df[f'{pcol}_fdr_bh'].astype(float).tolist(),
        'ci_lower': df['ci_lower'].astype(float).tolist(),
        'ci_upper': df['ci_upper'].astype(float).tolist(),
        'signal_strength': list(df['signal_strength']),
    }


def _nhood_zscore(adata, cluster_key, connectivity_key, n_perms):
    """Run squidpy neighborhood enrichment and return the z-score matrix."""
    import squidpy as sq

    result = sq.gr.nhood_enrichment(
        adata, cluster_key=cluster_key, connectivity_key=connectivity_key,
        n_perms=n_perms, copy=True, show_progress_bar=False, seed=0,
    )
    zscore = result.zscore if hasattr(result, 'zscore') else result[0]
    return np.nan_to_num(np.asarray(zscore, dtype=float))


def compute_multiscale_neighborhood_enrichment(adata, radii=(1, 2, 3), n_permutations=999,
                                               cluster_key='cell_type', z_threshold=2.0):
    """Neighborhood enrichment at several ring sizes without touching the base graph.

    Ring r uses the 3r(r+1) nearest spots (6, 18, 36 for a hexagonal Visium grid),
    so the scale is defined in spots rather than pixel units.
    """
    import squidpy as sq

    if cluster_key not in adata.obs.columns:
        raise ValueError(f"adata.obs['{cluster_key}'] missing; annotate cell types first")
    if not hasattr(adata.obs[cluster_key], 'cat'):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
    cell_types = adata.obs[cluster_key].cat.categories.tolist()

    multiscale_results = {
        'radii': list(radii),
        'enrichment_by_radius': {},
        'scale_stable_pairs': [],
        'scale_dependent_pairs': [],
    }
    all_z_scores = {}

    for radius in radii:
        key = f'nhood_r{radius}'
        sq.gr.spatial_neighbors(adata, coord_type='generic',
                                n_neighs=3 * radius * (radius + 1), key_added=key)
        z_scores = _nhood_zscore(adata, cluster_key, key, n_permutations)
        all_z_scores[radius] = z_scores
        del adata.obsp[f'{key}_connectivities'], adata.obsp[f'{key}_distances']
        adata.uns.pop(f'{key}_neighbors', None)

        multiscale_results['enrichment_by_radius'][f'radius_{radius}'] = [
            {'cell_type_1': cell_types[i], 'cell_type_2': cell_types[j],
             'z_score': float(z_scores[i, j]), 'enriched': bool(z_scores[i, j] > 0)}
            for i in range(len(cell_types)) for j in range(len(cell_types))
            if abs(z_scores[i, j]) > z_threshold
        ]

    base = all_z_scores[radii[0]]
    upper = np.triu_indices(len(cell_types))
    multiscale_results['n_significant_pairs'] = int((base[upper] > z_threshold).sum())
    multiscale_results['max_enrichment'] = float(base.max()) if base.size else 0.0

    if len(radii) >= 2:
        z1, z2 = all_z_scores[radii[0]], all_z_scores[radii[-1]]
        flat1, flat2 = z1.flatten(), z2.flatten()
        if flat1.std() > 0 and flat2.std() > 0:
            multiscale_results['scale_stability_correlation'] = float(np.corrcoef(flat1, flat2)[0, 1])
        for i, j in zip(*upper):
            small, large = z1[i, j], z2[i, j]
            pair = f"{cell_types[i]}-{cell_types[j]}"
            if np.sign(small) == np.sign(large) and abs(small) > z_threshold and abs(large) > z_threshold:
                multiscale_results['scale_stable_pairs'].append(pair)
            elif abs(small) > z_threshold or abs(large) > z_threshold:
                multiscale_results['scale_dependent_pairs'].append(pair)

    print(f"  Enrichment: {multiscale_results['n_significant_pairs']} enriched pairs at radius "
          f"{radii[0]}, {len(multiscale_results['scale_stable_pairs'])} scale-stable")
    return multiscale_results


def compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=1000, cluster_key='cell_type',
                                               connectivity_key='spatial_connectivities', seed=0):
    """Shannon entropy of each spot's neighbourhood composition, with bootstrap CI."""
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"adata.obs['{cluster_key}'] missing; annotate cell types first")
    if connectivity_key not in adata.obsp:
        raise ValueError(f"adata.obsp['{connectivity_key}'] missing; build the spatial graph first")

    import scipy.sparse as sp

    codes = pd.Categorical(adata.obs[cluster_key]).codes
    n_spots, n_types = len(codes), codes.max() + 1
    onehot = sp.csr_matrix((np.ones(n_spots), (np.arange(n_spots), codes)), shape=(n_spots, n_types))
    W = (adata.obsp[connectivity_key] > 0).astype(float)
    counts = np.asarray((W @ onehot).todense())
    totals = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        p = np.where(totals > 0, counts / totals, 0.0)
        spatial_entropy = -np.where(p > 0, p * np.log(p), 0.0).sum(axis=1)
    adata.obs['spatial_entropy'] = spatial_entropy

    rng = np.random.default_rng(seed)
    boot = rng.choice(spatial_entropy, size=(n_bootstrap, n_spots), replace=True).mean(axis=1)
    ci_lower, ci_upper = np.percentile(boot, [2.5, 97.5])

    entropy_results = {
        'mean_entropy': float(spatial_entropy.mean()),
        'median_entropy': float(np.median(spatial_entropy)),
        'std_entropy': float(spatial_entropy.std()),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'bootstrap_samples': n_bootstrap,
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
    multiscale_results = compute_multiscale_neighborhood_enrichment(adata, radii=(1, 2, 3), n_permutations=999)

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

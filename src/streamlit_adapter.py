"""
Streamlit Adapter for V2 Pipeline

Provides interface functions expected by Streamlit app that call V2 pipeline components.
"""

import sys
from pathlib import Path
import tempfile
import json

import math
import numpy as np
import scanpy as sc
# squidpy imported lazily inside calculate_spatial_heterogeneity to avoid
# cudf/RAPIDS initialization on GPUs with compute capability < 7.0 (e.g. P100)


def _safe_float(value, default=0.0):
    """Return default if value is None, NaN, or infinite."""
    try:
        v = float(value)
        return default if (math.isnan(v) or math.isinf(v)) else v
    except (TypeError, ValueError):
        return default

# Import V2 pipeline functions
from src.spatial_analysis.uncertainty_spatial_analysis import (
    detect_doublets_scrublet,
    compute_morans_i_with_uncertainty,
    compute_multiscale_neighborhood_enrichment,
    compute_spatial_entropy_with_bootstrapping,
    assess_annotation_confidence
)


def annotate_spatial_regions(adata, resolution=0.5, use_markers=True, tissue="Unknown"):
    """
    Annotate spatial regions with cell types and clusters.

    Wrapper for V2 pipeline stages: QC, clustering, and annotation.
    Returns (adata, metrics_dict).
    """
    # Make a copy to avoid modifying original
    adata = adata.copy()

    # Check if data is already preprocessed
    is_preprocessed = (
        'highly_variable' in adata.var.columns or
        'X_pca' in adata.obsm or
        adata.n_vars < 3000  # Likely already filtered
    )

    if not is_preprocessed:
        # Stage 0: QC and preprocessing
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        # Calculate QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Filter based on QC
        adata = adata[adata.obs['pct_counts_mt'] < 20].copy()

        # Normalize and log-transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Feature selection
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)

        # Save log1p matrix before scaling (needed for marker scoring)
        adata.layers['log1p_norm'] = adata.X.copy()

        # Scale data
        sc.pp.scale(adata, max_value=10)

    # PCA (check if already computed)
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution, key_added='spatial_region')

    # Tissue-agnostic annotation:
    #   Stage 1: CellTypist Immune_All_High on immune-enriched spots only
    #            (avoids "Epithelial cells" catch-all on non-immune spots)
    #   Stage 2: z-score marker scoring for all spots across pan-tissue panels
    #   Merge: CellTypist wins on immune-enriched spots; markers fill the rest
    annotated = False

    if use_markers:
        try:
            import celltypist
            from celltypist import models
            import pandas as pd

            # Build log1p counts for marker scoring (not scaled)
            adata_ct = adata.copy()
            if adata_ct.X.min() < 0:
                if 'log1p_norm' in adata_ct.layers:
                    adata_ct.X = adata_ct.layers['log1p_norm'].copy()
                elif adata_ct.raw is not None:
                    raw_adata = adata_ct.raw.to_adata()
                    raw_adata.obs = adata_ct.obs.copy()
                    adata_ct = raw_adata
                else:
                    sc.pp.normalize_total(adata_ct, target_sum=1e4)
                    sc.pp.log1p(adata_ct)

            X = np.array(adata_ct.X.toarray() if hasattr(adata_ct.X, 'toarray') else adata_ct.X)
            var_names = list(adata_ct.var_names)

            # Pan-tissue marker panels (works for any unknown tissue)
            compartment_markers = {
                'Plasma_cells':      ['JCHAIN', 'MZB1', 'SDC1', 'CD38', 'TNFRSF17'],
                'Macrophages':       ['CD68', 'CSF1R', 'MRC1', 'MSR1', 'MARCO'],
                'T_cells':           ['CD3D', 'CD3E', 'TRAC', 'CD2', 'CD7'],
                'B_cells':           ['MS4A1', 'CD19', 'CD79A', 'CD79B', 'PAX5'],
                'NK_cells':          ['NKG7', 'GNLY', 'NCAM1', 'KLRD1', 'FCGR3A'],
                'Luminal_secretory': ['SCGB1A1', 'SCGB2A2', 'LTF', 'PIGR', 'SCGB3A1'],
                'Luminal_HR':        ['ESR1', 'PGR', 'FOXA1', 'GATA3', 'AR'],
                'Myoepithelial':     ['TP63', 'KRT14', 'KRT5', 'ACTA2', 'MYLK'],
                'Stromal_CAF':       ['COL1A1', 'COL1A2', 'FAP', 'PDPN', 'POSTN'],
                'Endothelial':       ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'ENG'],
                'Epithelial':        ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
                'Neuronal':          ['MAP2', 'RBFOX3', 'SYP', 'SNAP25', 'GAD1'],
                'Hepatocyte':        ['ALB', 'APOA1', 'APOB', 'CYP3A4', 'HP'],
            }

            # z-score each gene for specificity, then score each spot per compartment
            gene_std = X.std(axis=0); gene_std[gene_std == 0] = 1
            X_z = (X - X.mean(axis=0)) / gene_std

            marker_label = np.full(adata_ct.n_obs, 'Unknown', dtype=object)
            marker_score = np.full(adata_ct.n_obs, -np.inf)

            for comp, genes in compartment_markers.items():
                present = [g for g in genes if g in var_names]
                if len(present) < 2:
                    continue
                idx = [var_names.index(g) for g in present]
                scores = X_z[:, idx].mean(axis=1)
                mask = scores > marker_score
                marker_label[mask] = comp
                marker_score[mask] = scores[mask]

            # Stage 1: Identify immune-enriched spots by pan-immune marker signal
            pan_immune_genes = [g for g in ['PTPRC', 'CD3D', 'CD19', 'MS4A1', 'CD14',
                                             'CSF1R', 'NKG7', 'JCHAIN'] if g in var_names]
            if pan_immune_genes:
                idx = [var_names.index(g) for g in pan_immune_genes]
                immune_signal = X_z[:, idx].mean(axis=1)
                immune_mask = immune_signal >= np.percentile(immune_signal, 70)
            else:
                immune_mask = np.zeros(adata_ct.n_obs, dtype=bool)

            # Run CellTypist only on immune-enriched spots
            final_labels = marker_label.copy()
            if immune_mask.sum() > 10:
                try:
                    print(f"CellTypist: annotating {immune_mask.sum()} immune spots...")
                    model_imm = models.Model.load(model='Immune_All_High.pkl')
                    adata_imm = adata_ct[immune_mask].copy()
                    pred = celltypist.annotate(adata_imm, model=model_imm, majority_voting=False)
                    ct_labels = pred.predicted_labels['predicted_labels'].values.astype(str)
                    ct_probs = pred.probability_matrix.max(axis=1).values

                    # Only accept CellTypist labels that are NOT its non-immune catch-alls
                    NON_IMMUNE = {'Epithelial cells', 'Endothelial cells', 'Unassigned'}
                    for i, (spot_idx) in enumerate(np.where(immune_mask)[0]):
                        if ct_labels[i] not in NON_IMMUNE and ct_probs[i] > 0.4:
                            final_labels[spot_idx] = ct_labels[i]
                    print(f"CellTypist assigned {(~np.isin(final_labels[immune_mask], list(NON_IMMUNE))).sum()} immune spots")
                except Exception as e:
                    print(f"CellTypist immune pass failed: {e}")

            adata.obs['cell_type'] = pd.Categorical(final_labels)
            n_types = len(adata.obs['cell_type'].unique())
            print(f"Annotation complete: {n_types} cell types")
            annotated = True

        except Exception as e:
            print(f"Annotation failed: {e}. Using cluster labels.")

    if not annotated:
        adata.obs['cell_type'] = adata.obs['spatial_region'].astype('category')

    # Compute confidence scores (returns dict, not adata)
    confidence_metrics = assess_annotation_confidence(adata, confidence_key='leiden_confidence')

    # Extract metrics
    metrics = {
        'n_spots': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_clusters': int(len(adata.obs['spatial_region'].unique())),
        'n_cell_types': int(len(adata.obs['cell_type'].unique())),
        'mean_confidence': _safe_float(confidence_metrics.get('mean_confidence'), 0.8),
        'cluster_distribution': adata.obs['spatial_region'].value_counts().to_dict(),
        # Used by report prompt function
        'cell_type_counts': adata.obs['cell_type'].value_counts().to_dict()
    }

    return adata, metrics


def calculate_spatial_heterogeneity(adata):
    """
    Calculate spatial heterogeneity metrics.

    Wrapper for V2 pipeline spatial analysis functions.
    Returns metrics_dict.
    """
    import squidpy as sq  # lazy import — avoids cudf init on low-CC GPUs
    metrics = {}

    # Ensure spatial graph exists
    if 'spatial_neighbors' not in adata.uns:
        if 'spatial' in adata.obsm:
            sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
        else:
            # Fallback: use PCA neighbors
            sc.pp.neighbors(adata, n_neighbors=6)

    # Stage 1: Moran's I (spatial autocorrelation)
    try:
        morans_results = compute_morans_i_with_uncertainty(adata, n_genes=50, n_permutations=100)

        if morans_results is not None:
            # Results is dict of lists - compute summary stats
            morans_i_values = morans_results.get('morans_i', [])
            p_values = morans_results.get('p_value', [])

            # Count significant genes (p < 0.05)
            n_significant = sum(1 for p in p_values if p < 0.05) if p_values else 0

            metrics['morans_i'] = {
                'mean': float(np.mean(morans_i_values)) if morans_i_values else 0.0,
                'significant_genes': int(n_significant),
                'p_value_mean': float(np.mean(p_values)) if p_values else 1.0
            }
        else:
            metrics['morans_i'] = {'mean': 0.0, 'significant_genes': 0, 'p_value_mean': 1.0}

    except Exception as e:
        print(f"Moran's I failed: {e}")
        metrics['morans_i'] = {'mean': 0.0, 'significant_genes': 0, 'p_value_mean': 1.0}

    # Stage 2: Neighborhood enrichment (spatial clustering)
    try:
        enrich_results = compute_multiscale_neighborhood_enrichment(
            adata,
            radii=[1, 2, 3],
            n_permutations=100
        )

        # Extract summary metrics
        metrics['neighborhood_enrichment'] = {
            'n_enriched_pairs': int(enrich_results.get('n_significant_pairs', 0)),
            'max_enrichment': float(enrich_results.get('max_enrichment', 1.0))
        }
    except Exception as e:
        print(f"Neighborhood enrichment failed: {e}")
        metrics['neighborhood_enrichment'] = {'n_enriched_pairs': 0, 'max_enrichment': 1.0}

    # Stage 3: Spatial entropy (diversity)
    try:
        entropy_results = compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=100)
        metrics['spatial_entropy'] = {
            'mean': float(entropy_results['mean_entropy']),
            'std': float(entropy_results['std_entropy']),
            'ci_lower': float(entropy_results['ci_lower']),
            'ci_upper': float(entropy_results['ci_upper'])
        }
    except Exception as e:
        print(f"Spatial entropy failed: {e}")
        metrics['spatial_entropy'] = {'mean': 0.0, 'std': 0.0, 'ci_lower': 0.0, 'ci_upper': 0.0}

    return metrics

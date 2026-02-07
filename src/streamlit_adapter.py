"""
Streamlit Adapter for V2 Pipeline

Provides interface functions expected by Streamlit app that call V2 pipeline components.
"""

import sys
from pathlib import Path
import tempfile
import json

import numpy as np
import scanpy as sc
import squidpy as sq

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

        # Scale data
        sc.pp.scale(adata, max_value=10)

    # PCA (check if already computed)
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution, key_added='spatial_region')

    # Marker-based annotation (if requested)
    if use_markers:
        try:
            # Load PanglaoDB markers
            markers_path = Path(__file__).parent.parent / 'data' / 'PanglaoDB_markers_27_Mar_2020.tsv'

            if markers_path.exists():
                import pandas as pd
                markers_df = pd.read_csv(markers_path, sep='\t')

                # Create marker dict for annotation
                marker_dict = {}
                for _, row in markers_df.iterrows():
                    cell_type = row['cell type']
                    gene = row['official gene symbol']

                    if cell_type not in marker_dict:
                        marker_dict[cell_type] = []
                    marker_dict[cell_type].append(gene)

                # Score marker genes for each cluster
                cluster_annotations = {}

                for cluster in adata.obs['spatial_region'].unique():
                    cluster_cells = adata[adata.obs['spatial_region'] == cluster]

                    # Get mean expression for this cluster
                    mean_expr = cluster_cells.X.mean(axis=0)
                    if hasattr(mean_expr, 'A1'):
                        mean_expr = mean_expr.A1

                    # Score each cell type
                    best_score = 0
                    best_type = f"Cluster_{cluster}"

                    for cell_type, genes in marker_dict.items():
                        # Find marker genes present in dataset
                        present_markers = [g for g in genes if g in adata.var_names]

                        if present_markers:
                            marker_indices = [adata.var_names.get_loc(g) for g in present_markers]
                            score = mean_expr[marker_indices].mean()

                            if score > best_score:
                                best_score = score
                                best_type = cell_type

                    cluster_annotations[cluster] = best_type

                # Map clusters to cell types (as categorical)
                adata.obs['cell_type'] = adata.obs['spatial_region'].map(cluster_annotations).astype('category')

            else:
                # Fallback: use cluster names
                adata.obs['cell_type'] = adata.obs['spatial_region'].astype('category')

        except Exception as e:
            print(f"Marker annotation failed: {e}. Using cluster labels.")
            adata.obs['cell_type'] = adata.obs['spatial_region'].astype('category')
    else:
        # Just use cluster labels as cell types
        adata.obs['cell_type'] = adata.obs['spatial_region'].astype('category')

    # Compute confidence scores (returns dict, not adata)
    confidence_metrics = assess_annotation_confidence(adata, confidence_key='leiden_confidence')

    # Extract metrics
    metrics = {
        'n_spots': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_clusters': int(len(adata.obs['spatial_region'].unique())),
        'n_cell_types': int(len(adata.obs['cell_type'].unique())),
        'mean_confidence': confidence_metrics.get('mean_confidence', 0.8),
        'cluster_distribution': adata.obs['spatial_region'].value_counts().to_dict()
    }

    return adata, metrics


def calculate_spatial_heterogeneity(adata):
    """
    Calculate spatial heterogeneity metrics.

    Wrapper for V2 pipeline spatial analysis functions.
    Returns metrics_dict.
    """
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

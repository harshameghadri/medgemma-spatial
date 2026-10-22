"""
Pipeline entry points used by the Streamlit app.

Produces the feature contract consumed by the report prompts and the LoRA
training data (see ``build_features``):

    {
      "annotation": {n_spots, n_genes, n_clusters, n_cell_types, mean_confidence,
                     annotation_method, cell_type_counts, cluster_distribution},
      "spatial_heterogeneity": {morans_i_mean, morans_i_n_significant,
                                morans_i_top_genes, n_enriched_pairs, max_enrichment,
                                entropy_mean, entropy_ci_lower, entropy_ci_upper,
                                gene_signatures},
      "uncertainty": {mean_prediction_entropy},
      "errors": [str, ...]
    }

A metric that could not be computed is ``None`` and the reason is appended to
``errors``; it is never replaced with a plausible-looking default.
"""

import math

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from src.spatial_analysis.uncertainty_spatial_analysis import (
    compute_morans_i_with_uncertainty,
    compute_multiscale_neighborhood_enrichment,
    compute_spatial_entropy_with_bootstrapping,
)

COMPARTMENT_MARKERS = {
    'Plasma_cells':      ['JCHAIN', 'MZB1', 'SDC1', 'CD38', 'TNFRSF17'],
    'Macrophages':       ['CD68', 'CSF1R', 'MRC1', 'MSR1', 'MARCO'],
    'T_cells':           ['CD3D', 'CD3E', 'TRAC', 'CD2', 'CD7'],
    'B_cells':           ['MS4A1', 'CD19', 'CD79A', 'CD79B', 'PAX5'],
    'NK_cells':          ['NKG7', 'GNLY', 'NCAM1', 'KLRD1', 'FCGR3A'],
    'Dendritic_cells':   ['CLEC9A', 'FCER1A', 'CD1C', 'ITGAX', 'HLA-DRA'],
    'Mast_cells':        ['KIT', 'TPSAB1', 'CPA3', 'MS4A2', 'FCER1G'],
    'Epithelial':        ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'Luminal_secretory': ['SCGB1A1', 'SCGB2A2', 'LTF', 'PIGR', 'SCGB3A1'],
    'Luminal_HR':        ['ESR1', 'PGR', 'FOXA1', 'GATA3', 'AR'],
    'Myoepithelial':     ['TP63', 'KRT14', 'KRT5', 'ACTA2', 'MYLK'],
    'Goblet_cells':      ['MUC2', 'TFF3', 'FCGBP', 'SPDEF', 'AGR2'],
    'Enteroendocrine':   ['CHGA', 'CHGB', 'NEUROD1', 'PAX4', 'ISL1'],
    'Tuft_cells':        ['POU2F3', 'DCLK1', 'TRPM5', 'CHAT', 'GFI1B'],
    'CAF_generic':       ['COL1A1', 'COL1A2', 'FAP', 'PDPN', 'POSTN'],
    'CAF_myCAF':         ['ACTA2', 'MYH11', 'TAGLN', 'TPM2', 'PDGFRB'],
    'CAF_iCAF':          ['IL6', 'CXCL12', 'PDGFRA', 'HAS1', 'CFD'],
    'Endothelial':       ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'ENG'],
    'Neuronal':          ['MAP2', 'RBFOX3', 'SYP', 'SNAP25', 'GAD1'],
    'Hepatocyte':        ['ALB', 'APOA1', 'APOB', 'CYP3A4', 'HP'],
}
PAN_IMMUNE_GENES = ['PTPRC', 'CD3D', 'CD19', 'MS4A1', 'CD14', 'CSF1R', 'NKG7', 'JCHAIN']
CELLTYPIST_NON_IMMUNE = {'Epithelial_cells', 'Endothelial_cells', 'Unassigned'}


def _safe_float(value, default=None):
    """Return default if value is None, NaN, or infinite."""
    try:
        v = float(value)
        return default if (math.isnan(v) or math.isinf(v)) else v
    except (TypeError, ValueError):
        return default


def looks_like_raw_counts(X, n_sample: int = 20000) -> bool:
    """True when the matrix holds non-negative integers (unnormalised UMI counts)."""
    values = X.data if sp.issparse(X) else np.asarray(X).ravel()
    values = values[:n_sample]
    if values.size == 0:
        return False
    return bool(np.all(values >= 0) and np.allclose(values, np.round(values)))


def _has_negative(X) -> bool:
    return bool((X.min() if sp.issparse(X) else np.min(X)) < 0)


def restore_full_matrix(adata):
    """Swap a scaled HVG-only X for the full matrix in adata.raw, keeping spot metadata."""
    full = adata.raw.to_adata()
    full.obs = adata.obs.copy()
    for key in adata.obsm:
        full.obsm[key] = adata.obsm[key]
    for key in adata.obsp:
        full.obsp[key] = adata.obsp[key]
    full.uns = adata.uns
    if 'highly_variable' in adata.var.columns:
        full.var['highly_variable'] = full.var_names.isin(adata.var_names[adata.var['highly_variable']])
    print(f"Using adata.raw ({full.n_vars} genes) instead of scaled X ({adata.n_vars} genes)")
    return full


def preprocess(adata, errors):
    """Normalise raw counts if needed and compute HVGs + PCA on a scaled HVG copy."""
    if adata.raw is not None and _has_negative(adata.X) and adata.raw.n_vars >= adata.n_vars:
        adata = restore_full_matrix(adata)

    if looks_like_raw_counts(adata.X):
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        if adata.var['mt'].any():
            keep = adata.obs['pct_counts_mt'] < 20
            if keep.sum() < 0.5 * adata.n_obs:
                errors.append(f"QC: {int((~keep).sum())}/{adata.n_obs} spots exceed 20% mito; mito filter skipped")
            else:
                adata._inplace_subset_obs(keep.values)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers['log1p_norm'] = adata.X.copy()
    elif 'log1p_norm' not in adata.layers:
        if not _has_negative(adata.X):
            adata.layers['log1p_norm'] = adata.X.copy()
        else:
            errors.append("Input X is scaled (negative values) and has no log1p layer; marker scoring may be unreliable")

    if 'highly_variable' not in adata.var.columns:
        layer = 'log1p_norm' if 'log1p_norm' in adata.layers else None
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars), flavor='seurat', layer=layer)

    if 'X_pca' not in adata.obsm:
        hv = adata[:, adata.var['highly_variable'].values].copy()
        if 'log1p_norm' in hv.layers:
            hv.X = hv.layers['log1p_norm']
        sc.pp.scale(hv, max_value=10)
        n_comps = max(2, min(50, hv.n_vars - 1, hv.n_obs - 1))
        sc.tl.pca(hv, n_comps=n_comps, svd_solver='arpack')
        adata.obsm['X_pca'] = hv.obsm['X_pca']
    return adata


def score_compartments(adata):
    """Mean z-scored marker expression per compartment; returns (labels, confidence, entropy)."""
    matrix = adata.layers['log1p_norm'] if 'log1p_norm' in adata.layers else adata.X
    var_index = {g: i for i, g in enumerate(adata.var_names)}
    panels = {
        comp: [var_index[g] for g in genes if g in var_index]
        for comp, genes in COMPARTMENT_MARKERS.items()
    }
    panels = {comp: idx for comp, idx in panels.items() if len(idx) >= 2}
    if not panels:
        raise ValueError("No marker panel has >=2 genes present in var_names")

    genes = sorted({i for idx in panels.values() for i in idx})
    col = {g: k for k, g in enumerate(genes)}
    X = matrix[:, genes]
    X = np.asarray(X.toarray() if sp.issparse(X) else X, dtype=float)
    std = X.std(axis=0)
    std[std == 0] = 1.0
    Z = (X - X.mean(axis=0)) / std

    names = list(panels)
    scores = np.column_stack([Z[:, [col[i] for i in panels[c]]].mean(axis=1) for c in names])
    labels = np.array(names, dtype=object)[scores.argmax(axis=1)]

    shifted = scores - scores.max(axis=1, keepdims=True)
    prob = np.exp(shifted)
    prob /= prob.sum(axis=1, keepdims=True)
    confidence = prob.max(axis=1)
    entropy = -(prob * np.log(np.clip(prob, 1e-12, None))).sum(axis=1)
    return labels, confidence, entropy


def refine_immune_with_celltypist(adata, labels, confidence, errors):
    """Relabel immune-enriched spots with CellTypist Immune_All_High where confident."""
    import celltypist
    from celltypist import models

    present = [g for g in PAN_IMMUNE_GENES if g in adata.var_names]
    if not present or 'log1p_norm' not in adata.layers:
        return labels, confidence
    signal = adata[:, present].layers['log1p_norm']
    signal = np.asarray(signal.mean(axis=1)).ravel()
    immune_mask = signal >= np.percentile(signal, 70)
    if immune_mask.sum() <= 10:
        return labels, confidence

    import anndata as ad

    adata_imm = ad.AnnData(adata.layers['log1p_norm'][immune_mask], var=adata.var[[]].copy())
    adata_imm.obs_names = adata.obs_names[immune_mask]
    try:
        model = models.Model.load(model='Immune_All_High.pkl')
        pred = celltypist.annotate(adata_imm, model=model, majority_voting=False)
    except Exception as e:
        errors.append(f"CellTypist immune refinement skipped: {e}")
        return labels, confidence

    ct_labels = np.char.replace(pred.predicted_labels['predicted_labels'].values.astype(str), ' ', '_')
    ct_probs = pred.probability_matrix.max(axis=1).values
    accept = np.array([lab not in CELLTYPIST_NON_IMMUNE for lab in ct_labels]) & (ct_probs > 0.4)
    idx = np.where(immune_mask)[0][accept]
    labels = labels.copy()
    confidence = confidence.copy()
    labels[idx] = ct_labels[accept]
    confidence[idx] = ct_probs[accept]
    print(f"CellTypist relabelled {accept.sum()} of {immune_mask.sum()} immune-enriched spots")
    return labels, confidence


def annotate_spatial_regions(adata, resolution=0.5, use_markers=True, tissue="Unknown",
                             use_celltypist=True):
    """QC, cluster and annotate spots. Returns (adata, annotation_metrics)."""
    adata = adata.copy()
    errors = []
    adata = preprocess(adata, errors)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40, adata.obsm['X_pca'].shape[1]))
    sc.tl.leiden(adata, resolution=resolution, key_added='spatial_region',
                 flavor='igraph', n_iterations=2, directed=False)

    method = 'leiden_clusters'
    if use_markers:
        try:
            labels, confidence, entropy = score_compartments(adata)
            method = 'marker_zscore'
            if use_celltypist:
                try:
                    n_errors = len(errors)
                    labels, confidence = refine_immune_with_celltypist(adata, labels, confidence, errors)
                    if len(errors) == n_errors:
                        method = 'marker_zscore+celltypist'
                except ImportError as e:
                    errors.append(f"CellTypist unavailable: {e}")
            adata.obs['cell_type'] = pd.Categorical(labels)
            adata.obs['annotation_confidence'] = confidence
            adata.obs['prediction_entropy'] = entropy
        except Exception as e:
            errors.append(f"Marker annotation failed, using Leiden clusters: {e}")

    if 'cell_type' not in adata.obs.columns:
        adata.obs['cell_type'] = adata.obs['spatial_region'].astype('category')

    mean_confidence = None
    if 'annotation_confidence' in adata.obs.columns:
        mean_confidence = _safe_float(adata.obs['annotation_confidence'].mean())

    metrics = {
        'n_spots': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'n_clusters': int(adata.obs['spatial_region'].nunique()),
        'n_cell_types': int(adata.obs['cell_type'].nunique()),
        'mean_confidence': mean_confidence,
        'annotation_method': method,
        'cluster_distribution': {str(k): int(v) for k, v in adata.obs['spatial_region'].value_counts().items()},
        'cell_type_counts': {str(k): int(v) for k, v in adata.obs['cell_type'].value_counts().items() if v > 0},
        'errors': errors,
    }
    print(f"Annotation complete ({method}): {metrics['n_cell_types']} cell types")
    return adata, metrics


def ensure_spatial_graph(adata, n_neighs=6):
    """Build the base 6-NN spatial graph after validating coordinates."""
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] missing; spatial statistics need spot coordinates")
    coords = np.asarray(adata.obsm['spatial'], dtype=float)
    if not np.isfinite(coords).all():
        raise ValueError(f"{int((~np.isfinite(coords)).any(axis=1).sum())} spots have NaN/inf coordinates")
    if 'spatial_connectivities' not in adata.obsp:
        import squidpy as sq

        sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=n_neighs)


def calculate_spatial_heterogeneity(adata, n_genes=50, n_permutations=100):
    """Moran's I, neighbourhood enrichment, neighbourhood entropy and gene signatures."""
    errors = []
    metrics = {
        'morans_i_mean': None,
        'morans_i_n_significant': None,
        'morans_i_top_genes': [],
        'n_enriched_pairs': None,
        'max_enrichment': None,
        'entropy_mean': None,
        'entropy_ci_lower': None,
        'entropy_ci_upper': None,
        'gene_signatures': {},
        'errors': errors,
    }

    try:
        ensure_spatial_graph(adata)
    except Exception as e:
        errors.append(f"Spatial graph: {e}")
        return metrics

    try:
        morans = compute_morans_i_with_uncertainty(adata, n_genes=n_genes)
        metrics['morans_i_mean'] = _safe_float(np.mean(morans['morans_i']))
        metrics['morans_i_n_significant'] = int(sum(p < 0.05 for p in morans['p_value_fdr']))
        metrics['morans_i_top_genes'] = morans['genes'][:5]
    except Exception as e:
        errors.append(f"Moran's I: {e}")

    try:
        enrich = compute_multiscale_neighborhood_enrichment(
            adata, radii=(1, 2, 3), n_permutations=n_permutations)
        metrics['n_enriched_pairs'] = enrich['n_significant_pairs']
        metrics['max_enrichment'] = _safe_float(enrich['max_enrichment'])
    except Exception as e:
        errors.append(f"Neighborhood enrichment: {e}")

    try:
        entropy = compute_spatial_entropy_with_bootstrapping(adata, n_bootstrap=200)
        metrics['entropy_mean'] = entropy['mean_entropy']
        metrics['entropy_ci_lower'] = entropy['ci_lower']
        metrics['entropy_ci_upper'] = entropy['ci_upper']
    except Exception as e:
        errors.append(f"Spatial entropy: {e}")

    try:
        from src.spatial_analysis.gene_signature_scoring import run_full_signature_analysis

        sig_results, _ = run_full_signature_analysis(
            adata, cell_type_col='cell_type' if 'cell_type' in adata.obs.columns else 'spatial_region')
        metrics['gene_signatures'] = sig_results
    except Exception as e:
        errors.append(f"Gene signatures: {e}")

    for msg in errors:
        print(f"WARNING: {msg}")
    return metrics


def build_features(adata, annotation_metrics, spatial_metrics):
    """Assemble the feature contract passed to prompts and stored with reports."""
    annotation = {k: v for k, v in annotation_metrics.items() if k != 'errors'}
    spatial = {k: v for k, v in spatial_metrics.items() if k != 'errors'}
    pred_entropy = None
    if 'prediction_entropy' in adata.obs.columns:
        pred_entropy = _safe_float(adata.obs['prediction_entropy'].mean())
    return {
        'annotation': annotation,
        'spatial_heterogeneity': spatial,
        'uncertainty': {'mean_prediction_entropy': pred_entropy},
        'errors': list(annotation_metrics.get('errors', [])) + list(spatial_metrics.get('errors', [])),
    }

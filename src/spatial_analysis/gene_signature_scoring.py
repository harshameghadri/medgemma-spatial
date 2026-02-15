"""
Gene signature scoring and spatial boundary detection for Visium spatial transcriptomics.

Enables spot-level gene expression queries:
  - Score each spot for biologically meaningful gene signatures
  - Detect tumor-normal and compartment boundary spots
  - Find top genes enriched at spatial interfaces (e.g. NOTUM+ epithelial, IGHG1+ plasma)
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')


SPATIAL_SIGNATURES: Dict[str, List[str]] = {
    # Epithelial boundary / tumor-normal interface
    'tumor_boundary':    ['NOTUM', 'CLDN3', 'CLDN4', 'EPCAM', 'KRT20', 'LGR5', 'AXIN2'],
    # IGHG1*3+ plasma cells near tumor margins (key immunotherapy biomarker)
    'plasma_ighg':       ['IGHG1', 'IGHG3', 'IGHG4', 'MZB1', 'JCHAIN', 'IGKC'],
    # T cell exhaustion
    'exhaustion':        ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX'],
    # Activated / myofibroblastic CAF
    'caf_activation':    ['FAP', 'POSTN', 'ACTA2', 'PDPN', 'S100A4', 'LRRC15'],
    # Hypoxic core
    'hypoxia':           ['HIF1A', 'VEGFA', 'SLC2A1', 'LDHA', 'CA9', 'BNIP3'],
    # Proliferating spots
    'proliferation':     ['MKI67', 'TOP2A', 'CCNB1', 'CDK1', 'PCNA', 'STMN1'],
    # Epithelial-to-mesenchymal transition
    'emt':               ['VIM', 'CDH2', 'FN1', 'SNAI1', 'ZEB1', 'TWIST1'],
    # Tertiary lymphoid structures (TLS)
    'tls':               ['CXCL13', 'BCL6', 'CD21', 'LAMP3', 'CCL19', 'SELL'],
    # M2 macrophage / immunosuppressive myeloid
    'm2_macrophage':     ['CD163', 'MRC1', 'IL10', 'TGFB1', 'ARG1', 'CCL18'],
    # Complement / innate immune activation
    'complement':        ['C1QA', 'C1QB', 'C1QC', 'C3', 'CFB', 'CFD'],
    # Goblet / secretory maturation (colon-specific)
    'goblet_maturation': ['MUC2', 'TFF3', 'FCGBP', 'AGR2', 'SPDEF', 'CLCA1'],
    # Wnt pathway activity (colon cancer driver)
    'wnt_activity':      ['AXIN2', 'LGR5', 'ASCL2', 'NOTUM', 'EPHB2', 'CD44'],
}


def score_gene_signatures(
    adata: sc.AnnData,
    signatures: Optional[Dict[str, List[str]]] = None,
    ctrl_size: int = 50,
    use_raw: bool = False,
) -> sc.AnnData:
    """Score each spot for predefined biological gene signatures.

    Adds columns `sig_{name}` to adata.obs for each signature.
    Uses scanpy score_genes (mean expression - background control set).
    """
    if signatures is None:
        signatures = SPATIAL_SIGNATURES

    adata_score = adata.copy()

    # Ensure log-normalised counts for scoring
    if adata_score.X.min() < 0:
        if 'log1p_norm' in adata_score.layers:
            adata_score.X = adata_score.layers['log1p_norm'].copy()
        elif adata_score.raw is not None:
            tmp = adata_score.raw.to_adata()
            tmp.obs = adata_score.obs.copy()
            adata_score = tmp
        else:
            sc.pp.normalize_total(adata_score, target_sum=1e4)
            sc.pp.log1p(adata_score)

    scored_sigs = {}
    for sig_name, gene_list in signatures.items():
        present = [g for g in gene_list if g in adata_score.var_names]
        if len(present) < 2:
            adata.obs[f'sig_{sig_name}'] = 0.0
            scored_sigs[sig_name] = {'n_genes_available': len(present), 'skipped': True}
            continue

        ctrl = min(ctrl_size, len(present) * 5)
        try:
            sc.tl.score_genes(adata_score, gene_list=present, ctrl_size=ctrl,
                              score_name=f'sig_{sig_name}')
            adata.obs[f'sig_{sig_name}'] = adata_score.obs[f'sig_{sig_name}'].values
            scored_sigs[sig_name] = {
                'n_genes_available': len(present),
                'genes_used': present,
                'mean_score': float(adata.obs[f'sig_{sig_name}'].mean()),
                'skipped': False,
            }
        except Exception as e:
            adata.obs[f'sig_{sig_name}'] = 0.0
            scored_sigs[sig_name] = {'error': str(e), 'skipped': True}

    adata.uns['signature_scoring'] = scored_sigs
    return adata


def detect_spatial_boundaries(
    adata: sc.AnnData,
    cell_type_col: str = 'cell_type',
    n_neighbors: int = 6,
    boundary_percentile: float = 80.0,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """Identify spots at cell-type boundaries (tumor-normal / compartment interfaces).

    A boundary spot has a high fraction of spatially adjacent neighbors with
    a DIFFERENT cell type label. Returns a boolean mask and summary DataFrame.
    """
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Column '{cell_type_col}' not found in adata.obs")

    # Build spatial graph if not present
    import squidpy as sq
    if 'spatial_connectivities' not in adata.obsp:
        if 'spatial' in adata.obsm:
            sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=n_neighbors)
        else:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)

    conn = adata.obsp.get('spatial_connectivities') or adata.obsp.get('connectivities')
    labels = adata.obs[cell_type_col].values

    n_spots = adata.n_obs
    heterogeneity_score = np.zeros(n_spots)
    dominant_neighbor_type = np.full(n_spots, 'Unknown', dtype=object)

    for i in range(n_spots):
        row = conn[i]
        neighbor_indices = row.indices if hasattr(row, 'indices') else np.where(np.asarray(row).flatten() > 0)[0]
        if len(neighbor_indices) == 0:
            continue
        neighbor_labels = labels[neighbor_indices]
        n_different = np.sum(neighbor_labels != labels[i])
        heterogeneity_score[i] = n_different / len(neighbor_indices)

        # Most common neighbor cell type (excluding self type)
        other = neighbor_labels[neighbor_labels != labels[i]]
        if len(other) > 0:
            vals, counts = np.unique(other, return_counts=True)
            dominant_neighbor_type[i] = vals[counts.argmax()]

    threshold = np.percentile(heterogeneity_score, boundary_percentile)
    boundary_mask = heterogeneity_score >= threshold

    adata.obs['boundary_score'] = heterogeneity_score
    adata.obs['is_boundary'] = boundary_mask
    adata.obs['dominant_neighbor_type'] = dominant_neighbor_type

    boundary_df = pd.DataFrame({
        'barcode': adata.obs_names[boundary_mask],
        'cell_type': labels[boundary_mask],
        'boundary_score': heterogeneity_score[boundary_mask],
        'dominant_neighbor_type': dominant_neighbor_type[boundary_mask],
    }).sort_values('boundary_score', ascending=False).reset_index(drop=True)

    return boundary_mask, boundary_df


def score_boundary_genes(
    adata: sc.AnnData,
    boundary_mask: np.ndarray,
    cell_type_filter: Optional[str] = None,
    n_top: int = 20,
    method: str = 'wilcoxon',
) -> pd.DataFrame:
    """Differential gene expression: boundary spots vs core spots.

    Finds genes enriched at spatial interfaces — e.g. NOTUM+ epithelial near margins,
    IGHG1+ plasma cells near tumor-normal boundaries.
    """
    if boundary_mask.sum() < 5:
        return pd.DataFrame(columns=['gene', 'score', 'pval', 'logfc', 'group'])

    adata_test = adata.copy()

    # Optionally restrict to one cell type (e.g. epithelial only)
    if cell_type_filter and 'cell_type' in adata_test.obs.columns:
        ct_mask = adata_test.obs['cell_type'].str.lower().str.contains(
            cell_type_filter.lower(), na=False
        )
        if ct_mask.sum() < 10:
            cell_type_filter = None  # fall back to all spots
        else:
            adata_test = adata_test[ct_mask].copy()
            boundary_mask = adata_test.obs.get('is_boundary', pd.Series(False)).values

    if boundary_mask.sum() < 5 or (~boundary_mask).sum() < 5:
        return pd.DataFrame(columns=['gene', 'score', 'pval', 'logfc', 'group'])

    # Ensure log-normalised expression
    if adata_test.X.min() < 0:
        if 'log1p_norm' in adata_test.layers:
            adata_test.X = adata_test.layers['log1p_norm'].copy()
        else:
            sc.pp.normalize_total(adata_test, target_sum=1e4)
            sc.pp.log1p(adata_test)

    adata_test.obs['zone'] = np.where(boundary_mask, 'boundary', 'core')

    try:
        sc.tl.rank_genes_groups(
            adata_test, groupby='zone', groups=['boundary'],
            reference='core', method=method, n_genes=n_top
        )
        result = sc.get.rank_genes_groups_df(adata_test, group='boundary')
        result = result.rename(columns={'names': 'gene', 'pvals_adj': 'pval',
                                         'logfoldchanges': 'logfc'})
        result['group'] = cell_type_filter or 'all'
        return result[['gene', 'scores', 'pval', 'logfc', 'group']].head(n_top)

    except Exception as e:
        return pd.DataFrame({'gene': [str(e)], 'scores': [0], 'pval': [1],
                             'logfc': [0], 'group': ['error']})


def get_signature_spatial_stats(
    adata: sc.AnnData,
    sig_name: str,
    cell_type_col: str = 'cell_type',
    top_pct: float = 10.0,
) -> dict:
    """Spatial statistics for a scored gene signature.

    Returns mean score, top-N% threshold, Moran's I (spatial clustering),
    and co-occurrence with cell types (which cell types score highest).
    """
    col = f'sig_{sig_name}'
    if col not in adata.obs.columns:
        return {'error': f'{col} not found — run score_gene_signatures first'}

    scores = adata.obs[col].values
    threshold = np.percentile(scores, 100 - top_pct)
    high_mask = scores >= threshold

    stats: dict = {
        'sig_name': sig_name,
        'mean_score': float(scores.mean()),
        'std_score': float(scores.std()),
        'top_pct_threshold': float(threshold),
        'n_high_spots': int(high_mask.sum()),
    }

    # Cell type enrichment in top-scoring spots
    if cell_type_col in adata.obs.columns and high_mask.sum() > 0:
        ct_counts = adata.obs.loc[high_mask, cell_type_col].value_counts()
        ct_pct = ct_counts / high_mask.sum() * 100
        stats['top_cell_types'] = ct_pct.head(5).to_dict()

    # Moran's I for the signature score (is it spatially organised?)
    try:
        import squidpy as sq
        if 'spatial_connectivities' not in adata.obsp and 'spatial' in adata.obsm:
            sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

        adata.obs[f'_morans_tmp_{sig_name}'] = scores
        sq.gr.spatial_autocorr(adata, mode='moran',
                               genes=[f'_morans_tmp_{sig_name}'],
                               layer=None, use_raw=False)
        mi_key = 'moranI' if 'moranI' in adata.uns else 'spatial_autocorr'
        if mi_key in adata.uns:
            mi_df = adata.uns[mi_key]
            gene_key = f'_morans_tmp_{sig_name}'
            if gene_key in mi_df.index:
                stats['morans_i'] = float(mi_df.loc[gene_key, 'I'])
                stats['morans_pval'] = float(mi_df.loc[gene_key, 'pval_norm'])
        del adata.obs[f'_morans_tmp_{sig_name}']
    except Exception:
        pass

    return stats


def run_full_signature_analysis(
    adata: sc.AnnData,
    signatures: Optional[Dict[str, List[str]]] = None,
    cell_type_col: str = 'cell_type',
) -> dict:
    """Convenience wrapper: score signatures + detect boundaries + get stats.

    Returns a dict ready to merge into the pipeline features JSON under
    the 'gene_signatures' key.
    """
    results: dict = {}

    # 1. Score all signatures
    adata = score_gene_signatures(adata, signatures=signatures)

    # 2. Spatial stats per signature
    sigs_to_use = signatures or SPATIAL_SIGNATURES
    sig_stats = {}
    for sig_name in sigs_to_use:
        col = f'sig_{sig_name}'
        if col in adata.obs.columns and adata.obs[col].abs().sum() > 0:
            sig_stats[sig_name] = get_signature_spatial_stats(adata, sig_name, cell_type_col)
    results['signature_stats'] = sig_stats

    # 3. Boundary detection (requires cell_type column)
    if cell_type_col in adata.obs.columns:
        try:
            boundary_mask, boundary_df = detect_spatial_boundaries(
                adata, cell_type_col=cell_type_col
            )
            results['n_boundary_spots'] = int(boundary_mask.sum())
            results['boundary_top_transitions'] = (
                boundary_df.groupby(['cell_type', 'dominant_neighbor_type'])
                .size()
                .reset_index(name='n')
                .sort_values('n', ascending=False)
                .head(10)
                .to_dict(orient='records')
            )

            # 4. Top boundary genes (all cell types combined)
            boundary_genes_df = score_boundary_genes(adata, boundary_mask, n_top=20)
            if not boundary_genes_df.empty:
                results['top_boundary_genes'] = boundary_genes_df[['gene', 'logfc', 'pval']].head(10).to_dict(orient='records')

            # 5. Epithelial-specific boundary genes (finds NOTUM-like signals)
            epi_mask = adata.obs.get('is_boundary', pd.Series(False)).values & (
                adata.obs.get(cell_type_col, pd.Series('')).str.contains('pithelial|Goblet|Tuft|Entero', na=False).values
            )
            if epi_mask.sum() >= 5:
                epi_genes = score_boundary_genes(adata, epi_mask, cell_type_filter='Epithelial', n_top=10)
                if not epi_genes.empty:
                    results['epithelial_boundary_genes'] = epi_genes[['gene', 'logfc']].to_dict(orient='records')

        except Exception as e:
            results['boundary_error'] = str(e)

    return results, adata

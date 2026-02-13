"""
Marker gene-based cell type annotation using PanglaoDB.
"""

import pandas as pd
import numpy as np
import scanpy as sc
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')


def load_panglaodb_markers(marker_file: str = "data/PanglaoDB_markers_27_Mar_2020.tsv",
                           species: str = "Hs",
                           min_sensitivity: float = 0.5,
                           min_specificity: float = 0.5) -> pd.DataFrame:
    """Load and filter PanglaoDB marker database."""
    print(f"Loading PanglaoDB markers from: {marker_file}")
    markers = pd.read_csv(marker_file, sep='\t')

    markers = markers[markers['species'].str.contains(species, na=False)]
    markers = markers[markers['canonical marker'] == 1]

    sens_col = f'sensitivity_{species.lower()}'
    spec_col = f'specificity_{species.lower()}'

    if sens_col in markers.columns and spec_col in markers.columns:
        markers[sens_col] = markers[sens_col].fillna(0.0)
        markers[spec_col] = markers[spec_col].fillna(0.0)
        markers = markers[
            (markers[sens_col] >= min_sensitivity) |
            (markers[spec_col] >= min_specificity)
        ]

    print(f"  Loaded {len(markers):,} markers for {markers['cell type'].nunique()} cell types")
    return markers


def get_tissue_specific_markers(markers_df: pd.DataFrame,
                                tissue: str = None) -> Dict[str, List[str]]:
    """Extract cell type → marker gene dictionary, optionally filtered by tissue."""
    if tissue:
        markers_df = markers_df[
            markers_df['organ'].str.contains(tissue, case=False, na=False) |
            markers_df['organ'].isna()
        ]

    marker_dict = {}
    for cell_type, group in markers_df.groupby('cell type'):
        marker_dict[cell_type] = group['official gene symbol'].unique().tolist()

    return marker_dict


def score_clusters_for_markers(adata,
                                marker_dict: Dict[str, List[str]],
                                cluster_key: str = 'spatial_region',
                                use_raw: bool = False) -> pd.DataFrame:
    """Score each cluster for each cell type using marker gene expression."""
    print(f"\nScoring {adata.obs[cluster_key].nunique()} clusters for {len(marker_dict)} cell types...")

    expr_matrix = adata.raw.to_adata() if (use_raw and adata.raw is not None) else adata

    cluster_means = pd.DataFrame(
        index=expr_matrix.var_names,
        columns=expr_matrix.obs[cluster_key].unique()
    )

    for cluster in expr_matrix.obs[cluster_key].unique():
        cluster_cells = expr_matrix.obs[cluster_key] == cluster
        cluster_means[cluster] = np.array(expr_matrix[cluster_cells].X.mean(axis=0)).flatten()

    scores = pd.DataFrame(
        index=cluster_means.columns,
        columns=marker_dict.keys(),
        dtype=float
    )

    for cell_type, genes in marker_dict.items():
        available_genes = [g for g in genes if g in cluster_means.index]
        scores[cell_type] = cluster_means.loc[available_genes].mean(axis=0) if available_genes else 0.0

    scores = (scores - scores.mean()) / (scores.std() + 1e-10)
    return scores


def assign_cell_type_labels(adata,
                             score_matrix: pd.DataFrame,
                             cluster_key: str = 'spatial_region',
                             min_score: float = 0.5,
                             top_n_markers: int = 3) -> pd.DataFrame:
    """Assign biological labels to clusters based on marker scores."""
    print(f"\nAssigning cell type labels (min_score={min_score})...")

    annotations = []

    for cluster in score_matrix.index:
        cluster_scores = score_matrix.loc[cluster].sort_values(ascending=False)
        top_types = cluster_scores.head(top_n_markers)

        if top_types.iloc[0] >= min_score:
            primary_type = top_types.index[0]
            primary_score = top_types.iloc[0]
        else:
            primary_type = "Unknown"
            primary_score = 0.0

        cluster_size = (adata.obs[cluster_key] == cluster).sum()
        annotations.append({
            'cluster': cluster,
            'cell_type': primary_type,
            'score': primary_score,
            'top_3_types': ', '.join([f"{k}({v:.2f})" for k, v in top_types.items()]),
            'n_spots': cluster_size
        })

    annotation_df = pd.DataFrame(annotations)

    print(f"\n  Annotated {len(annotation_df)} clusters:")
    for _, row in annotation_df.iterrows():
        print(f"    Cluster {row['cluster']}: {row['cell_type']} "
              f"(score={row['score']:.2f}, n={row['n_spots']:,})")

    return annotation_df


def annotate_clusters_with_markers(adata,
                                    marker_file: str = "data/PanglaoDB_markers_27_Mar_2020.tsv",
                                    tissue: str = None,
                                    cluster_key: str = 'spatial_region',
                                    species: str = "Hs") -> Tuple:
    """Complete workflow: load markers, score clusters, assign labels."""
    print("=" * 80)
    print("MARKER-BASED CELL TYPE ANNOTATION")
    print("=" * 80)

    markers = load_panglaodb_markers(marker_file, species=species)
    marker_dict_full = get_tissue_specific_markers(markers, tissue=None)
    print(f"\n  Total cell types in database: {len(marker_dict_full)}")

    priority_types = {
        'T cells', 'B cells', 'Macrophages', 'Fibroblasts',
        'Epithelial cells', 'Endothelial cells', 'Smooth muscle cells',
        'Goblet cells', 'Enterocytes', 'Plasma cells', 'Dendritic cells',
        'NK cells', 'Neutrophils', 'Mast cells', 'Pericytes'
    }

    marker_dict = {k: v for k, v in marker_dict_full.items() if k in priority_types}

    if len(marker_dict) == 0:
        print("  ⚠ No priority types found, using canonical markers")
        marker_dict = CANONICAL_MARKERS

    print(f"  Using {len(marker_dict)} cell type marker sets")

    score_matrix = score_clusters_for_markers(adata, marker_dict, cluster_key=cluster_key)
    annotation_df = assign_cell_type_labels(adata, score_matrix, cluster_key=cluster_key)

    cluster_to_celltype = dict(zip(annotation_df['cluster'], annotation_df['cell_type']))
    adata.obs['cell_type_predicted'] = adata.obs[cluster_key].map(cluster_to_celltype)

    n_annotated = (annotation_df['cell_type'] != 'Unknown').sum()
    n_total = len(annotation_df)
    pct_annotated = 100 * n_annotated / n_total

    metrics = {
        'n_clusters': n_total,
        'n_annotated': int(n_annotated),
        'pct_annotated': float(pct_annotated),
        'cell_type_counts': adata.obs['cell_type_predicted'].value_counts().to_dict()
    }

    print(f"\n{'=' * 80}")
    print(f"✅ Annotation complete: {n_annotated}/{n_total} clusters ({pct_annotated:.1f}%)")
    print(f"{'=' * 80}")

    return adata, annotation_df, score_matrix, metrics


def run_marker_gene_analysis(adata,
                              cluster_key: str = 'spatial_region',
                              n_genes: int = 10,
                              method: str = 'wilcoxon') -> Dict:
    """Identify differentially expressed marker genes per cluster."""
    print(f"\nIdentifying marker genes per cluster (method={method})...")

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, n_genes=n_genes)

    top_genes = {}
    for cluster in adata.obs[cluster_key].unique():
        genes = sc.get.rank_genes_groups_df(adata, group=cluster).head(n_genes)
        top_genes[cluster] = genes['names'].tolist()
        print(f"  Cluster {cluster} top genes: {', '.join(genes['names'].head(5).tolist())}")

    return top_genes


# Canonical marker gene sets for major cell types (fallback when PanglaoDB unavailable)
CANONICAL_MARKERS = {
    'T cells': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4'],
    'B cells': ['CD19', 'CD79A', 'CD79B', 'MS4A1'],
    'Plasma cells': ['IGHG1', 'IGHG3', 'MZB1', 'SDC1'],
    'Macrophages': ['CD68', 'CD163', 'C1QA', 'C1QB'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC9A'],
    'NK cells': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
    'Neutrophils': ['S100A8', 'S100A9', 'CSF3R'],
    'Mast cells': ['KIT', 'TPSAB1', 'CPA3'],
    'Fibroblasts': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN'],
    'Endothelial cells': ['PECAM1', 'VWF', 'CDH5'],
    'Epithelial cells': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
    'Goblet cells': ['MUC2', 'TFF3', 'FCGBP'],
    'Enterocytes': ['FABP1', 'APOA1', 'APOA4'],
    'Smooth muscle cells': ['ACTA2', 'MYH11', 'TAGLN'],
    'Pericytes': ['RGS5', 'PDGFRB', 'CSPG4'],
}

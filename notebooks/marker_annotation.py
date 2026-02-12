#!/usr/bin/env python3
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
    """
    Load and filter PanglaoDB marker database.

    Args:
        marker_file: Path to PanglaoDB TSV file
        species: "Hs" (human) or "Mm" (mouse)
        min_sensitivity: Minimum sensitivity score
        min_specificity: Minimum specificity score

    Returns:
        Filtered marker dataframe
    """
    print(f"Loading PanglaoDB markers from: {marker_file}")
    markers = pd.read_csv(marker_file, sep='\t')

    # Filter for species
    markers = markers[markers['species'].str.contains(species, na=False)]

    # Filter for canonical markers (high confidence)
    markers = markers[markers['canonical marker'] == 1]

    # Filter by sensitivity and specificity
    sens_col = f'sensitivity_{species.lower()}'
    spec_col = f'specificity_{species.lower()}'

    if sens_col in markers.columns and spec_col in markers.columns:
        # Handle NaN values in sensitivity/specificity columns
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
    """
    Extract cell type marker gene dictionary.

    Args:
        markers_df: PanglaoDB markers dataframe
        tissue: Tissue type filter (e.g., "Colon", "Intestine")

    Returns:
        Dictionary mapping cell types to gene lists
    """
    if tissue:
        markers_df = markers_df[
            markers_df['organ'].str.contains(tissue, case=False, na=False) |
            markers_df['organ'].isna()
        ]

    marker_dict = {}
    for cell_type, group in markers_df.groupby('cell type'):
        genes = group['official gene symbol'].unique().tolist()
        marker_dict[cell_type] = genes

    return marker_dict


def score_clusters_for_markers(adata,
                               marker_dict: Dict[str, List[str]],
                               cluster_key: str = 'spatial_region',
                               use_raw: bool = False) -> pd.DataFrame:
    """
    Score each cluster for each cell type using marker genes.

    Args:
        adata: AnnData object with clustering results
        marker_dict: Cell type → gene list mapping
        cluster_key: Column in adata.obs with cluster labels
        use_raw: Use adata.raw for expression

    Returns:
        DataFrame with scores (clusters × cell types)
    """
    print(f"\nScoring {adata.obs[cluster_key].nunique()} clusters for {len(marker_dict)} cell types...")

    # Get expression matrix
    if use_raw and adata.raw is not None:
        expr_matrix = adata.raw.to_adata()
    else:
        expr_matrix = adata

    # Calculate mean expression per cluster
    cluster_means = pd.DataFrame(
        index=expr_matrix.var_names,
        columns=expr_matrix.obs[cluster_key].unique()
    )

    for cluster in expr_matrix.obs[cluster_key].unique():
        cluster_cells = expr_matrix.obs[cluster_key] == cluster
        cluster_means[cluster] = np.array(expr_matrix[cluster_cells].X.mean(axis=0)).flatten()

    # Score each cluster for each cell type
    scores = pd.DataFrame(
        index=cluster_means.columns,
        columns=marker_dict.keys(),
        dtype=float
    )

    for cell_type, genes in marker_dict.items():
        # Find genes that exist in dataset
        available_genes = [g for g in genes if g in cluster_means.index]

        if len(available_genes) > 0:
            # Mean expression of marker genes
            scores[cell_type] = cluster_means.loc[available_genes].mean(axis=0)
        else:
            scores[cell_type] = 0.0

    # Normalize scores (z-score per cell type)
    scores = (scores - scores.mean()) / (scores.std() + 1e-10)

    return scores


def assign_cell_type_labels(adata,
                            score_matrix: pd.DataFrame,
                            cluster_key: str = 'spatial_region',
                            min_score: float = 0.5,
                            top_n_markers: int = 3) -> pd.DataFrame:
    """
    Assign biological labels to clusters based on marker scores.

    Args:
        adata: AnnData object
        score_matrix: Cluster × cell type score matrix
        cluster_key: Column with cluster labels
        min_score: Minimum z-score to consider
        top_n_markers: Number of top cell types to report

    Returns:
        Cluster annotation dataframe
    """
    print(f"\nAssigning cell type labels (min_score={min_score})...")

    annotations = []

    for cluster in score_matrix.index:
        # Get top scoring cell types
        cluster_scores = score_matrix.loc[cluster].sort_values(ascending=False)
        top_types = cluster_scores.head(top_n_markers)

        # Check if top score meets threshold
        if top_types.iloc[0] >= min_score:
            primary_type = top_types.index[0]
            primary_score = top_types.iloc[0]
        else:
            primary_type = "Unknown"
            primary_score = 0.0

        # Get cluster size
        cluster_size = (adata.obs[cluster_key] == cluster).sum()

        annotations.append({
            'cluster': cluster,
            'cell_type': primary_type,
            'score': primary_score,
            'top_3_types': ', '.join([f"{k}({v:.2f})" for k, v in top_types.items()]),
            'n_spots': cluster_size
        })

    annotation_df = pd.DataFrame(annotations)

    # Print summary
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
    """
    Complete workflow: Load markers, score clusters, assign labels.

    Args:
        adata: AnnData object with clustering
        marker_file: Path to PanglaoDB markers
        tissue: Tissue filter (e.g., "Colon")
        cluster_key: Cluster column name
        species: "Hs" or "Mm"

    Returns:
        (updated adata, annotation_df, score_matrix)
    """
    print("="*80)
    print("MARKER-BASED CELL TYPE ANNOTATION")
    print("="*80)

    # Load markers
    markers = load_panglaodb_markers(marker_file, species=species)

    # Get tissue-specific markers
    marker_dict_full = get_tissue_specific_markers(markers, tissue=None)  # Get all markers first
    print(f"\n  Total cell types in database: {len(marker_dict_full)}")

    # Use canonical markers for common cell types
    # Don't filter by tissue-specific organ, use universal cell type markers
    priority_types = {
        'T cells', 'B cells', 'Macrophages', 'Fibroblasts',
        'Epithelial cells', 'Endothelial cells', 'Smooth muscle cells',
        'Goblet cells', 'Enterocytes', 'Plasma cells', 'Dendritic cells',
        'NK cells', 'Neutrophils', 'Mast cells', 'Pericytes'
    }

    marker_dict = {k: v for k, v in marker_dict_full.items() if k in priority_types}

    # If no priority types found, use CANONICAL_MARKERS as fallback
    if len(marker_dict) == 0:
        print(f"  ⚠ No priority types found in database, using canonical markers")
        marker_dict = CANONICAL_MARKERS

    print(f"  Using {len(marker_dict)} cell type marker sets")

    # Score clusters
    score_matrix = score_clusters_for_markers(adata, marker_dict, cluster_key=cluster_key)

    # Assign labels
    annotation_df = assign_cell_type_labels(adata, score_matrix, cluster_key=cluster_key)

    # Add to adata.obs
    cluster_to_celltype = dict(zip(annotation_df['cluster'], annotation_df['cell_type']))
    adata.obs['cell_type_predicted'] = adata.obs[cluster_key].map(cluster_to_celltype)

    # Calculate metrics
    n_annotated = (annotation_df['cell_type'] != 'Unknown').sum()
    n_total = len(annotation_df)
    pct_annotated = 100 * n_annotated / n_total

    metrics = {
        'n_clusters': n_total,
        'n_annotated': int(n_annotated),
        'pct_annotated': float(pct_annotated),
        'cell_type_counts': adata.obs['cell_type_predicted'].value_counts().to_dict()
    }

    print(f"\n{'='*80}")
    print(f"✅ Annotation complete: {n_annotated}/{n_total} clusters ({pct_annotated:.1f}%)")
    print(f"{'='*80}")

    return adata, annotation_df, score_matrix, metrics


def run_marker_gene_analysis(adata,
                            cluster_key: str = 'spatial_region',
                            n_genes: int = 10,
                            method: str = 'wilcoxon') -> pd.DataFrame:
    """
    Identify differentially expressed marker genes per cluster.

    Args:
        adata: AnnData object
        cluster_key: Cluster column
        n_genes: Number of top genes to report
        method: Statistical test ('wilcoxon', 't-test')

    Returns:
        Top marker genes per cluster
    """
    print(f"\nIdentifying marker genes per cluster (method={method})...")

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, n_genes=n_genes)

    # Extract top genes per cluster
    top_genes = {}
    for cluster in adata.obs[cluster_key].unique():
        genes = sc.get.rank_genes_groups_df(adata, group=cluster).head(n_genes)
        top_genes[cluster] = genes['names'].tolist()

        print(f"  Cluster {cluster} top genes: {', '.join(genes['names'].head(5).tolist())}")

    return top_genes


# Predefined marker gene sets for major cell types
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
    'Pericytes': ['RGS5', 'PDGFRB', 'CSPG4']
}


if __name__ == "__main__":
    print("Marker annotation module loaded.")
    print(f"Available canonical markers for {len(CANONICAL_MARKERS)} cell types")

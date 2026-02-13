#!/usr/bin/env python3
"""
Enhanced Spatial Statistics Analysis (Scanpy-only version)
Fallback version without squidpy dependency for compatibility
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist, pdist
from scipy.stats import entropy
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors
import json
import warnings
import sys
import os

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, facecolor='white')


def calculate_spatial_entropy(adata, cluster_key='cell_type', n_neighbors=30):
    """Calculate Shannon entropy of cell type diversity in local neighborhoods."""
    coords = adata.obsm['spatial']
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1, algorithm='ball_tree').fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    entropy_values = []
    for neighbors in indices:
        neighbor_types = adata.obs[cluster_key].iloc[neighbors[1:]].values
        unique, counts = np.unique(neighbor_types, return_counts=True)
        probs = counts / counts.sum()
        ent = entropy(probs)
        entropy_values.append(ent)

    return np.array(entropy_values)


def calculate_nn_distances(adata, cluster_key='cell_type'):
    """Calculate nearest neighbor distances within cell types."""
    coords = adata.obsm['spatial']
    cell_types = adata.obs[cluster_key].values

    nn_distances = np.full(len(adata), np.nan)

    for ct in np.unique(cell_types):
        mask = cell_types == ct
        if mask.sum() < 2:
            continue

        ct_coords = coords[mask]
        dist_matrix = cdist(ct_coords, ct_coords)
        np.fill_diagonal(dist_matrix, np.inf)
        min_distances = np.min(dist_matrix, axis=1)

        nn_distances[mask] = min_distances

    return nn_distances


def calculate_cluster_compactness(adata, cluster_key='leiden'):
    """Calculate spatial compactness metrics for clusters."""
    coords = adata.obsm['spatial']
    clusters = adata.obs[cluster_key].values

    compactness = {}

    for cluster in np.unique(clusters):
        mask = clusters == cluster
        cluster_coords = coords[mask]

        if len(cluster_coords) < 2:
            continue

        centroid = cluster_coords.mean(axis=0)
        distances = np.linalg.norm(cluster_coords - centroid, axis=1)

        pairwise_dist = pdist(cluster_coords)

        convex_hull_area = 0
        try:
            if len(cluster_coords) >= 3:
                hull = ConvexHull(cluster_coords)
                convex_hull_area = hull.volume
        except:
            pass

        compactness[str(cluster)] = {
            'n_spots': int(mask.sum()),
            'mean_distance_to_centroid': float(distances.mean()),
            'std_distance_to_centroid': float(distances.std()),
            'max_pairwise_distance': float(pairwise_dist.max()),
            'convex_hull_area': float(convex_hull_area),
            'density': float(mask.sum() / (convex_hull_area + 1)),
            'compactness_score': float(distances.mean() / (distances.std() + 1))
        }

    return compactness


def calculate_morans_i_manual(adata, gene_list=None, n_neighbors=6):
    """Calculate Moran's I manually without squidpy."""
    from sklearn.neighbors import NearestNeighbors

    coords = adata.obsm['spatial']
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1, algorithm='ball_tree').fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    W = np.zeros((len(adata), len(adata)))
    for i, neighbors in enumerate(indices):
        for j in neighbors[1:]:
            W[i, j] = 1

    W_sum = W.sum()

    if gene_list is None:
        if 'highly_variable' in adata.var.columns:
            gene_list = adata.var_names[adata.var['highly_variable']].tolist()[:100]
        else:
            gene_list = adata.var_names[:100].tolist()

    results = []
    for gene in gene_list:
        if gene not in adata.var_names:
            continue

        x = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
        x_mean = x.mean()
        x_centered = x - x_mean

        numerator = 0
        for i in range(len(adata)):
            for j in range(len(adata)):
                numerator += W[i, j] * x_centered[i] * x_centered[j]

        denominator = (x_centered ** 2).sum()

        morans_i = (len(adata) / W_sum) * (numerator / denominator) if denominator > 0 else 0

        results.append({
            'gene': gene,
            'I': morans_i
        })

    return pd.DataFrame(results).set_index('gene')


def calculate_neighborhood_enrichment(adata, cluster_key='cell_type', n_neighbors=6):
    """Calculate neighborhood enrichment without squidpy."""
    from sklearn.neighbors import NearestNeighbors
    from scipy.stats import chi2_contingency

    coords = adata.obsm['spatial']
    cell_types = adata.obs[cluster_key].values
    unique_types = np.unique(cell_types)

    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1, algorithm='ball_tree').fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    enrichment_matrix = np.zeros((len(unique_types), len(unique_types)))

    for i, ct1 in enumerate(unique_types):
        for j, ct2 in enumerate(unique_types):
            ct1_mask = cell_types == ct1
            ct1_indices = np.where(ct1_mask)[0]

            if len(ct1_indices) == 0:
                continue

            neighbor_counts = 0
            total_neighbors = 0

            for idx in ct1_indices:
                neighbors = indices[idx][1:]
                neighbor_types = cell_types[neighbors]
                neighbor_counts += (neighbor_types == ct2).sum()
                total_neighbors += len(neighbors)

            observed = neighbor_counts
            expected = total_neighbors * (cell_types == ct2).sum() / len(cell_types)

            if expected > 0:
                enrichment_matrix[i, j] = (observed - expected) / np.sqrt(expected)

    return enrichment_matrix, unique_types


def generate_clinical_summary(metrics):
    """Generate clinical interpretation of spatial statistics."""
    summary = []

    summary.append("SPATIAL ANALYSIS SUMMARY")
    summary.append("=" * 50)
    summary.append(f"\nSample: {metrics['sample_info']['n_spots']} spots, {metrics['sample_info']['n_cell_types']} cell types")

    summary.append("\n1. SPATIAL HETEROGENEITY:")
    entropy_interp = metrics['spatial_entropy']['overall']['interpretation']
    summary.append(f"   - Overall heterogeneity: {entropy_interp.upper()}")
    summary.append(f"   - Median entropy: {metrics['spatial_entropy']['overall']['median']:.3f}")

    summary.append("\n2. CELL TYPE SPATIAL PATTERNS:")
    for ct, stats in metrics['nearest_neighbor_distances'].items():
        pattern = stats['spatial_pattern']
        summary.append(f"   - {ct}: {pattern.replace('_', ' ').upper()}")

    if 'neighborhood_enrichment' in metrics and metrics['neighborhood_enrichment']:
        summary.append("\n3. NEIGHBORHOOD INTERACTIONS:")
        for ct, enrichment in metrics['neighborhood_enrichment'].items():
            if enrichment['n_enriched'] > 0:
                summary.append(f"   - {ct} enriched near: {', '.join(enrichment['enriched_neighbors'])}")

    if 'spatial_autocorrelation' in metrics and metrics['spatial_autocorrelation']:
        summary.append("\n4. SPATIALLY VARIABLE GENES:")
        n_sig = metrics['spatial_autocorrelation'].get('n_genes_tested', 0)
        top_genes = metrics['spatial_autocorrelation'].get('top_genes', [])[:5]
        summary.append(f"   - {n_sig} genes tested for spatial patterns")
        if top_genes:
            summary.append(f"   - Top spatial genes: {', '.join(top_genes)}")

    summary.append("\n5. CLUSTER ORGANIZATION:")
    for cluster, comp in metrics['cluster_compactness'].items():
        if comp['n_spots'] > 100:
            score = comp['compactness_score']
            org = "tightly organized" if score > 1.5 else "dispersed" if score < 0.8 else "moderate"
            summary.append(f"   - Cluster {cluster}: {org} ({comp['n_spots']} spots)")

    return "\n".join(summary)


def main(input_h5ad=None, output_dir=None):
    """Run comprehensive spatial statistics analysis."""

    if input_h5ad is None:
        input_h5ad = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs/annotated_visium_enhanced.h5ad'

    if output_dir is None:
        output_dir = '/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma/outputs'

    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading data from: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    print(f"Loaded: {adata.shape[0]} spots, {adata.shape[1]} genes")

    cell_types = adata.obs['cell_type'].unique()
    print(f"Cell types: {cell_types.tolist()}")

    print("\n1. Calculating neighborhood enrichment...")
    try:
        enrichment_matrix, unique_types = calculate_neighborhood_enrichment(adata, cluster_key='cell_type', n_neighbors=6)

        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(enrichment_matrix, xticklabels=unique_types, yticklabels=unique_types,
                   center=0, cmap='RdBu_r', ax=ax, cbar_kws={'label': 'Z-score'})
        ax.set_title('Neighborhood Enrichment')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/neighborhood_enrichment.png', dpi=300, bbox_inches='tight')
        plt.close()

        enrichment_summary = {}
        for i, ct1 in enumerate(unique_types):
            enriched_neighbors = []
            depleted_neighbors = []

            for j, ct2 in enumerate(unique_types):
                if i != j:
                    zscore = enrichment_matrix[i, j]

                    if zscore > 2:
                        enriched_neighbors.append(ct2)
                    elif zscore < -2:
                        depleted_neighbors.append(ct2)

            enrichment_summary[ct1] = {
                'enriched_neighbors': enriched_neighbors,
                'depleted_neighbors': depleted_neighbors,
                'n_enriched': len(enriched_neighbors),
                'n_depleted': len(depleted_neighbors)
            }
        print("   ✓ Neighborhood enrichment complete")
    except Exception as e:
        print(f"   ✗ Neighborhood enrichment failed: {e}")
        enrichment_summary = {}

    print("\n2. Calculating spatial entropy...")
    try:
        spatial_entropy = calculate_spatial_entropy(adata, cluster_key='cell_type', n_neighbors=30)
        adata.obs['spatial_entropy'] = spatial_entropy

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        sc.pl.spatial(adata, color='spatial_entropy', ax=axes[0], show=False, title='Spatial Entropy', cmap='viridis')
        axes[1].hist(spatial_entropy, bins=50, edgecolor='black')
        axes[1].set_xlabel('Spatial Entropy')
        axes[1].set_ylabel('Count')
        axes[1].axvline(np.median(spatial_entropy), color='red', linestyle='--', label=f'Median: {np.median(spatial_entropy):.2f}')
        axes[1].legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/spatial_entropy.png', dpi=300, bbox_inches='tight')
        plt.close()

        entropy_stats = {
            'mean': float(np.mean(spatial_entropy)),
            'median': float(np.median(spatial_entropy)),
            'std': float(np.std(spatial_entropy)),
            'min': float(np.min(spatial_entropy)),
            'max': float(np.max(spatial_entropy)),
            'interpretation': 'high' if np.median(spatial_entropy) > 1.0 else 'moderate' if np.median(spatial_entropy) > 0.5 else 'low'
        }

        entropy_by_type = {}
        for ct in cell_types:
            mask = adata.obs['cell_type'] == ct
            entropy_by_type[ct] = {
                'mean_entropy': float(np.mean(spatial_entropy[mask])),
                'std_entropy': float(np.std(spatial_entropy[mask]))
            }
        print("   ✓ Spatial entropy complete")
    except Exception as e:
        print(f"   ✗ Spatial entropy failed: {e}")
        entropy_stats = {}
        entropy_by_type = {}

    print("\n3. Calculating nearest neighbor distances...")
    try:
        nn_distances = calculate_nn_distances(adata, cluster_key='cell_type')
        adata.obs['nn_distance'] = nn_distances

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        nn_df = pd.DataFrame({
            'cell_type': adata.obs['cell_type'],
            'nn_distance': nn_distances
        }).dropna()

        nn_df.boxplot(column='nn_distance', by='cell_type', ax=axes[0])
        axes[0].set_title('Nearest Neighbor Distance by Cell Type')
        axes[0].set_xlabel('Cell Type')
        axes[0].set_ylabel('Distance (pixels)')
        plt.sca(axes[0])
        plt.xticks(rotation=45, ha='right')

        axes[1].hist(nn_distances[~np.isnan(nn_distances)], bins=50, edgecolor='black')
        axes[1].set_xlabel('Nearest Neighbor Distance')
        axes[1].set_ylabel('Count')
        axes[1].axvline(np.nanmedian(nn_distances), color='red', linestyle='--', label=f'Median: {np.nanmedian(nn_distances):.1f}')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(f'{output_dir}/nearest_neighbor_distances.png', dpi=300, bbox_inches='tight')
        plt.close()

        nn_stats = {}
        for ct in cell_types:
            mask = adata.obs['cell_type'] == ct
            ct_distances = nn_distances[mask]
            ct_distances = ct_distances[~np.isnan(ct_distances)]

            if len(ct_distances) > 0:
                nn_stats[ct] = {
                    'mean_distance': float(np.mean(ct_distances)),
                    'median_distance': float(np.median(ct_distances)),
                    'std_distance': float(np.std(ct_distances)),
                    'spatial_pattern': 'tightly_clustered' if np.median(ct_distances) < 50 else 'dispersed' if np.median(ct_distances) > 150 else 'moderate'
                }
        print("   ✓ Nearest neighbor distances complete")
    except Exception as e:
        print(f"   ✗ Nearest neighbor distances failed: {e}")
        nn_stats = {}

    print("\n4. Running spatial autocorrelation (Moran's I)...")
    try:
        morans_df = calculate_morans_i_manual(adata, gene_list=None, n_neighbors=6)
        morans_df = morans_df.sort_values('I', ascending=False)

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        axes[0].hist(morans_df['I'], bins=50, edgecolor='black')
        axes[0].axvline(morans_df['I'].median(), color='red', linestyle='--', label=f"Median: {morans_df['I'].median():.3f}")
        axes[0].set_xlabel("Moran's I")
        axes[0].set_ylabel('Count')
        axes[0].set_title("Distribution of Moran's I")
        axes[0].legend()

        top_genes = morans_df.head(20)
        y_pos = np.arange(len(top_genes))
        axes[1].barh(y_pos, top_genes['I'])
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels(top_genes.index)
        axes[1].set_xlabel("Moran's I")
        axes[1].set_title('Top 20 Spatially Variable Genes')
        axes[1].invert_yaxis()

        plt.tight_layout()
        plt.savefig(f'{output_dir}/spatial_autocorrelation_extended.png', dpi=300, bbox_inches='tight')
        plt.close()

        autocorr_stats = {
            'n_genes_tested': len(morans_df),
            'n_significant': int((morans_df['I'] > 0.3).sum()),
            'mean_morans_i': float(morans_df['I'].mean()),
            'median_morans_i': float(morans_df['I'].median()),
            'top_genes': morans_df.index.tolist()[:10],
            'top_morans_i': [float(x) for x in morans_df['I'].head(10).tolist()]
        }
        print("   ✓ Spatial autocorrelation complete")
    except Exception as e:
        print(f"   ✗ Spatial autocorrelation failed: {e}")
        autocorr_stats = {}

    print("\n5. Calculating cluster compactness...")
    try:
        cluster_compactness = calculate_cluster_compactness(adata, cluster_key='leiden')

        compactness_df = pd.DataFrame(cluster_compactness).T

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        axes[0, 0].bar(compactness_df.index, compactness_df['mean_distance_to_centroid'])
        axes[0, 0].set_xlabel('Cluster')
        axes[0, 0].set_ylabel('Mean Distance to Centroid')
        axes[0, 0].set_title('Cluster Spread')
        axes[0, 0].tick_params(axis='x', rotation=45)

        axes[0, 1].bar(compactness_df.index, compactness_df['density'])
        axes[0, 1].set_xlabel('Cluster')
        axes[0, 1].set_ylabel('Density (spots/area)')
        axes[0, 1].set_title('Cluster Density')
        axes[0, 1].tick_params(axis='x', rotation=45)

        axes[1, 0].scatter(compactness_df['convex_hull_area'], compactness_df['n_spots'], s=100, alpha=0.6)
        for idx, cluster in enumerate(compactness_df.index):
            axes[1, 0].annotate(cluster, (compactness_df['convex_hull_area'].iloc[idx], compactness_df['n_spots'].iloc[idx]))
        axes[1, 0].set_xlabel('Convex Hull Area')
        axes[1, 0].set_ylabel('Number of Spots')

        axes[1, 1].bar(compactness_df.index, compactness_df['compactness_score'])
        axes[1, 1].set_xlabel('Cluster')
        axes[1, 1].set_ylabel('Compactness Score')
        axes[1, 1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(f'{output_dir}/cluster_compactness.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("   ✓ Cluster compactness complete")
    except Exception as e:
        print(f"   ✗ Cluster compactness failed: {e}")
        cluster_compactness = {}

    print("\n6. Generating comprehensive metrics JSON...")
    comprehensive_metrics = {
        'sample_info': {
            'n_spots': int(adata.shape[0]),
            'n_genes': int(adata.shape[1]),
            'n_cell_types': len(cell_types),
            'cell_types': cell_types.tolist(),
            'n_clusters': int(adata.obs['leiden'].nunique())
        },
        'neighborhood_enrichment': enrichment_summary,
        'spatial_entropy': {
            'overall': entropy_stats,
            'by_cell_type': entropy_by_type
        },
        'nearest_neighbor_distances': nn_stats,
        'spatial_autocorrelation': autocorr_stats,
        'cluster_compactness': cluster_compactness,
        'cell_type_distribution': adata.obs['cell_type'].value_counts().to_dict(),
        'cluster_distribution': adata.obs['leiden'].value_counts().to_dict()
    }

    output_path = f'{output_dir}/spatial_statistics_enhanced.json'
    with open(output_path, 'w') as f:
        json.dump(comprehensive_metrics, f, indent=2)
    print(f"   ✓ Saved to: {output_path}")

    print("\n7. Generating clinical summary...")
    clinical_summary = generate_clinical_summary(comprehensive_metrics)
    summary_path = f'{output_dir}/clinical_spatial_summary.txt'
    with open(summary_path, 'w') as f:
        f.write(clinical_summary)
    print(f"   ✓ Saved to: {summary_path}")
    print("\n" + clinical_summary)

    print("\n8. Saving updated h5ad...")
    output_h5ad = f'{output_dir}/annotated_visium_spatial_stats.h5ad'
    adata.write_h5ad(output_h5ad)
    print(f"   ✓ Saved to: {output_h5ad}")

    print("\n" + "="*50)
    print("SPATIAL STATISTICS ANALYSIS COMPLETE")
    print("="*50)
    print(f"\nGenerated files in {output_dir}:")
    print("  - spatial_statistics_enhanced.json")
    print("  - clinical_spatial_summary.txt")
    print("  - neighborhood_enrichment.png")
    print("  - spatial_entropy.png")
    print("  - nearest_neighbor_distances.png")
    print("  - spatial_autocorrelation_extended.png")
    print("  - cluster_compactness.png")
    print("  - annotated_visium_spatial_stats.h5ad")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_h5ad = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else None
        main(input_h5ad, output_dir)
    else:
        main()

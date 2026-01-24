#!/usr/bin/env python
"""Execute Scanpy baseline spatial transcriptomics analysis"""

if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    
    import scanpy as sc
    import squidpy as sq
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    import json
    from datetime import datetime
    import psutil
    import os
    
    print(f"scanpy {sc.__version__}, squidpy {sq.__version__}, anndata {ad.__version__}")
    
    SEED = 42
    np.random.seed(SEED)
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)
    
    PROJECT_ROOT = Path(__file__).parent.parent
    DATA_DIR = PROJECT_ROOT / "data" / "sample"
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print(f"Data: {DATA_DIR}, Output: {OUTPUT_DIR}")
    
    # Load
    h5_file = list(DATA_DIR.glob("*.h5"))[0]
    print(f"\n{'='*60}\nLOADING: {h5_file.name}\n{'='*60}")
    adata = sc.read_10x_h5(h5_file)
    adata.var_names_make_unique()
    print(f"Dataset: {adata.shape} spots")
    
    # Spatial coordinates
    print(f"\n{'='*60}\nSPATIAL COORDINATES\n{'='*60}")
    spatial_dir = DATA_DIR / "spatial"
    spatial_coords = pd.read_csv(
        spatial_dir / "tissue_positions_list.csv",
        header=None,
        names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    )
    spatial_coords = spatial_coords.set_index('barcode').loc[adata.obs_names]
    adata.obsm['spatial'] = spatial_coords[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
    adata.obs['in_tissue'] = spatial_coords['in_tissue'].values
    
    import PIL
    tissue_img = PIL.Image.open(spatial_dir / "tissue_hires_image.png")
    adata.uns['spatial'] = {
        'Visium_Human_Breast_Cancer': {
            'images': {'hires': np.array(tissue_img)},
            'scalefactors': {'tissue_hires_scalef': 1.0, 'spot_diameter_fullres': 89.43}
        }
    }
    print(f"Spatial: {adata.obsm['spatial'].shape}, In tissue: {adata.obs['in_tissue'].sum()}")
    
    # QC
    print(f"\n{'='*60}\nQC & FILTERING\n{'='*60}")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    print(f"Before: {adata.n_obs} spots, {adata.n_vars} genes")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After: {adata.n_obs} spots, {adata.n_vars} genes")
    
    # Normalize
    print(f"\n{'='*60}\nNORMALIZATION & PCA\n{'='*60}")
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    n_hvg = adata.var['highly_variable'].sum()
    print(f"HVG: {n_hvg}")
    
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', random_state=SEED)
    
    # Clustering
    print(f"\n{'='*60}\nCLUSTERING\n{'='*60}")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    print(f"Clusters: {len(cluster_counts)}")
    print(cluster_counts)
    
    # Spatial analysis
    print(f"\n{'='*60}\nSPATIAL ANALYSIS\n{'='*60}")
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
    print(f"Spatial graph: {adata.obsp['spatial_connectivities'].shape}")
    
    print("Computing Moran's I (100 genes, 100 permutations)...")
    top_hvgs = adata.var_names[adata.var['highly_variable']][:100]
    sq.gr.spatial_autocorr(adata, mode='moran', genes=top_hvgs, n_perms=100, n_jobs=-1)
    morans_i = adata.uns['moranI'].sort_values('I', ascending=False)
    print(f"Mean Moran's I: {morans_i['I'].mean():.4f}")
    print(f"Significant genes: {(morans_i['pval_norm'] < 0.05).sum()}")
    
    # Visualizations
    print(f"\n{'='*60}\nVISUALIZATIONS\n{'='*60}")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    sq.pl.spatial_scatter(adata, color='leiden', ax=axes[0,0], title='Clusters', size=1.5, legend_loc='right margin')
    sq.pl.spatial_scatter(adata, color='total_counts', ax=axes[0,1], title='Counts', size=1.5, cmap='viridis')
    top_gene = morans_i.index[0]
    sq.pl.spatial_scatter(adata, color=top_gene, ax=axes[1,0], title=f'{top_gene} (I={morans_i.loc[top_gene, "I"]:.3f})', size=1.5, cmap='Reds')
    sq.pl.spatial_scatter(adata, color='n_genes_by_counts', ax=axes[1,1], title='Genes', size=1.5, cmap='Blues')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "spatial_tissue_overview.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("✅ Spatial plots saved")
    
    # Co-occurrence (simplified - squidpy 1.5 has complex output)
    print("Computing co-occurrence...")
    sq.gr.co_occurrence(adata, cluster_key='leiden')

    # Create simple co-occurrence summary from cluster neighbors
    from collections import Counter
    co_occur_pairs = []
    for i in range(len(adata)):
        cluster_i = adata.obs['leiden'].iloc[i]
        neighbors = adata.obsp['spatial_connectivities'][i].nonzero()[1]
        for j in neighbors:
            cluster_j = adata.obs['leiden'].iloc[j]
            if cluster_i != cluster_j:
                pair = tuple(sorted([str(cluster_i), str(cluster_j)]))
                co_occur_pairs.append(pair)

    co_occur_counts = Counter(co_occur_pairs)
    print(f"✅ Co-occurrence computed: {len(co_occur_counts)} cluster pairs")

    # Create a simple summary plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    if len(co_occur_counts) > 0:
        pairs = list(co_occur_counts.keys())[:20]  # Top 20
        counts = [co_occur_counts[p] for p in pairs]
        pair_labels = [f"{p[0]}-{p[1]}" for p in pairs]

        ax.barh(pair_labels, counts, color='steelblue')
        ax.set_xlabel('Co-occurrence Count')
        ax.set_title('Top Cluster Co-occurrences')
        ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "spatial_cooccurrence.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("✅ Co-occurrence plot saved")

    # Export JSON
    print(f"\n{'='*60}\nEXPORTING JSON\n{'='*60}")
    cluster_info = {str(c): {
        "count": int((adata.obs['leiden'] == c).sum()),
        "mean_genes": float(adata.obs.loc[adata.obs['leiden'] == c, 'n_genes_by_counts'].mean())
    } for c in adata.obs['leiden'].unique()}

    morans_dict = {gene: {
        "morans_i": float(morans_i.loc[gene, 'I']),
        "p_value": float(morans_i.loc[gene, 'pval_norm'])
    } for gene in morans_i.head(10).index}

    # Create co-occurrence summary for JSON
    max_cooccur = []
    for cluster in adata.obs['leiden'].unique():
        neighbors = []
        cluster_spots = adata.obs['leiden'] == cluster
        for i in np.where(cluster_spots)[0]:
            neighbor_indices = adata.obsp['spatial_connectivities'][i].nonzero()[1]
            neighbors.extend(adata.obs['leiden'].iloc[neighbor_indices].tolist())

        if neighbors:
            neighbor_counts = Counter(neighbors)
            neighbor_counts.pop(cluster, None)  # Remove self
            if neighbor_counts:
                most_common = neighbor_counts.most_common(1)[0]
                max_cooccur.append({
                    "cluster": str(cluster),
                    "most_colocalized_with": str(most_common[0]),
                    "cooccurrence_count": int(most_common[1])
                })
    
    features = {
        "metadata": {"analysis_date": datetime.now().isoformat(), "scanpy_version": sc.__version__,
                     "squidpy_version": sq.__version__, "dataset": h5_file.name, "random_seed": SEED,
                     "analysis_type": "spatial_transcriptomics"},
        "dataset_summary": {"n_spots": int(adata.n_obs), "n_genes": int(adata.n_vars),
                           "n_highly_variable_genes": int(n_hvg), "spots_in_tissue": int(adata.obs['in_tissue'].sum())},
        "qc_metrics": {"mean_genes_per_spot": float(adata.obs['n_genes_by_counts'].mean()),
                      "mean_counts_per_spot": float(adata.obs['total_counts'].mean())},
        "clustering": {"n_clusters": len(cluster_counts), "resolution": 0.5, "clusters": cluster_info},
        "spatial_statistics": {
            "spatial_neighbors": {"n_neighbors": 6, "coord_type": "generic"},
            "morans_i": {"mean": float(morans_i['I'].mean()), "std": float(morans_i['I'].std()),
                        "significant_genes_count": int((morans_i['pval_norm'] < 0.05).sum()),
                        "top_genes": morans_dict},
            "spatial_cooccurrence": max_cooccur
        }
    }
    
    output_file = OUTPUT_DIR / "scanpy_features_spatial.json"
    with open(output_file, 'w') as f:
        json.dump(features, f, indent=2)
    print(f"✅ JSON: {output_file}")
    
    # Save h5ad
    output_h5ad = OUTPUT_DIR / "processed_visium.h5ad"
    adata.write(output_h5ad)
    print(f"✅ H5AD: {output_h5ad}")
    
    # Summary
    memory_mb = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024
    print(f"\n{'='*60}\nSUMMARY\n{'='*60}")
    print(f"Spots: {adata.n_obs}, Genes: {adata.n_vars}, Clusters: {len(cluster_counts)}")
    print(f"Spatial analysis: ✅")
    print(f"Mean Moran's I: {morans_i['I'].mean():.4f}")
    print(f"Memory: {memory_mb:.0f} MB")
    print(f"\n✅ COMPLETE!")
    print("="*60)

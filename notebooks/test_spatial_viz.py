#!/usr/bin/env python
"""Test spatial visualization with corrected plotting code."""

if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')

    import scanpy as sc
    import squidpy as sq
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path
    import json
    from datetime import datetime
    from PIL import Image

    print(f"scanpy version: {sc.__version__}")
    print(f"squidpy version: {sq.__version__}")
    print(f"anndata version: {ad.__version__}")

    SEED = 42
    np.random.seed(SEED)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False, figsize=(8, 6))

    PROJECT_ROOT = Path("/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma")
    DATA_DIR = PROJECT_ROOT / "data" / "sample"
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")

    print("\n" + "="*60)
    print("LOADING DATA")
    print("="*60)

    h5_file = list(DATA_DIR.glob("*.h5"))[0]
    print(f"Loading: {h5_file.name}\n")

    adata = sc.read_10x_h5(h5_file)
    adata.var_names_make_unique()

    print(f"Dataset shape: {adata.shape}")
    print(f"Spots (observations): {adata.n_obs}")
    print(f"Genes (variables): {adata.n_vars}")

    print("\n" + "="*60)
    print("LOADING SPATIAL COORDINATES")
    print("="*60)

    spatial_dir = DATA_DIR / "spatial"
    tissue_positions = pd.read_csv(
        spatial_dir / "tissue_positions_list.csv",
        header=None,
        index_col=0,
        names=['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    )

    common_barcodes = adata.obs_names.intersection(tissue_positions.index)
    adata = adata[common_barcodes]
    tissue_positions = tissue_positions.loc[common_barcodes]

    spatial_coords = tissue_positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
    adata.obsm['spatial'] = spatial_coords

    print(f"Spatial coordinates loaded: {adata.obsm['spatial'].shape}")
    print(f"Spots in tissue: {adata.n_obs} / {len(tissue_positions)}")

    tissue_img = Image.open(spatial_dir / "tissue_hires_image.png")
    with open(spatial_dir / "scalefactors_json.json", 'r') as f:
        scalefactors = json.load(f)

    library_id = 'breast_cancer_visium'
    adata.uns['spatial'] = {
        library_id: {
            'images': {'hires': np.array(tissue_img)},
            'scalefactors': {
                'tissue_hires_scalef': scalefactors['tissue_hires_scalef'],
                'spot_diameter_fullres': scalefactors['spot_diameter_fullres']
            }
        }
    }

    print(f"Library ID: {library_id}")
    print(f"Tissue image shape: {adata.uns['spatial'][library_id]['images']['hires'].shape}")

    print("\n" + "="*60)
    print("QUALITY CONTROL")
    print("="*60)

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    print(f"Mean genes per spot: {int(adata.obs['n_genes_by_counts'].mean())}")
    print(f"Mean counts per spot: {int(adata.obs['total_counts'].mean())}")
    print(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")

    print("\n" + "="*60)
    print("FILTERING")
    print("="*60)

    print(f"Before filtering: {adata.n_obs} spots, {adata.n_vars} genes")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After filtering: {adata.n_obs} spots, {adata.n_vars} genes")

    print("\n" + "="*60)
    print("NORMALIZATION & PREPROCESSING")
    print("="*60)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

    adata.raw = adata
    adata = adata[:, adata.var['highly_variable']]

    sc.pp.scale(adata, max_value=10)

    print("\n" + "="*60)
    print("DIMENSIONALITY REDUCTION")
    print("="*60)

    sc.tl.pca(adata, n_comps=50, random_state=SEED)
    sc.pp.neighbors(adata, n_pcs=40, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)

    print("\n" + "="*60)
    print("CLUSTERING")
    print("="*60)

    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    print(f"\nClusters found: {adata.obs['leiden'].nunique()}")
    print(f"\nCluster sizes:\n{adata.obs['leiden'].value_counts()}")

    print("\n" + "="*60)
    print("SPATIAL NEIGHBOR GRAPH & MORAN'S I")
    print("="*60)

    print("Building spatial neighbor graph...")
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

    print(f"Spatial connectivities: {adata.obsp['spatial_connectivities'].shape}")
    print(f"Spatial distances: {adata.obsp['spatial_distances'].shape}")

    print("\nComputing Moran's I for spatial autocorrelation...")
    top_hvgs = adata.var_names[:100]

    sq.gr.spatial_autocorr(
        adata,
        mode='moran',
        genes=top_hvgs,
        n_perms=100,
        n_jobs=1
    )

    morans_df = adata.uns['moranI']
    significant_genes = morans_df[morans_df['pval_norm'] < 0.05].sort_values('I', ascending=False)

    print(f"\nSignificant spatial genes (p<0.05): {len(significant_genes)}")
    if len(significant_genes) > 0:
        print(f"Top spatial gene: {significant_genes.index[0]} (Moran's I = {significant_genes.iloc[0]['I']:.2f})")

    print("\n" + "="*60)
    print("SPATIAL VISUALIZATION ON TISSUE - CORRECTED")
    print("="*60)

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='leiden',
        ax=axes[0, 0],
        title='Spatial Clusters on Tissue',
        size=1.5,
        img=True,
        img_res_key='hires',
        img_alpha=0.5,
        alpha=0.8,
        legend_loc='right margin',
        frameon=False
    )

    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='total_counts',
        ax=axes[0, 1],
        title='Total Counts per Spot',
        size=1.5,
        img=True,
        img_res_key='hires',
        img_alpha=0.5,
        alpha=0.8,
        cmap='viridis',
        frameon=False
    )

    if len(significant_genes) > 0:
        top_gene = significant_genes.index[0]
        sq.pl.spatial_scatter(
            adata,
            library_id=library_id,
            color=top_gene,
            ax=axes[1, 0],
            title=f'Top Spatial Gene: {top_gene}',
            size=1.5,
            img=True,
            img_res_key='hires',
            img_alpha=0.5,
            alpha=0.8,
            cmap='Reds',
            frameon=False
        )

    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='n_genes_by_counts',
        ax=axes[1, 1],
        title='Genes Detected per Spot',
        size=1.5,
        img=True,
        img_res_key='hires',
        img_alpha=0.5,
        alpha=0.8,
        cmap='plasma',
        frameon=False
    )

    plt.tight_layout()
    output_path = OUTPUT_DIR / "spatial_tissue_overview_CORRECTED.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_path}")
    plt.close()

    print("\n" + "="*60)
    print("SPATIAL CO-OCCURRENCE ANALYSIS")
    print("="*60)

    co_occur_pairs = []
    for i in range(len(adata)):
        cluster_i = adata.obs['leiden'].iloc[i]
        neighbors = adata.obsp['spatial_connectivities'][i].nonzero()[1]
        for j in neighbors:
            cluster_j = adata.obs['leiden'].iloc[j]
            if cluster_i != cluster_j:
                pair = tuple(sorted([str(cluster_i), str(cluster_j)]))
                co_occur_pairs.append(pair)

    from collections import Counter
    pair_counts = Counter(co_occur_pairs)

    co_occur_data = []
    for cluster in adata.obs['leiden'].unique():
        cluster_pairs = [(pair, count) for pair, count in pair_counts.items() if str(cluster) in pair]
        if cluster_pairs:
            max_pair = max(cluster_pairs, key=lambda x: x[1])
            other_cluster = [c for c in max_pair[0] if c != str(cluster)][0]
            co_occur_data.append({
                'cluster': str(cluster),
                'most_colocalized_with': other_cluster,
                'cooccurrence_count': max_pair[1]
            })

    print(f"Co-occurrence patterns identified: {len(co_occur_data)}")

    print("\n" + "="*60)
    print("EXPORTING SPATIAL FEATURES TO JSON")
    print("="*60)

    cluster_info = {}
    for cluster in adata.obs['leiden'].unique():
        cluster_spots = adata[adata.obs['leiden'] == cluster]
        cluster_info[str(cluster)] = {
            'count': int(len(cluster_spots)),
            'mean_genes': float(cluster_spots.obs['n_genes_by_counts'].mean())
        }

    morans_top = {}
    for gene in significant_genes.head(10).index:
        morans_top[gene] = {
            'morans_i': float(significant_genes.loc[gene, 'I']),
            'p_value': float(significant_genes.loc[gene, 'pval_norm'])
        }

    features = {
        'metadata': {
            'analysis_date': datetime.now().isoformat(),
            'scanpy_version': sc.__version__,
            'squidpy_version': sq.__version__,
            'dataset': h5_file.name,
            'random_seed': SEED,
            'analysis_type': 'spatial_transcriptomics'
        },
        'dataset_summary': {
            'n_spots': int(adata.n_obs),
            'n_genes': int(adata.raw.n_vars),
            'n_highly_variable_genes': 2000,
            'spots_in_tissue': int(adata.n_obs)
        },
        'qc_metrics': {
            'mean_genes_per_spot': float(adata.obs['n_genes_by_counts'].mean()),
            'mean_counts_per_spot': float(adata.obs['total_counts'].mean())
        },
        'clustering': {
            'n_clusters': int(adata.obs['leiden'].nunique()),
            'resolution': 0.5,
            'clusters': cluster_info
        },
        'spatial_statistics': {
            'spatial_neighbors': {
                'n_neighbors': 6,
                'coord_type': 'generic'
            },
            'morans_i': {
                'mean': float(morans_df['I'].mean()),
                'std': float(morans_df['I'].std()),
                'significant_genes_count': int(len(significant_genes)),
                'top_genes': morans_top
            },
            'spatial_cooccurrence': co_occur_data
        }
    }

    json_path = OUTPUT_DIR / "scanpy_features_spatial_CORRECTED.json"
    with open(json_path, 'w') as f:
        json.dump(features, f, indent=2)

    print(f"\nSaved spatial features: {json_path}")

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nOutputs generated:")
    print(f"1. {output_path}")
    print(f"2. {json_path}")
    print("\nCorrected visualization should now show clusters overlaid on tissue!")

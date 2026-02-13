#!/usr/bin/env python
"""Test proper Visium loading with sc.read_visium()."""

if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')

    import scanpy as sc
    import squidpy as sq
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    print(f"scanpy version: {sc.__version__}")
    print(f"squidpy version: {sq.__version__}")

    SEED = 42
    np.random.seed(SEED)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

    PROJECT_ROOT = Path("/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma")
    DATA_DIR = PROJECT_ROOT / "data" / "sample"
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("\n" + "="*60)
    print("LOADING VISIUM DATA PROPERLY")
    print("="*60)

    # USE SCANPY'S BUILT-IN VISIUM READER!
    adata = sc.read_visium(
        DATA_DIR,
        count_file='Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5',
        load_images=True
    )
    adata.var_names_make_unique()

    print(f"Dataset loaded: {adata.shape}")
    print(f"Spatial coordinates: {adata.obsm['spatial'].shape}")

    # Get the library_id that was auto-assigned
    library_id = list(adata.uns['spatial'].keys())[0]
    print(f"Library ID: {library_id}")

    print("\n" + "="*60)
    print("QUALITY CONTROL")
    print("="*60)

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    print(f"Mean genes per spot: {int(adata.obs['n_genes_by_counts'].mean())}")
    print(f"Mean counts per spot: {int(adata.obs['total_counts'].mean())}")

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
    print("DIMENSIONALITY REDUCTION & CLUSTERING")
    print("="*60)

    sc.tl.pca(adata, n_comps=50, random_state=SEED)
    sc.pp.neighbors(adata, n_pcs=40, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    print(f"Clusters found: {adata.obs['leiden'].nunique()}")

    print("\n" + "="*60)
    print("SPATIAL NEIGHBOR GRAPH")
    print("="*60)

    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
    print(f"Spatial graph built: {adata.obsp['spatial_connectivities'].shape}")

    print("\n" + "="*60)
    print("MORAN'S I")
    print("="*60)

    top_hvgs = adata.var_names[:100]
    sq.gr.spatial_autocorr(adata, mode='moran', genes=top_hvgs, n_perms=100, n_jobs=1)

    morans_df = adata.uns['moranI'].sort_values('I', ascending=False)
    significant_genes = morans_df[morans_df['pval_norm'] < 0.05]

    print(f"Significant spatial genes: {len(significant_genes)}")
    if len(significant_genes) > 0:
        print(f"Top gene: {significant_genes.index[0]} (Moran's I = {significant_genes.iloc[0]['I']:.2f})")

    print("\n" + "="*60)
    print("SPATIAL VISUALIZATION - PROPER ALIGNMENT")
    print("="*60)

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    # Test with Scanpy's spatial plotting first
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='leiden',
        ax=axes[0, 0],
        title='Scanpy: Clusters on Tissue',
        size=1.5,
        show=False
    )

    # Now try Squidpy
    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='leiden',
        ax=axes[0, 1],
        title='Squidpy: Clusters on Tissue',
        size=1.5,
        img=True,
        img_res_key='hires',
        img_alpha=0.5,
        alpha=0.8,
        frameon=False
    )

    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='total_counts',
        ax=axes[1, 0],
        title='Scanpy: Total Counts',
        size=1.5,
        cmap='viridis',
        show=False
    )

    if len(significant_genes) > 0:
        top_gene = significant_genes.index[0]
        sc.pl.spatial(
            adata,
            library_id=library_id,
            color=top_gene,
            ax=axes[1, 1],
            title=f'Scanpy: {top_gene}',
            size=1.5,
            cmap='Reds',
            show=False
        )

    plt.tight_layout()
    output_path = OUTPUT_DIR / "spatial_tissue_PROPER.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    print("\n" + "="*60)
    print("COMPLETE!")
    print("="*60)
    print(f"\nCheck {output_path} - spots should now be ON the tissue, not around it!")

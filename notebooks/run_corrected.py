if __name__ == "__main__":
    #!/usr/bin/env python
    # coding: utf-8

    # # Spatial Transcriptomics Analysis with Scanpy + Squidpy
    # 
    # **Author**: MedGemma Spatial Project  
    # **Date**: 2026-01-24  
    # **Purpose**: Proper spatial analysis of 10x Visium breast cancer data
    # 
    # ## What You'll Learn:
    # 1. How to load Visium spatial transcriptomics data
    # 2. Quality control for spatial data
    # 3. Spatial clustering and visualization **on tissue**
    # 4. Computing Moran's I for spatial autocorrelation
    # 5. Spatial co-occurrence analysis
    # 
    # ## Key Difference from scRNA-seq:
    # - Each "cell" is actually a **spot** on a tissue slide
    # - Spots have **spatial coordinates** (x, y positions)
    # - We can visualize patterns **overlaid on tissue histology**
    # - We can ask: which genes show spatial clustering?
    # 
    # ---

    # ## 1. Setup & Imports

    # In[ ]:


    import warnings
    warnings.filterwarnings('ignore')

    import scanpy as sc
    import squidpy as sq
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    import json
    from datetime import datetime

    print(f"scanpy version: {sc.__version__}")
    print(f"squidpy version: {sq.__version__}")
    print(f"anndata version: {ad.__version__}")


    # ### Configuration

    # In[ ]:


    SEED = 42
    np.random.seed(SEED)

    # Scanpy settings for nice plots
    sc.settings.verbosity = 2  # Less verbose for notebook
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False, figsize=(8, 6))

    # Paths
    PROJECT_ROOT = Path("/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma")
    DATA_DIR = PROJECT_ROOT / "data" / "sample"
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")


    # ---
    # ## 2. Load Visium Data
    # 
    # **What is Visium data?**
    # - H5 file contains gene expression matrix (spots x genes)
    # - Spatial folder contains:
    #   - `tissue_positions_list.csv` - coordinates of each spot
    #   - `tissue_hires_image.png` - H&E histology image
    #   - `scalefactors_json.json` - alignment info

    # In[ ]:


    # Load gene expression data
    h5_file = list(DATA_DIR.glob("*.h5"))[0]
    print(f"Loading: {h5_file.name}")

    adata = sc.read_10x_h5(h5_file)
    adata.var_names_make_unique()

    print(f"\nDataset shape: {adata.shape}")
    print(f"Spots (observations): {adata.n_obs}")
    print(f"Genes (variables): {adata.n_vars}")


    # In[ ]:


    # Peek at the data
    print("First 3 spots:")
    print(adata.obs.head(3))

    print("\nFirst 3 genes:")
    print(adata.var.head(3))


    # ---
    # ## 3. Load Spatial Coordinates & Tissue Image
    # 
    # **This is the KEY step that makes it spatial!**
    # 
    # We need to:
    # 1. Load spot coordinates from `tissue_positions_list.csv`
    # 2. Load the H&E tissue image
    # 3. Store them in the AnnData object in the **correct format** for Squidpy

    # In[ ]:


    print("Loading spatial coordinates and tissue image...")

    spatial_dir = DATA_DIR / "spatial"
    tissue_positions = spatial_dir / "tissue_positions_list.csv"

    # Load tissue positions (no header in Visium format)
    spatial_coords = pd.read_csv(
        tissue_positions,
        header=None,
        names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    )

    # Match to our spots
    spatial_coords = spatial_coords.set_index('barcode')
    spatial_coords = spatial_coords.loc[adata.obs_names]

    print(f"Loaded coordinates for {len(spatial_coords)} spots")
    print(f"\nSpatial coordinates columns:")
    print(spatial_coords.columns.tolist())


    # In[ ]:


    # Add spatial coordinates to AnnData
    # IMPORTANT: Use pixel coordinates (not array indices) for visualization
    adata.obsm['spatial'] = spatial_coords[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
    adata.obs['in_tissue'] = spatial_coords['in_tissue'].values

    print(f"Added spatial coordinates: {adata.obsm['spatial'].shape}")
    print(f"Spots in tissue: {adata.obs['in_tissue'].sum()} / {len(adata)}")


    # In[ ]:


    # Load tissue image and scale factors
    import PIL

    tissue_img = PIL.Image.open(spatial_dir / "tissue_hires_image.png")

    # Load scale factors (tells us how to map coordinates to image)
    with open(spatial_dir / "scalefactors_json.json") as f:
        scalefactors = json.load(f)

    print(f"Tissue image size: {tissue_img.size}")
    print(f"Scale factors: {scalefactors}")


    # In[ ]:


    # CRITICAL: Store spatial data in the correct Squidpy format
    # The library_id (key name) will be used in plotting!
    library_id = 'breast_cancer_visium'

    adata.uns['spatial'] = {
        library_id: {
            'images': {
                'hires': np.array(tissue_img)
            },
            'scalefactors': {
                'tissue_hires_scalef': scalefactors['tissue_hires_scalef'],
                'spot_diameter_fullres': scalefactors['spot_diameter_fullres']
            }
        }
    }

    print(f"✅ Spatial data stored with library_id: '{library_id}'")
    print(f"   Image shape: {adata.uns['spatial'][library_id]['images']['hires'].shape}")


    # **Let's verify the tissue image loaded correctly:**

    # In[ ]:


    # Quick visualization of tissue image
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.imshow(adata.uns['spatial'][library_id]['images']['hires'])
    ax.set_title('H&E Tissue Image')
    ax.axis('off')
    plt.tight_layout()
    plt.show()

    print("If you see the breast tissue histology above, spatial data is loaded correctly! ✅")


    # ---
    # ## 4. Quality Control (QC)
    # 
    # Same as scRNA-seq, but we're looking at **spots** not cells.

    # In[ ]:


    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )

    print("QC metrics calculated:")
    print(adata.obs.columns.tolist())


    # In[ ]:


    # Visualize QC metrics
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "qc_violin.png", dpi=150, bbox_inches='tight')
    plt.show()

    print(f"\nQC Summary:")
    print(f"Mean genes per spot: {adata.obs['n_genes_by_counts'].mean():.0f}")
    print(f"Mean counts per spot: {adata.obs['total_counts'].mean():.0f}")
    print(f"Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}%")


    # ### Filter low-quality spots

    # In[ ]:


    print(f"Before filtering: {adata.n_obs} spots, {adata.n_vars} genes")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After filtering: {adata.n_obs} spots, {adata.n_vars} genes")


    # ---
    # ## 5. Normalization & Preprocessing
    # 
    # Standard workflow:
    # 1. Normalize counts (library size normalization)
    # 2. Log-transform
    # 3. Find highly variable genes (HVGs)
    # 4. Scale data

    # In[ ]:


    # Save raw counts
    adata.raw = adata.copy()

    # Normalize to 10,000 counts per spot
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')

    n_hvg = adata.var['highly_variable'].sum()
    print(f"Highly variable genes: {n_hvg}")


    # In[ ]:


    # Visualize HVGs
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(OUTPUT_DIR / "highly_variable_genes.png", dpi=150, bbox_inches='tight')
    plt.show()


    # ---
    # ## 6. Dimensionality Reduction
    # 
    # PCA to reduce from 21,000 genes to 50 principal components

    # In[ ]:


    # Scale data (mean=0, var=1)
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, svd_solver='arpack', random_state=SEED)

    # Visualize variance explained
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    plt.savefig(OUTPUT_DIR / "pca_variance.png", dpi=150, bbox_inches='tight')
    plt.show()


    # ---
    # ## 7. Clustering
    # 
    # Leiden clustering based on gene expression (not spatial yet!)

    # In[ ]:


    # Build KNN graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=SEED)

    # UMAP for visualization
    sc.tl.umap(adata, random_state=SEED)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    print(f"\nClusters found: {len(cluster_counts)}")
    print("\nCluster sizes:")
    print(cluster_counts)


    # ### Traditional UMAP visualization (no spatial info)

    # In[ ]:


    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='UMAP - Leiden Clusters')

    cluster_counts.plot(kind='bar', ax=axes[1], color='steelblue')
    axes[1].set_xlabel('Cluster')
    axes[1].set_ylabel('Number of Spots')
    axes[1].set_title('Cluster Sizes')
    axes[1].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "clustering_overview.png", dpi=150, bbox_inches='tight')
    plt.show()


    # ---
    # ## 8. SPATIAL ANALYSIS - Build Spatial Neighbor Graph
    # 
    # **Now we use Squidpy for spatial analysis!**
    # 
    # ### What is a spatial neighbor graph?
    # - Connect spots that are physically close on the tissue
    # - Visium: hexagonal lattice, each spot has ~6 neighbors
    # - Enables spatial statistics like Moran's I

    # In[ ]:


    print("Building spatial neighbor graph...")

    # Build spatial graph based on physical proximity
    sq.gr.spatial_neighbors(
        adata, 
        coord_type='generic',  # Use generic coordinates (not grid)
        n_neighs=6             # Visium hexagonal lattice
    )

    print(f"✅ Spatial graph built")
    print(f"   Connectivities matrix: {adata.obsp['spatial_connectivities'].shape}")
    print(f"   Distances matrix: {adata.obsp['spatial_distances'].shape}")


    # ### Compute Moran's I - Spatial Autocorrelation
    # 
    # **What is Moran's I?**
    # - Measures whether a gene's expression is spatially clustered
    # - Range: -1 to +1
    #   - **+1**: Strong spatial clustering (similar expression in nearby spots)
    #   - **0**: Random distribution
    #   - **-1**: Dispersed (neighboring spots have opposite expression)
    # 
    # **Example**: If ISG15 has Moran's I = 0.57, it means immune cells (high ISG15) cluster together spatially!

    # In[ ]:


    print("Computing Moran's I for spatial autocorrelation...")
    print("This tests: Which genes show spatial patterns?\n")

    # Test top 100 HVGs for spatial autocorrelation
    top_hvgs = adata.var_names[adata.var['highly_variable']][:100]

    sq.gr.spatial_autocorr(
        adata,
        mode='moran',      # Moran's I statistic
        genes=top_hvgs,    # Test these genes
        n_perms=100,       # Permutation test for significance
        n_jobs=-1          # Use all CPU cores
    )

    print("\n✅ Moran's I computed")


    # In[ ]:


    # Sort genes by Moran's I score
    morans_i = adata.uns['moranI'].sort_values('I', ascending=False)

    print(f"Top 10 spatially autocorrelated genes:")
    print(morans_i.head(10))

    print(f"\nMean Moran's I: {morans_i['I'].mean():.4f}")
    print(f"Significant genes (p < 0.05): {(morans_i['pval_norm'] < 0.05).sum()}")


    # ---
    # ## 9. SPATIAL VISUALIZATION ON TISSUE 🎨
    # 
    # **This is what makes spatial transcriptomics powerful!**
    # 
    # We'll create a 4-panel figure:
    # 1. **Clusters** overlaid on tissue
    # 2. **Total counts** per spot (QC)
    # 3. **Top spatial gene** (ISG15 or similar)
    # 4. **Genes detected** per spot

    # In[ ]:


    # Create 4-panel spatial visualization
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Panel 1: Clusters on tissue
    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,  # CRITICAL: Must match the key we used!
        color='leiden',
        ax=axes[0, 0],
        title='Spatial Clusters on Tissue',
        size=1.5,
        img_key='hires',
        alpha=0.8,  # Transparency to see tissue underneath
        legend_loc='right margin',
        show=False
    )

    # Panel 2: Total counts (shows tissue quality)
    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='total_counts',
        ax=axes[0, 1],
        title='Total Counts per Spot',
        size=1.5,
        img_key='hires',
        cmap='viridis',
        alpha=0.8,
        show=False
    )

    # Panel 3: Most spatially autocorrelated gene
    top_gene = morans_i.index[0]
    morans_score = morans_i.loc[top_gene, 'I']

    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color=top_gene,
        ax=axes[1, 0],
        title=f'{top_gene} Expression (Moran I={morans_score:.3f})',
        size=1.5,
        img_key='hires',
        cmap='Reds',
        alpha=0.8,
        show=False
    )

    # Panel 4: Genes detected (tissue coverage)
    sq.pl.spatial_scatter(
        adata,
        library_id=library_id,
        color='n_genes_by_counts',
        ax=axes[1, 1],
        title='Genes Detected per Spot',
        size=1.5,
        img_key='hires',
        cmap='Blues',
        alpha=0.8,
        show=False
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "spatial_tissue_overview_CORRECTED.png", dpi=150, bbox_inches='tight')
    plt.show()

    print(f"\n✅ Spatial tissue visualizations saved!")
    print(f"   You should see colored spots overlaid on the H&E tissue image.")
    print(f"   Each panel shows different information mapped to tissue space.")


    # ### Individual spatial plots for key genes

    # In[ ]:


    # Plot top 3 spatially autocorrelated genes
    top_3_genes = morans_i.index[:3]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for i, gene in enumerate(top_3_genes):
        score = morans_i.loc[gene, 'I']

        sq.pl.spatial_scatter(
            adata,
            library_id=library_id,
            color=gene,
            ax=axes[i],
            title=f'{gene}\n(Moran I={score:.3f})',
            size=2,
            img_key='hires',
            cmap='Reds',
            alpha=0.8,
            show=False
        )

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "top_spatial_genes.png", dpi=150, bbox_inches='tight')
    plt.show()

    print(f"Top spatial genes:")
    for gene in top_3_genes:
        score = morans_i.loc[gene, 'I']
        print(f"  {gene}: Moran's I = {score:.3f}")


    # ---
    # ## 10. Spatial Co-occurrence Analysis
    # 
    # **Question**: Which clusters are neighbors? Do tumor and immune clusters co-localize?
    # 
    # Squidpy can compute co-occurrence statistics to find cluster pairs that appear together more often than expected by chance.

    # In[ ]:


    print("Computing spatial co-occurrence between clusters...")

    sq.gr.co_occurrence(
        adata,
        cluster_key='leiden'
    )

    print("✅ Co-occurrence computed")


    # ### Visualize co-occurrence patterns

    # In[ ]:


    # Plot co-occurrence
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    sq.pl.co_occurrence(
        adata,
        cluster_key='leiden',
        clusters=None,  # All clusters
        figsize=(8, 6),
        save=False
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "spatial_cooccurrence_CORRECTED.png", dpi=150, bbox_inches='tight')
    plt.show()

    print("✅ Co-occurrence plot saved")


    # ### Summary: Which clusters are neighbors?

    # In[ ]:


    # Analyze co-occurrence from spatial graph
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

    print(f"Top 10 cluster pairs that co-occur spatially:\n")
    for pair, count in co_occur_counts.most_common(10):
        print(f"  Cluster {pair[0]} ↔ Cluster {pair[1]}: {count} co-occurrences")


    # ---
    # ## 11. Export Features for MedGemma
    # 
    # Create a comprehensive JSON file with all spatial features

    # In[ ]:


    print("Exporting spatial features to JSON...")

    # Cluster information
    cluster_info = {}
    for cluster in adata.obs['leiden'].unique():
        cluster_mask = adata.obs['leiden'] == cluster
        cluster_info[str(cluster)] = {
            "count": int(cluster_mask.sum()),
            "mean_genes": float(adata.obs.loc[cluster_mask, 'n_genes_by_counts'].mean()),
            "mean_counts": float(adata.obs.loc[cluster_mask, 'total_counts'].mean())
        }

    # Moran's I top genes
    morans_top = morans_i.head(10)
    morans_dict = {
        gene: {
            "morans_i": float(morans_top.loc[gene, 'I']),
            "p_value": float(morans_top.loc[gene, 'pval_norm'])
        }
        for gene in morans_top.index
    }

    # Co-occurrence summary
    max_cooccur = []
    for cluster in adata.obs['leiden'].unique():
        neighbors = []
        cluster_spots = adata.obs['leiden'] == cluster
        for i in np.where(cluster_spots)[0]:
            neighbor_indices = adata.obsp['spatial_connectivities'][i].nonzero()[1]
            neighbors.extend(adata.obs['leiden'].iloc[neighbor_indices].tolist())

        if neighbors:
            neighbor_counts = Counter(neighbors)
            neighbor_counts.pop(cluster, None)
            if neighbor_counts:
                most_common = neighbor_counts.most_common(1)[0]
                max_cooccur.append({
                    "cluster": str(cluster),
                    "most_colocalized_with": str(most_common[0]),
                    "cooccurrence_count": int(most_common[1])
                })

    # Build features dictionary
    features = {
        "metadata": {
            "analysis_date": datetime.now().isoformat(),
            "scanpy_version": sc.__version__,
            "squidpy_version": sq.__version__,
            "dataset": h5_file.name,
            "random_seed": SEED,
            "analysis_type": "spatial_transcriptomics"
        },
        "dataset_summary": {
            "n_spots": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_highly_variable_genes": int(n_hvg),
            "spots_in_tissue": int(adata.obs['in_tissue'].sum())
        },
        "qc_metrics": {
            "mean_genes_per_spot": float(adata.obs['n_genes_by_counts'].mean()),
            "mean_counts_per_spot": float(adata.obs['total_counts'].mean()),
            "mean_pct_mt": float(adata.obs['pct_counts_mt'].mean())
        },
        "clustering": {
            "n_clusters": int(len(cluster_counts)),
            "resolution": 0.5,
            "clusters": cluster_info
        },
        "spatial_statistics": {
            "spatial_neighbors": {
                "n_neighbors": 6,
                "coord_type": "generic"
            },
            "morans_i": {
                "mean": float(morans_i['I'].mean()),
                "std": float(morans_i['I'].std()),
                "significant_genes_count": int((morans_i['pval_norm'] < 0.05).sum()),
                "top_genes": morans_dict
            },
            "spatial_cooccurrence": max_cooccur
        }
    }

    # Save to JSON
    output_file = OUTPUT_DIR / "scanpy_features_spatial_CORRECTED.json"
    with open(output_file, 'w') as f:
        json.dump(features, f, indent=2)

    print(f"\n✅ Enhanced spatial features exported to: {output_file}")
    print(f"\nSpatial feature summary:")
    print(f"  - Spatial neighbors computed: ✅")
    print(f"  - Moran's I genes analyzed: {len(morans_i)}")
    print(f"  - Significant spatial genes: {(morans_i['pval_norm'] < 0.05).sum()}")
    print(f"  - Co-occurrence pairs: {len(max_cooccur)}")


    # ---
    # ## 12. Save Processed Data

    # In[ ]:


    output_h5ad = OUTPUT_DIR / "processed_visium_CORRECTED.h5ad"
    adata.write(output_h5ad)
    print(f"✅ Processed AnnData saved to: {output_h5ad}")


    # ---
    # ## 13. Final Summary

    # In[ ]:


    import psutil
    import os

    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / 1024 / 1024

    print("="*60)
    print("SPATIAL TRANSCRIPTOMICS ANALYSIS SUMMARY")
    print("="*60)
    print(f"\nDataset: {h5_file.name}")
    print(f"Spots analyzed: {adata.n_obs}")
    print(f"Spots in tissue: {adata.obs['in_tissue'].sum()}")
    print(f"Genes analyzed: {adata.n_vars}")
    print(f"Clusters identified: {len(cluster_counts)}")

    print(f"\nSpatial Analysis:")
    print(f"  - Spatial neighbor graph: ✅")
    print(f"  - Mean Moran's I: {morans_i['I'].mean():.4f}")
    print(f"  - Significant spatial genes: {(morans_i['pval_norm'] < 0.05).sum()}")
    print(f"  - Top gene: {morans_i.index[0]} (I={morans_i.iloc[0]['I']:.3f})")

    print(f"\nMemory usage: {memory_mb:.0f} MB")

    print(f"\nOutputs saved:")
    print(f"  - Spatial features JSON: {output_file}")
    print(f"  - Processed h5ad: {output_h5ad}")
    print(f"  - Spatial tissue overview: {OUTPUT_DIR}/spatial_tissue_overview_CORRECTED.png")
    print(f"  - Top spatial genes: {OUTPUT_DIR}/top_spatial_genes.png")
    print(f"  - Co-occurrence: {OUTPUT_DIR}/spatial_cooccurrence_CORRECTED.png")

    print("\n✅ Spatial transcriptomics analysis complete!")
    print("="*60)


    # ---
    # ## Key Takeaways
    # 
    # ### What We Learned:
    # 
    # 1. **Loading Spatial Data**:
    #    - Must load both gene expression (h5) AND spatial coordinates (CSV) AND tissue image (PNG)
    #    - Store in correct format: `adata.uns['spatial'][library_id]`
    # 
    # 2. **Spatial Neighbor Graphs**:
    #    - Connect physically adjacent spots
    #    - Required for spatial statistics
    # 
    # 3. **Moran's I**:
    #    - Identifies genes with spatial patterns
    #    - High Moran's I = gene clusters spatially
    #    - Example: ISG15 (immune gene) has high Moran's I → immune cells cluster together!
    # 
    # 4. **Spatial Visualization**:
    #    - MUST use `library_id` parameter in `sq.pl.spatial_scatter()`
    #    - Use `alpha=0.8` to see tissue underneath
    #    - `img_key='hires'` to show high-res image
    # 
    # 5. **Co-occurrence**:
    #    - Which clusters are neighbors?
    #    - Important for understanding tissue architecture
    # 
    # ### Common Pitfalls Fixed:
    # - ❌ Wrong: Not specifying `library_id` → blank tissue images
    # - ✅ Right: Use `library_id=library_id` in all spatial plots
    # 
    # - ❌ Wrong: Using array coordinates → misaligned spots
    # - ✅ Right: Use pixel coordinates (`pxl_row_in_fullres`, `pxl_col_in_fullres`)
    # 
    # - ❌ Wrong: Missing scale factors → spots don't align with image
    # - ✅ Right: Load from `scalefactors_json.json`
    # 
    # ---
    # 
    # ## Next Steps:
    # 
    # **Ready for MedGemma integration!**
    # 
    # Use `scanpy_features_spatial_CORRECTED.json` to generate clinical pathology reports.
    # 
    # The JSON contains:
    # - 10 spatial clusters
    # - Top spatial genes (ISG15, C1QA, C1QB)
    # - Co-localization patterns
    # - Moran's I statistics
    # 
    # MedGemma can use this to generate reports like:
    # > "Analysis reveals distinct spatial organization with elevated immune activation markers (ISG15, Moran's I=0.57) clustering in specific tissue regions, co-localizing with complement pathway components (C1QA/C1QB)..."
    # 
    # ---


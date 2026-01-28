#!/usr/bin/env python
"""Deep dive into cell type analysis with marker validation and interface detection."""

if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')

    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    from scipy.sparse import csr_matrix

    print(f"scanpy version: {sc.__version__}")

    SEED = 42
    np.random.seed(SEED)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

    PROJECT_ROOT = Path("/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma")
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("\n" + "="*60)
    print("LOADING ANNOTATED DATA")
    print("="*60)

    h5ad_path = OUTPUT_DIR / "annotated_visium.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Annotated data not found at {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded: {adata.shape}")
    print(f"Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    library_id = list(adata.uns['spatial'].keys())[0]

    print("\n" + "="*60)
    print("MARKER GENE VALIDATION")
    print("="*60)

    marker_genes = {
        'Luminal': ['ESR1', 'PGR', 'KRT8', 'KRT18', 'GATA3', 'SCGB2A2'],
        'Plasma_B': ['CD79A', 'CD79B', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'MZB1', 'SDC1', 'JCHAIN']
    }

    available_markers = {}
    for category, genes in marker_genes.items():
        available = [g for g in genes if g in adata.raw.var_names]
        available_markers[category] = available
        print(f"\n{category} markers available: {len(available)}/{len(genes)}")
        if available:
            print(f"  {', '.join(available)}")

    all_available_markers = []
    for genes in available_markers.values():
        all_available_markers.extend(genes[:3])

    if len(all_available_markers) > 0:
        print("\nCreating marker gene spatial plots...")
        n_markers = min(len(all_available_markers), 6)
        markers_to_plot = all_available_markers[:n_markers]

        ncols = 3
        nrows = int(np.ceil(n_markers / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5*nrows))
        if n_markers == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for idx, gene in enumerate(markers_to_plot):
            sc.pl.spatial(
                adata,
                library_id=library_id,
                color=gene,
                ax=axes[idx],
                title=f'{gene}',
                size=1.3,
                cmap='Reds',
                show=False
            )

        for idx in range(n_markers, len(axes)):
            axes[idx].axis('off')

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / 'marker_genes_spatial.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {OUTPUT_DIR / 'marker_genes_spatial.png'}")
        plt.close()

        print("\nQuantitative marker validation:")
        validation_results = []

        for cell_type in adata.obs['cell_type'].unique():
            mask = adata.obs['cell_type'] == cell_type

            if len(available_markers['Luminal']) > 0:
                luminal_expr = adata.raw[mask, available_markers['Luminal']].X.mean()
            else:
                luminal_expr = 0

            if len(available_markers['Plasma_B']) > 0:
                plasma_expr = adata.raw[mask, available_markers['Plasma_B']].X.mean()
            else:
                plasma_expr = 0

            validation_results.append({
                'Cell Type': cell_type,
                'Luminal (mean)': f"{luminal_expr:.3f}",
                'Plasma (mean)': f"{plasma_expr:.3f}",
                'Spots': int(mask.sum())
            })

        validation_df = pd.DataFrame(validation_results)
        print(validation_df.to_string(index=False))

    print("\n" + "="*60)
    print("TUMOR-IMMUNE INTERFACE DETECTION")
    print("="*60)

    if 'spatial_connectivities' in adata.obsp:
        spatial_graph = adata.obsp['spatial_connectivities']

        is_luminal = adata.obs['cell_type'].str.contains('Lumm')
        is_plasma = adata.obs['cell_type'] == 'plasma_IgG'

        interface_spots = []

        for spot_idx in range(adata.n_obs):
            neighbors = spatial_graph[spot_idx].nonzero()[1]

            if len(neighbors) == 0:
                interface_spots.append(False)
                continue

            has_luminal_neighbor = is_luminal.iloc[neighbors].any()
            has_plasma_neighbor = is_plasma.iloc[neighbors].any()

            is_interface = has_luminal_neighbor and has_plasma_neighbor
            interface_spots.append(is_interface)

        adata.obs['tumor_immune_interface'] = interface_spots

        n_interface = sum(interface_spots)
        pct_interface = n_interface / adata.n_obs * 100

        print(f"\nInterface spots: {n_interface} ({pct_interface:.1f}%)")

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        sc.pl.spatial(
            adata,
            library_id=library_id,
            color='cell_type',
            ax=axes[0],
            title='Cell Types',
            size=1.3,
            show=False
        )

        sc.pl.spatial(
            adata,
            library_id=library_id,
            color='tumor_immune_interface',
            ax=axes[1],
            title=f'Tumor-Immune Interface ({n_interface} spots)',
            size=1.5,
            palette=['lightgray', 'red'],
            show=False
        )

        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / 'tumor_immune_interface.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {OUTPUT_DIR / 'tumor_immune_interface.png'}")
        plt.close()

    print("\n" + "="*60)
    print("CONFIDENCE ANALYSIS")
    print("="*60)

    threshold = 0.5
    low_conf = adata.obs['cell_type_confidence'] < threshold

    print(f"Low confidence spots (<{threshold}): {low_conf.sum()} ({low_conf.sum()/adata.n_obs*100:.1f}%)")

    conf_by_type = adata.obs.groupby('cell_type')['cell_type_confidence'].mean()
    print("\nMean confidence by cell type:")
    for ct, conf in conf_by_type.items():
        print(f"  {ct}: {conf:.3f}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='cell_type_confidence',
        ax=axes[0, 0],
        title='Prediction Confidence',
        size=1.3,
        cmap='viridis',
        show=False
    )

    adata.obs['low_confidence'] = low_conf.astype(str)
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='low_confidence',
        ax=axes[0, 1],
        title=f'Low Confidence Regions (<{threshold})',
        size=1.3,
        palette=['green', 'orange'],
        show=False
    )

    axes[1, 0].hist(adata.obs['cell_type_confidence'], bins=50, color='skyblue', edgecolor='black')
    axes[1, 0].axvline(threshold, color='red', linestyle='--', label=f'Threshold={threshold}')
    axes[1, 0].set_xlabel('Confidence Score')
    axes[1, 0].set_ylabel('Number of Spots')
    axes[1, 0].set_title('Confidence Distribution')
    axes[1, 0].legend()

    conf_by_type.plot(kind='barh', ax=axes[1, 1], color='steelblue')
    axes[1, 1].set_xlabel('Mean Confidence')
    axes[1, 1].set_title('Average Confidence by Cell Type')
    axes[1, 1].axvline(threshold, color='red', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'confidence_analysis.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR / 'confidence_analysis.png'}")
    plt.close()

    print("\n" + "="*60)
    print("DIFFERENTIAL EXPRESSION ANALYSIS")
    print("="*60)

    print("Running Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', n_genes=50)

    print("\nTop 10 marker genes per cell type:")
    for cell_type in adata.obs['cell_type'].cat.categories:
        print(f"\n{cell_type}:")
        genes = sc.get.rank_genes_groups_df(adata, group=cell_type).head(10)
        for idx, row in genes.iterrows():
            print(f"  {row['names']:15s} score={row['scores']:6.2f} padj={row['pvals_adj']:.2e}")

    print("\n" + "="*60)
    print("SAVING ENHANCED RESULTS")
    print("="*60)

    import json

    deg_summary = {}
    for cell_type in adata.obs['cell_type'].cat.categories:
        genes_df = sc.get.rank_genes_groups_df(adata, group=cell_type).head(20)
        deg_summary[cell_type] = genes_df['names'].tolist()

    enhanced_summary = {
        "cell_type_stats": {
            "total_spots": int(adata.n_obs),
            "n_cell_types": int(adata.obs['cell_type'].nunique()),
            "median_confidence": float(adata.obs['cell_type_confidence'].median()),
            "low_confidence_spots": int(low_conf.sum()),
            "low_confidence_pct": float(low_conf.sum() / adata.n_obs * 100)
        },
        "tumor_immune_interface": {
            "n_interface_spots": int(adata.obs.get('tumor_immune_interface', pd.Series([False]*adata.n_obs)).sum()),
            "interface_pct": float(adata.obs.get('tumor_immune_interface', pd.Series([False]*adata.n_obs)).sum() / adata.n_obs * 100)
        },
        "top_markers_per_celltype": deg_summary,
        "cell_type_composition": {k: int(v) for k, v in adata.obs['cell_type'].value_counts().items()},
        "confidence_by_celltype": {k: float(v) for k, v in adata.obs.groupby('cell_type')['cell_type_confidence'].mean().items()}
    }

    json_path = OUTPUT_DIR / "cell_type_enhanced_summary.json"
    with open(json_path, 'w') as f:
        json.dump(enhanced_summary, f, indent=2)
    print(f"Saved: {json_path}")

    h5ad_updated = OUTPUT_DIR / "annotated_visium_enhanced.h5ad"
    adata.write(h5ad_updated)
    print(f"Saved: {h5ad_updated}")
    print(f"Size: {h5ad_updated.stat().st_size / (1024**2):.1f} MB")

    print("\n" + "="*60)
    print("COMPLETE!")
    print("="*60)
    print("\nKey Findings:")
    print(f"  - {adata.obs['cell_type'].nunique()} cell types identified")
    if 'tumor_immune_interface' in adata.obs.columns:
        print(f"  - {adata.obs['tumor_immune_interface'].sum()} tumor-immune interface spots ({adata.obs['tumor_immune_interface'].sum()/adata.n_obs*100:.1f}%)")
    print(f"  - {low_conf.sum()} low confidence spots ({low_conf.sum()/adata.n_obs*100:.1f}%)")
    print(f"  - Marker genes validated for {len([v for v in available_markers.values() if v])} categories")

#!/usr/bin/env python
"""Run CellTypist annotation on processed Visium data."""

if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')

    import scanpy as sc
    import celltypist
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path

    print(f"scanpy version: {sc.__version__}")
    print(f"celltypist version: {celltypist.__version__}")

    SEED = 42
    np.random.seed(SEED)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

    PROJECT_ROOT = Path("/Users/sriharshameghadri/randomAIProjects/kaggle/medGemma")
    OUTPUT_DIR = PROJECT_ROOT / "outputs"
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("\n" + "="*60)
    print("LOADING PROCESSED DATA")
    print("="*60)

    h5ad_path = OUTPUT_DIR / "processed_visium.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Processed data not found at {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded: {adata.shape}")
    print(f"Clusters: {adata.obs['leiden'].nunique()}")

    print("\n" + "="*60)
    print("DOWNLOADING CELLTYPIST MODEL")
    print("="*60)

    model_name = 'Cells_Adult_Breast.pkl'
    print(f"Model: {model_name}")
    celltypist.models.download_models(model=model_name, force_update=False)

    model = celltypist.models.Model.load(model=model_name)
    print(f"Cell types: {len(model.cell_types)}")
    print(f"Features: {len(model.features)}")

    print("\n" + "="*60)
    print("PREPARING DATA")
    print("="*60)

    if adata.raw is None:
        raise ValueError("adata.raw is missing. Re-run 01_spatial_analysis.ipynb")

    adata_for_celltypist = adata.raw.to_adata()
    print(f"Prepared data: {adata_for_celltypist.shape}")

    print("\n" + "="*60)
    print("RUNNING CELLTYPIST ANNOTATION")
    print("="*60)
    print("This may take 30-60 seconds...")

    predictions = celltypist.annotate(
        adata_for_celltypist,
        model=model_name,
        majority_voting=True
    )

    print("Prediction complete!")

    adata_annotated = predictions.to_adata()

    print(f"\nAnnotation columns:")
    print(f"  predicted_labels: {adata_annotated.obs['predicted_labels'].nunique()} types")
    print(f"  majority_voting: {adata_annotated.obs['majority_voting'].nunique()} types")
    print(f"  conf_score: {adata_annotated.obs['conf_score'].min():.2f} - {adata_annotated.obs['conf_score'].max():.2f}")

    print("\n" + "="*60)
    print("TRANSFERRING ANNOTATIONS")
    print("="*60)

    adata.obs['cell_type'] = adata_annotated.obs['majority_voting'].values
    adata.obs['cell_type_raw'] = adata_annotated.obs['predicted_labels'].values
    adata.obs['cell_type_confidence'] = adata_annotated.obs['conf_score'].values

    print(f"Cell types identified: {adata.obs['cell_type'].nunique()}")

    print("\n" + "="*60)
    print("ANALYZING CELL TYPE COMPOSITION")
    print("="*60)

    cell_type_counts = adata.obs['cell_type'].value_counts()
    print("\nTop 15 cell types:")
    print(cell_type_counts.head(15))

    print(f"\nMedian confidence: {adata.obs['cell_type_confidence'].median():.2f}")
    print(f"Mean confidence: {adata.obs['cell_type_confidence'].mean():.2f}")

    print("\n" + "="*60)
    print("CREATING VISUALIZATIONS")
    print("="*60)

    library_id = list(adata.uns['spatial'].keys())[0]

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='cell_type',
        ax=axes[0, 0],
        title='Cell Types (Majority Voting)',
        size=1.3,
        show=False,
        legend_loc=None
    )

    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='cell_type_confidence',
        ax=axes[0, 1],
        title='Prediction Confidence',
        size=1.3,
        cmap='viridis',
        show=False
    )

    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='leiden',
        ax=axes[1, 0],
        title='Leiden Clusters',
        size=1.3,
        show=False
    )

    most_common = cell_type_counts.index[0]
    adata.obs['is_most_common'] = (adata.obs['cell_type'] == most_common).astype(str)
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color='is_most_common',
        ax=axes[1, 1],
        title=f'Highlighted: {most_common}',
        size=1.3,
        show=False,
        palette=['lightgray', 'red']
    )

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'spatial_cell_types.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

    print("\n" + "="*60)
    print("CLUSTER VS CELL TYPE ANALYSIS")
    print("="*60)

    crosstab = pd.crosstab(adata.obs['leiden'], adata.obs['cell_type'], normalize='index') * 100

    print("\nDominant cell type per cluster:")
    for cluster in crosstab.index:
        dominant_type = crosstab.loc[cluster].idxmax()
        percentage = crosstab.loc[cluster, dominant_type]
        print(f"  Cluster {cluster}: {dominant_type} ({percentage:.1f}%)")

    print("\n" + "="*60)
    print("SAVING ANNOTATED DATA")
    print("="*60)

    annotated_h5ad = OUTPUT_DIR / "annotated_visium.h5ad"
    adata.write(annotated_h5ad)
    print(f"Saved: {annotated_h5ad}")
    print(f"Size: {annotated_h5ad.stat().st_size / (1024**2):.1f} MB")

    import json

    cell_type_summary = {
        "total_spots": int(adata.n_obs),
        "n_cell_types": int(adata.obs['cell_type'].nunique()),
        "median_confidence": float(adata.obs['cell_type_confidence'].median()),
        "mean_confidence": float(adata.obs['cell_type_confidence'].mean()),
        "cell_type_counts": {k: int(v) for k, v in cell_type_counts.head(20).items()},
        "cell_type_proportions": {k: float(v) for k, v in (cell_type_counts.head(20) / adata.n_obs * 100).items()},
        "dominant_types_per_cluster": {
            str(cluster): crosstab.loc[cluster].idxmax()
            for cluster in crosstab.index
        }
    }

    json_path = OUTPUT_DIR / "cell_type_summary.json"
    with open(json_path, 'w') as f:
        json.dump(cell_type_summary, f, indent=2)
    print(f"Saved: {json_path}")

    print("\n" + "="*60)
    print("COMPLETE!")
    print("="*60)
    print(f"\nCheck outputs:")
    print(f"  - {output_path}")
    print(f"  - {annotated_h5ad}")
    print(f"  - {json_path}")

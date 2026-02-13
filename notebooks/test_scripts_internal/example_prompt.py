#!/usr/bin/env python3
"""
Example Clinical Prompt
View the prompt template that will be sent to MedGemma without running the model.

Usage:
    python example_prompt.py
"""

import json


def load_example_data():
    """Load actual data from outputs directory."""
    with open('../outputs/spatial_features.json', 'r') as f:
        spatial_data = json.load(f)

    with open('../outputs/cell_type_enhanced_summary.json', 'r') as f:
        celltype_data = json.load(f)

    return spatial_data, celltype_data


def create_clinical_prompt(spatial_data, celltype_data):
    """Generate clinical pathology prompt (same as run_medgemma.py)."""

    total_spots = celltype_data['cell_type_stats']['total_spots']
    composition = celltype_data['cell_type_composition']

    luminal_scgb = composition.get('LummHR-SCGB', 0)
    luminal_major = composition.get('LummHR-major', 0)
    plasma = composition.get('plasma_IgG', 0)

    luminal_pct = ((luminal_scgb + luminal_major) / total_spots) * 100
    plasma_pct = (plasma / total_spots) * 100

    interface_pct = celltype_data['tumor_immune_interface']['interface_pct']

    top_morans = list(spatial_data['spatial_statistics']['morans_i']['top_genes'].keys())[:5]

    top_luminal_markers = celltype_data['top_markers_per_celltype'].get('LummHR-SCGB', [])[:5]
    top_plasma_markers = celltype_data['top_markers_per_celltype'].get('plasma_IgG', [])[:5]

    prompt = f"""You are a board-certified pathologist reviewing spatial transcriptomics data from a breast cancer biopsy specimen. Based on the following molecular and spatial findings, generate a concise clinical pathology report (approximately 200 words).

SPECIMEN DATA:
- Total analyzed tissue spots: {total_spots}
- Spatial resolution: 10x Genomics Visium (55μm spots)

CELL TYPE COMPOSITION:
- Luminal epithelial cells (HR+): {luminal_pct:.1f}% of tissue
  * Primary subtype: LummHR-SCGB ({luminal_scgb} spots)
  * Secondary subtype: LummHR-major ({luminal_major} spots)
  * Key markers: {', '.join(top_luminal_markers)}

- Plasma cells (IgG+): {plasma_pct:.1f}% of tissue
  * Count: {plasma} spots
  * Key markers: {', '.join(top_plasma_markers)}

SPATIAL ORGANIZATION:
- Tumor-immune interface: {interface_pct:.1f}% of tissue shows direct contact between tumor and immune cells
- Spatial clusters identified: {spatial_data['dataset_info']['n_clusters']} distinct tissue regions
- Spatially variable genes (Moran's I > 0.1): {spatial_data['spatial_statistics']['morans_i']['n_significant_genes']} genes
  * Top spatially clustered genes: {', '.join(top_morans)}

MOLECULAR FEATURES:
- Quality metrics: {spatial_data['qc_metrics']['mean_genes_per_spot']:.0f} genes/spot, {spatial_data['qc_metrics']['mean_mt_percent']:.1f}% mitochondrial

Generate a clinical report including:
1. Diagnosis and tumor classification
2. Immune microenvironment assessment
3. Spatial organization patterns
4. Clinical implications for treatment
5. Brief statement on analytical limitations

Report should be professional, concise, and suitable for inclusion in patient medical records.

CLINICAL PATHOLOGY REPORT:"""

    return prompt


def main():
    print("=" * 80)
    print("Example Clinical Prompt for MedGemma")
    print("=" * 80)

    print("\nThis is the prompt that will be sent to MedGemma-4b-it.")
    print("No model loading required - just viewing the template.\n")

    try:
        spatial_features, celltype_summary = load_example_data()

        print("Data Summary:")
        print("-" * 80)
        print(f"Total spots: {spatial_features['dataset_info']['total_spots']}")
        print(f"Spatial clusters: {spatial_features['dataset_info']['n_clusters']}")
        print(f"Spatially variable genes: {spatial_features['spatial_statistics']['morans_i']['n_significant_genes']}")

        total = celltype_summary['cell_type_stats']['total_spots']
        print("\nCell Type Composition:")
        for celltype, count in celltype_summary['cell_type_composition'].items():
            pct = (count / total) * 100
            print(f"  {celltype}: {count} ({pct:.1f}%)")

        print(f"\nTumor-Immune Interface: {celltype_summary['tumor_immune_interface']['interface_pct']:.1f}%")

        prompt = create_clinical_prompt(spatial_features, celltype_summary)

        print("\n" + "=" * 80)
        print("GENERATED PROMPT")
        print("=" * 80)
        print(prompt)
        print("=" * 80)

        print(f"\nPrompt length: {len(prompt)} characters")
        print(f"Estimated tokens: ~{len(prompt.split())} tokens")

        print("\nTo generate actual report:")
        print("  jupyter notebook 04_medgemma_reports.ipynb")
        print("  python run_medgemma.py")

    except FileNotFoundError as e:
        print(f"\n✗ Data files not found: {e}")
        print("\nRun spatial analysis first to generate:")
        print("  - outputs/spatial_features.json")
        print("  - outputs/cell_type_enhanced_summary.json")


if __name__ == "__main__":
    main()

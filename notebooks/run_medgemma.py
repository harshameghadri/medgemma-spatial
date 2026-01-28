#!/usr/bin/env python3
"""
MedGemma Clinical Report Generator
Executable script for generating pathology reports from spatial transcriptomics data.

Usage:
    python run_medgemma.py --spatial spatial_features.json --celltype cell_type_summary.json
    python run_medgemma.py  # Uses default paths in outputs/
"""

import json
import torch
import argparse
from pathlib import Path
from datetime import datetime
from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig
import warnings
warnings.filterwarnings('ignore')


def load_analysis_data(spatial_path, celltype_path):
    """Load spatial features and cell type annotations."""
    with open(spatial_path, 'r') as f:
        spatial_data = json.load(f)

    with open(celltype_path, 'r') as f:
        celltype_data = json.load(f)

    return spatial_data, celltype_data


def load_medgemma_model(model_id="google/medgemma-4b-it"):
    """Load MedGemma-4b-it with 4-bit quantization."""
    quantization_config = BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_compute_dtype=torch.float16,
        bnb_4bit_quant_type="nf4",
        bnb_4bit_use_double_quant=True
    )

    print(f"Loading {model_id} with 4-bit quantization...")

    tokenizer = AutoTokenizer.from_pretrained(model_id)

    model = AutoModelForCausalLM.from_pretrained(
        model_id,
        quantization_config=quantization_config,
        device_map="auto",
        trust_remote_code=True
    )

    print("Model loaded successfully")
    return tokenizer, model


def create_clinical_prompt(spatial_data, celltype_data):
    """Generate clinical pathology prompt from spatial analysis."""

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


def generate_report(tokenizer, model, prompt, max_length=400):
    """Generate clinical report using MedGemma."""

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    print("Generating clinical report...")

    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_length,
            temperature=0.7,
            top_p=0.9,
            do_sample=True,
            pad_token_id=tokenizer.eos_token_id
        )

    full_output = tokenizer.decode(outputs[0], skip_special_tokens=True)

    report = full_output.split("CLINICAL PATHOLOGY REPORT:")[-1].strip()

    return report


def save_report(report, spatial_data, celltype_data, output_dir):
    """Save clinical report as TXT and structured JSON."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    txt_path = output_dir / f"clinical_report_{timestamp}.txt"
    with open(txt_path, 'w') as f:
        f.write("SPATIAL TRANSCRIPTOMICS CLINICAL PATHOLOGY REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Analysis Method: 10x Genomics Visium + MedGemma-4b-it\n")
        f.write(f"Specimen Type: Breast Cancer Biopsy\n\n")
        f.write("REPORT:\n")
        f.write("-" * 80 + "\n")
        f.write(report)
        f.write("\n" + "=" * 80 + "\n")

    report_data = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "model": "google/medgemma-4b-it",
            "quantization": "4-bit (NF4)",
            "word_count": len(report.split())
        },
        "clinical_report": report,
        "input_data_summary": {
            "total_spots": spatial_data['dataset_info']['total_spots'],
            "n_clusters": spatial_data['dataset_info']['n_clusters'],
            "n_spatially_variable_genes": spatial_data['spatial_statistics']['morans_i']['n_significant_genes'],
            "cell_type_composition": celltype_data['cell_type_composition'],
            "tumor_immune_interface_pct": celltype_data['tumor_immune_interface']['interface_pct']
        }
    }

    json_path = output_dir / f"clinical_report_{timestamp}.json"
    with open(json_path, 'w') as f:
        json.dump(report_data, f, indent=2)

    return txt_path, json_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate clinical pathology reports from spatial transcriptomics data"
    )
    parser.add_argument(
        '--spatial',
        type=str,
        default='../outputs/spatial_features.json',
        help='Path to spatial features JSON file'
    )
    parser.add_argument(
        '--celltype',
        type=str,
        default='../outputs/cell_type_enhanced_summary.json',
        help='Path to cell type summary JSON file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='../outputs',
        help='Output directory for reports'
    )
    parser.add_argument(
        '--model',
        type=str,
        default='google/medgemma-4b-it',
        help='HuggingFace model ID'
    )

    args = parser.parse_args()

    print("=" * 80)
    print("MedGemma Clinical Report Generator")
    print("=" * 80)

    print(f"\nPyTorch version: {torch.__version__}")
    print(f"MPS available: {torch.backends.mps.is_available()}")
    print(f"Device: {'mps' if torch.backends.mps.is_available() else 'cpu'}")

    print("\n1. Loading analysis data...")
    spatial_features, celltype_summary = load_analysis_data(args.spatial, args.celltype)

    print(f"   Total spots: {spatial_features['dataset_info']['total_spots']}")
    print(f"   Spatial clusters: {spatial_features['dataset_info']['n_clusters']}")

    total = celltype_summary['cell_type_stats']['total_spots']
    for celltype, count in celltype_summary['cell_type_composition'].items():
        pct = (count / total) * 100
        print(f"   {celltype}: {count} ({pct:.1f}%)")

    print("\n2. Loading MedGemma model...")
    tokenizer, model = load_medgemma_model(args.model)

    print("\n3. Creating clinical prompt...")
    prompt = create_clinical_prompt(spatial_features, celltype_summary)

    print("\n4. Generating clinical report...")
    clinical_report = generate_report(tokenizer, model, prompt)

    print("\n" + "=" * 80)
    print("GENERATED CLINICAL PATHOLOGY REPORT")
    print("=" * 80)
    print(clinical_report)
    print("=" * 80)
    print(f"\nWord count: {len(clinical_report.split())}")

    print("\n5. Saving reports...")
    txt_file, json_file = save_report(
        clinical_report,
        spatial_features,
        celltype_summary,
        args.output
    )

    print(f"\nReports saved:")
    print(f"  TXT: {txt_file}")
    print(f"  JSON: {json_file}")

    try:
        import psutil
        process = psutil.Process()
        memory_gb = process.memory_info().rss / 1024**3
        print(f"\nMemory usage: {memory_gb:.2f} GB")
        print(f"Status: {'✓ Within 64GB limit' if memory_gb < 32 else '⚠ High memory usage'}")
    except ImportError:
        print("\nInstall psutil for memory monitoring: pip install psutil")

    print("\n✓ Complete")


if __name__ == "__main__":
    main()

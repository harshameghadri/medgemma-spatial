#!/usr/bin/env python3
"""
MedGemma Spatial Transcriptomics - Streamlit Web Application

Production-ready spatial analysis and report generation interface.
"""

import streamlit as st
import sys
from pathlib import Path
import tempfile
import json
import time
from datetime import datetime

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import scanpy as sc
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    from src.streamlit_adapter import (
        annotate_spatial_regions,
        calculate_spatial_heterogeneity
    )
    from src.report_generation.medgemma_multimodal import (
        generate_multimodal_report,
        generate_textonly_fallback
    )
    from notebooks.medgemma_report_generator import create_anti_parroting_prompt
except ImportError as e:
    st.error(f"Missing dependencies: {e}")
    st.stop()


# Page config
st.set_page_config(
    page_title="MedGemma Spatial Transcriptomics",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)


def create_header():
    """Create app header with branding."""
    st.title("🔬 MedGemma Spatial Transcriptomics")
    st.markdown("""
    **AI-powered spatial transcriptomics analysis and pathology report generation**

    Upload Visium HD spatial transcriptomics data to generate:
    - Cell type annotations with marker-based identification
    - Spatial heterogeneity analysis
    - Clinical pathology reports (with optional H&E image integration)
    """)
    st.divider()


def create_sidebar():
    """Create sidebar with settings and info."""
    with st.sidebar:
        st.header("⚙️ Settings")

        use_multimodal = st.checkbox(
            "Use Multimodal (H&E + Text)",
            value=False,
            help="Requires MedGemma 1.5. Falls back to text-only if unavailable."
        )

        use_markers = st.checkbox(
            "Marker-based Annotation",
            value=True,
            help="Use PanglaoDB markers for cell type identification"
        )

        resolution = st.slider(
            "Leiden Resolution",
            min_value=0.1,
            max_value=2.0,
            value=0.5,
            step=0.1,
            help="Higher = more clusters"
        )

        st.divider()

        st.header("ℹ️ About")
        st.markdown("""
        **Pipeline:**
        1. Quality control & normalization
        2. Spatial clustering (Leiden)
        3. Marker-based cell type annotation
        4. Spatial heterogeneity analysis
        5. MedGemma report generation

        **Models:**
        - Scanpy/Squidpy for spatial analysis
        - PanglaoDB markers for annotation
        - MedGemma 1.5 for report generation

        **Data Privacy:**
        - Tissue-blind prompts (no identifying info)
        - Local processing only
        - No data stored on servers
        """)

        st.divider()

        st.markdown("""
        **GitHub:** [medgemma-spatial](https://github.com/harshameghadri/medgemma-spatial)

        Built with Claude Code
        """)

    return use_multimodal, use_markers, resolution


def load_data(uploaded_file):
    """Load H5AD file from upload."""
    with st.spinner("Loading spatial data..."):
        try:
            # Save to temp file
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp_path = tmp.name

            # Load with scanpy
            adata = sc.read_h5ad(tmp_path)

            st.success(f"✓ Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

            # Clean up
            Path(tmp_path).unlink()

            return adata

        except Exception as e:
            st.error(f"Failed to load file: {e}")
            return None


def run_spatial_analysis(adata, use_markers, resolution):
    """Run spatial analysis pipeline."""
    progress_bar = st.progress(0)
    status_text = st.empty()

    try:
        # Step 1: Annotation
        status_text.text("Step 1/3: Annotating spatial regions...")
        progress_bar.progress(10)

        adata, annot_metrics = annotate_spatial_regions(
            adata,
            resolution=resolution,
            use_markers=use_markers,
            tissue="Unknown"  # Tissue-blind
        )

        progress_bar.progress(40)

        # Step 2: Spatial heterogeneity
        status_text.text("Step 2/3: Calculating spatial heterogeneity...")

        spatial_metrics = calculate_spatial_heterogeneity(adata)

        progress_bar.progress(70)

        # Step 3: Extract features
        status_text.text("Step 3/3: Extracting features...")

        features = {
            'annotation': annot_metrics,
            'spatial_heterogeneity': spatial_metrics,
            'uncertainty': {
                'mean_prediction_entropy': float(adata.obs.get('prediction_entropy', [1.0]).mean())
                    if 'prediction_entropy' in adata.obs else 1.0
            }
        }

        progress_bar.progress(100)
        status_text.text("✓ Analysis complete!")

        return adata, features, True

    except Exception as e:
        st.error(f"Analysis failed: {e}")
        import traceback
        st.code(traceback.format_exc())
        return adata, None, False


def visualize_spatial_data(adata):
    """Create interactive spatial visualizations."""
    st.header("📊 Spatial Visualization")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Cell Type Distribution")

        if 'cell_type' in adata.obs:
            cell_counts = adata.obs['cell_type'].value_counts()

            fig = px.bar(
                x=cell_counts.index,
                y=cell_counts.values,
                labels={'x': 'Cell Type', 'y': 'Count'},
                title='Cell Type Composition'
            )
            fig.update_layout(xaxis_tickangle=-45, height=400)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No cell type annotations available")

    with col2:
        st.subheader("Cluster Distribution")

        if 'spatial_region' in adata.obs:
            cluster_counts = adata.obs['spatial_region'].value_counts()

            fig = px.pie(
                values=cluster_counts.values,
                names=cluster_counts.index,
                title='Spatial Clusters'
            )
            fig.update_traces(textposition='inside', textinfo='percent+label')
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No clustering available")

    # Metrics
    st.subheader("📈 Key Metrics")

    metrics_col1, metrics_col2, metrics_col3, metrics_col4 = st.columns(4)

    with metrics_col1:
        st.metric("Total Spots", f"{adata.n_obs:,}")

    with metrics_col2:
        st.metric("Genes", f"{adata.n_vars:,}")

    with metrics_col3:
        if 'cell_type' in adata.obs:
            n_types = len(adata.obs['cell_type'].unique())
            st.metric("Cell Types", n_types)
        else:
            st.metric("Cell Types", "N/A")

    with metrics_col4:
        if 'spatial_region' in adata.obs:
            n_clusters = len(adata.obs['spatial_region'].unique())
            st.metric("Clusters", n_clusters)
        else:
            st.metric("Clusters", "N/A")


def generate_report(adata, features, use_multimodal):
    """Generate MedGemma report."""
    st.header("📝 Pathology Report")

    with st.spinner("Generating report with MedGemma..."):
        try:
            start_time = time.time()

            if use_multimodal:
                report, metadata = generate_multimodal_report(
                    adata,
                    features,
                    device='cpu'  # Use CPU for demo
                )
            else:
                # Text-only using prompt
                prompt = create_anti_parroting_prompt(features)
                report = "[DEMO MODE] Report generation requires MedGemma model.\n\n"
                report += f"Prompt created successfully ({len(prompt)} chars).\n\n"
                report += "In production, this would generate a 200-word clinical pathology report "
                report += "correlating spatial patterns with tissue morphology."

                metadata = {
                    'mode': 'demo',
                    'prompt_length': len(prompt)
                }

            elapsed = time.time() - start_time

            st.success(f"✓ Report generated in {elapsed:.1f}s")

            # Display report
            st.markdown("### Generated Report")
            st.markdown(report)

            # Metadata
            with st.expander("📋 Report Metadata"):
                st.json(metadata)

            # Download button
            report_text = f"""MEDGEMMA SPATIAL TRANSCRIPTOMICS REPORT
{'='*80}

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Mode: {metadata.get('mode', 'unknown')}

REPORT:
{'-'*80}
{report}
{'='*80}

FEATURES:
{json.dumps(features, indent=2)}
"""

            st.download_button(
                label="📥 Download Report (TXT)",
                data=report_text,
                file_name=f"medgemma_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain"
            )

            return True

        except Exception as e:
            st.error(f"Report generation failed: {e}")
            import traceback
            st.code(traceback.format_exc())
            return False


def main():
    """Main application flow."""
    create_header()

    # Sidebar settings
    use_multimodal, use_markers, resolution = create_sidebar()

    # File upload
    st.header("📁 Upload Data")

    uploaded_file = st.file_uploader(
        "Upload Visium HD H5AD file",
        type=['h5ad'],
        help="Spatial transcriptomics data in AnnData H5AD format"
    )

    if uploaded_file is None:
        st.info("👆 Upload an H5AD file to begin analysis")

        # Show example
        with st.expander("📖 Example Data"):
            st.markdown("""
            **Expected format:** AnnData H5AD file with:
            - `adata.X`: Gene expression matrix
            - `adata.obs`: Spot metadata
            - `adata.var`: Gene metadata
            - `adata.obsm['spatial']`: Spatial coordinates (optional)
            - `adata.uns['spatial']`: H&E images (optional, for multimodal)

            **Sample data:**
            - [10x Genomics Visium datasets](https://www.10xgenomics.com/datasets)
            - Public spatial transcriptomics repositories
            """)

        st.stop()

    # Load data
    adata = load_data(uploaded_file)

    if adata is None:
        st.stop()

    # Run analysis button
    if st.button("🚀 Run Analysis", type="primary"):

        # Analysis
        adata, features, success = run_spatial_analysis(adata, use_markers, resolution)

        if not success:
            st.error("Analysis failed. Check errors above.")
            st.stop()

        # Store in session state
        st.session_state['adata'] = adata
        st.session_state['features'] = features

    # Show results if available
    if 'adata' in st.session_state and 'features' in st.session_state:

        adata = st.session_state['adata']
        features = st.session_state['features']

        # Visualizations
        visualize_spatial_data(adata)

        st.divider()

        # Generate report
        generate_report(adata, features, use_multimodal)


if __name__ == "__main__":
    main()

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path

st.set_page_config(
    page_title="RNA-seq Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS (theme-safe)
st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 20px;
        border-radius: 8px 8px 0 0;
    }
</style>
""", unsafe_allow_html=True)

# Paths
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"

# Sidebar
st.sidebar.title("🧬 RNA-seq Dashboard")
st.sidebar.markdown("---")

# Data source
data_source = st.sidebar.radio(
    "Data source",
    ["Demo dataset", "Airway (Himes et al.)", "Upload my data"],
    help="Use the built-in demo, the classic airway dataset, or upload your own counts + metadata files.",
)

@st.cache_data(show_spinner=False)
def load_demo():
    counts = pd.read_csv(DATA_DIR / "demo_counts.csv", index_col=0)
    metadata = pd.read_csv(DATA_DIR / "demo_metadata.csv")
    return counts, metadata

@st.cache_data(show_spinner=False)
def load_airway():
    counts = pd.read_csv(DATA_DIR / "airway_counts.csv", index_col=0)
    metadata = pd.read_csv(DATA_DIR / "airway_metadata.csv")
    return counts, metadata

if data_source == "Demo dataset":
    counts, metadata = load_demo()
    st.sidebar.success("Demo loaded: 2,500 genes × 12 samples")
elif data_source == "Airway (Himes et al.)":
    counts, metadata = load_airway()
    st.sidebar.success(f"Airway loaded: {counts.shape[0]:,} genes × {counts.shape[1]} samples")
else:
    st.sidebar.info("Upload TSV/CSV files with genes as rows and samples as columns.")
    counts_file = st.sidebar.file_uploader("Counts matrix", type=["csv", "tsv", "txt"])
    meta_file = st.sidebar.file_uploader("Metadata", type=["csv", "tsv", "txt"])

    if counts_file and meta_file:
        sep = "\t" if counts_file.name.endswith((".tsv", ".txt")) else ","
        counts = pd.read_csv(counts_file, index_col=0, sep=sep)
        sep_meta = "\t" if meta_file.name.endswith((".tsv", ".txt")) else ","
        metadata = pd.read_csv(meta_file, sep=sep_meta)
    else:
        st.warning("Please upload both files or switch to Demo dataset.")
        st.stop()

# Validate
required_meta_cols = {"Sample", "Condition"}
if not required_meta_cols.issubset(metadata.columns):
    st.error(f"Metadata must contain columns: {required_meta_cols}")
    st.stop()

missing_samples = set(metadata["Sample"]) - set(counts.columns)
if missing_samples:
    st.error(f"Samples in metadata not found in counts: {missing_samples}")
    st.stop()

# Subset counts to metadata samples
counts = counts[metadata["Sample"].tolist()]

# Sidebar parameters
st.sidebar.markdown("---")
st.sidebar.subheader("Analysis parameters")

conditions = metadata["Condition"].unique().tolist()
if len(conditions) < 2:
    st.error("Need at least 2 conditions in metadata.")
    st.stop()

group_a = st.sidebar.selectbox("Reference group", conditions, index=0)
group_b = st.sidebar.selectbox("Test group", [c for c in conditions if c != group_a], index=0)

pval_thresh = st.sidebar.slider("P-value threshold", 0.001, 0.2, 0.05, 0.001, format="%.3f")
log2fc_thresh = st.sidebar.slider("|log2FC| threshold", 0.0, 4.0, 0.5, 0.1)
use_padj = st.sidebar.checkbox("Use BH-adjusted p-value (FDR)", value=True, help="For small sample sizes (n<5), BH correction can be very conservative. Uncheck to use nominal p-values for exploratory analysis.")

st.sidebar.markdown("---")
st.sidebar.caption("v1.0 · Built with Streamlit & Plotly")

# Main header
st.title("🧬 RNA-seq Analysis Dashboard")
st.caption("Differential expression, quality control & pathway exploration")

# Metrics
from utils.analysis import cpm_normalize, log2_cpm, differential_expression, run_pca, get_de_genes, mock_enrichment
from utils.plots import plot_pca, plot_volcano, plot_heatmap, plot_ma, plot_sample_boxplots, plot_enrichment_bar

log_cpm = log2_cpm(counts)
cpm = cpm_normalize(counts)

@st.cache_data(show_spinner=False)
def run_de(counts, metadata, group_a, group_b):
    return differential_expression(counts, metadata, group_col="Condition", group_a=group_a, group_b=group_b, min_counts=50)

with st.spinner("Running differential expression analysis..."):
    de_results = run_de(counts, metadata, group_a, group_b)
    de_results = get_de_genes(de_results, pval_thresh=pval_thresh, log2fc_thresh=log2fc_thresh, use_padj=use_padj)

sig_up = ((de_results["padj"] < pval_thresh) & (de_results["log2FC"] > log2fc_thresh)).sum()
sig_down = ((de_results["padj"] < pval_thresh) & (de_results["log2FC"] < -log2fc_thresh)).sum()

m1, m2, m3, m4 = st.columns(4)
m1.metric("Genes", f"{counts.shape[0]:,}")
m2.metric("Samples", counts.shape[1])
m3.metric("Up-regulated", sig_up, delta_color="inverse")
m4.metric("Down-regulated", sig_down, delta_color="inverse")

st.markdown("---")

# Tabs
tab_overview, tab_qc, tab_pca, tab_de, tab_viz, tab_enrich = st.tabs([
    "📋 Overview", "📊 QC", "🔬 PCA", "📈 DE Analysis", "🎨 Visualizations", "🧪 Enrichment"
])

# --- Overview ---
with tab_overview:
    c1, c2 = st.columns([1, 1])
    with c1:
        st.subheader("Expression matrix")
        st.dataframe(
            counts.head(200).style.background_gradient(axis=1, cmap="YlGnBu").format("{:.0f}"),
            use_container_width=True,
            height=350,
        )
        st.caption(f"Showing top 200 of {counts.shape[0]} genes. Full matrix: {counts.shape[0]} × {counts.shape[1]}")
    with c2:
        st.subheader("Metadata")
        st.dataframe(metadata, use_container_width=True, hide_index=True)

    st.subheader("Top differentially expressed genes")
    top_de = de_results[de_results["Regulation"] != "Not significant"].head(20)
    if len(top_de) == 0:
        st.info("No significant DE genes found with current thresholds. Try relaxing them.")
    else:
        st.dataframe(
            top_de[["Gene", "log2FC", "pvalue", "padj", "Regulation"]].style.apply(
                lambda x: ["background: #fee2e2" if v == "Up-regulated" else "background: #dbeafe" if v == "Down-regulated" else "" for v in x],
                subset=["Regulation"],
            ),
            use_container_width=True,
            hide_index=True,
        )

# --- QC ---
with tab_qc:
    st.subheader("Sample expression distribution")
    fig_qc = plot_sample_boxplots(counts, metadata)
    st.plotly_chart(fig_qc, use_container_width=True)

    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Library sizes (total counts)")
        lib_sizes = counts.sum(axis=0).reset_index()
        lib_sizes.columns = ["Sample", "Total counts"]
        lib_sizes = lib_sizes.merge(metadata[["Sample", "Condition"]], on="Sample")
        fig_lib = px.bar(lib_sizes, x="Sample", y="Total counts", color="Condition", template="simple_white")
        fig_lib.update_layout(xaxis_tickangle=-45, height=400)
        st.plotly_chart(fig_lib, use_container_width=True)

    with c2:
        st.subheader("Detected genes per sample")
        detected = (counts > 0).sum(axis=0).reset_index()
        detected.columns = ["Sample", "Genes detected"]
        detected = detected.merge(metadata[["Sample", "Condition"]], on="Sample")
        fig_det = px.bar(detected, x="Sample", y="Genes detected", color="Condition", template="simple_white")
        fig_det.update_layout(xaxis_tickangle=-45, height=400)
        st.plotly_chart(fig_det, use_container_width=True)

# --- PCA ---
with tab_pca:
    st.subheader("Principal Component Analysis")
    pca_df, variance = run_pca(counts, n_components=2, top_var_genes=500)
    fig_pca = plot_pca(pca_df, metadata, variance)
    st.plotly_chart(fig_pca, use_container_width=True)
    st.info(f"PC1 explains {variance[0]:.1f}% and PC2 explains {variance[1]:.1f}% of total variance.")

# --- DE Analysis ---
with tab_de:
    st.subheader("Differential expression results")

    c1, c2 = st.columns([3, 1])
    with c2:
        search_gene = st.text_input("Search gene", placeholder="e.g. TP53")
    with c1:
        st.write("")

    if search_gene:
        filtered = de_results[de_results["Gene"].str.contains(search_gene, case=False, na=False)]
    else:
        filtered = de_results

    st.dataframe(
        filtered[["Gene", "log2FC", "pvalue", "padj", "Regulation"]].style.apply(
            lambda x: ["background: #fee2e2" if v == "Up-regulated" else "background: #dbeafe" if v == "Down-regulated" else "" for v in x],
            subset=["Regulation"],
        ),
        use_container_width=True,
        hide_index=True,
    )

    csv = de_results.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 Download full DE table (CSV)",
        data=csv,
        file_name=f"DE_{group_b}_vs_{group_a}.csv",
        mime="text/csv",
    )

# --- Visualizations ---
with tab_viz:
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Volcano plot")
        fig_volc = plot_volcano(de_results, pval_thresh=pval_thresh, log2fc_thresh=log2fc_thresh)
        st.plotly_chart(fig_volc, use_container_width=True)
    with c2:
        st.subheader("MA plot")
        fig_ma = plot_ma(de_results, pval_thresh=pval_thresh, use_padj=use_padj)
        st.plotly_chart(fig_ma, use_container_width=True)

    st.subheader("Heatmap (top DE genes)")
    top_n = st.slider("Number of genes", 10, 100, 40, 5)
    fig_hm = plot_heatmap(counts, de_results, metadata, top_n=top_n)
    st.plotly_chart(fig_hm, use_container_width=True)

# --- Enrichment ---
with tab_enrich:
    st.subheader("Pathway enrichment (demo)")
    st.info("This is a simulated enrichment for demonstration. Replace with gseapy.enrichr() for real analysis.")

    de_gene_list = de_results[de_results["Regulation"] != "Not significant"]["Gene"].tolist()
    bg_genes = de_results["Gene"].tolist()

    if len(de_gene_list) < 5:
        st.warning("Not enough DE genes for enrichment. Lower the thresholds.")
    else:
        enrich_df = mock_enrichment(de_gene_list, bg_genes, top_n=10)
        enrich_df["-log10(P)"] = -np.log10(enrich_df["P-value"])

        c1, c2 = st.columns([2, 3])
        with c1:
            st.dataframe(enrich_df[["Term", "P-value", "Adjusted P-value"]], use_container_width=True, hide_index=True)
        with c2:
            fig_enrich = plot_enrichment_bar(enrich_df)
            st.plotly_chart(fig_enrich, use_container_width=True)

st.markdown("---")
st.caption("Built by Marwan El Saabi · Data is for demonstration purposes only.")

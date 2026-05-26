"""
RNA-seq Plotting Module
-----------------------
Interactive Plotly visualizations for RNA-seq data.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plot_pca(pca_df: pd.DataFrame, metadata: pd.DataFrame, variance: np.ndarray, color_col: str = "Condition") -> go.Figure:
    """Interactive PCA scatter plot."""
    df = pca_df.merge(metadata, on="Sample")
    fig = px.scatter(
        df,
        x="PC1",
        y="PC2",
        color=color_col,
        hover_data=["Sample", "Batch", "RIN"],
        title="Principal Component Analysis",
        template="simple_white",
    )
    fig.update_traces(marker=dict(size=14, line=dict(width=1, color="DarkSlateGrey")))
    fig.update_layout(
        xaxis_title=f"PC1 ({variance[0]:.1f}% variance)",
        yaxis_title=f"PC2 ({variance[1]:.1f}% variance)",
        legend_title_text=color_col,
        height=550,
    )
    return fig


def plot_volcano(de_results: pd.DataFrame, pval_thresh: float = 0.05, log2fc_thresh: float = 1.0) -> go.Figure:
    """Interactive volcano plot."""
    df = de_results.copy()
    df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=1e-300))

    color_map = {
        "Up-regulated": "#EF4444",
        "Down-regulated": "#3B82F6",
        "Not significant": "#9CA3AF",
    }

    fig = px.scatter(
        df,
        x="log2FC",
        y="-log10(padj)",
        color="Regulation",
        color_discrete_map=color_map,
        hover_data=["Gene", "pvalue", "padj"],
        title="Differential Expression Volcano Plot",
        template="simple_white",
        opacity=0.8,
    )
    fig.update_traces(marker=dict(size=8))
    fig.add_hline(y=-np.log10(pval_thresh), line_dash="dash", line_color="grey", opacity=0.5)
    fig.add_vline(x=log2fc_thresh, line_dash="dash", line_color="grey", opacity=0.5)
    fig.add_vline(x=-log2fc_thresh, line_dash="dash", line_color="grey", opacity=0.5)
    fig.update_layout(
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10 Adjusted P-value",
        height=550,
    )
    return fig


def plot_heatmap(counts: pd.DataFrame, de_results: pd.DataFrame, metadata: pd.DataFrame,
                 top_n: int = 40, group_col: str = "Condition") -> go.Figure:
    """Heatmap of top differentially expressed genes."""
    sig_genes = de_results[de_results["Regulation"] != "Not significant"]["Gene"].head(top_n).tolist()
    if len(sig_genes) < 5:
        sig_genes = de_results["Gene"].head(top_n).tolist()

    from utils.analysis import log2_cpm
    log_cpm = log2_cpm(counts)
    mat = log_cpm.loc[sig_genes]

    # Z-score per gene
    mat_z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)

    # Reorder columns by condition
    meta_sorted = metadata.sort_values(by=[group_col, "Sample"])
    mat_z = mat_z[meta_sorted["Sample"].tolist()]

    fig = px.imshow(
        mat_z.values,
        x=mat_z.columns,
        y=mat_z.index,
        color_continuous_scale="RdBu_r",
        aspect="auto",
        title=f"Z-score Heatmap (Top {len(sig_genes)} DE Genes)",
        template="simple_white",
    )
    fig.update_layout(height=max(400, len(sig_genes) * 18))
    return fig


def plot_ma(de_results: pd.DataFrame, pval_thresh: float = 0.05, use_padj: bool = True) -> go.Figure:
    """MA plot: mean expression vs log2 fold change."""
    df = de_results.copy()
    df["mean_expr"] = (df["mean_A"] + df["mean_B"]) / 2
    df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=1e-300))

    pcol = "padj" if use_padj else "pvalue"
    df["Significant"] = df[pcol] < pval_thresh
    df["Color"] = df["Significant"].map({True: "#10B981", False: "#9CA3AF"})

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["mean_expr"],
        y=df["log2FC"],
        mode="markers",
        marker=dict(color=df["Color"], size=7, opacity=0.7),
        text=df["Gene"],
        hovertemplate="<b>%{text}</b><br>Mean expr: %{x:.2f}<br>log2FC: %{y:.2f}<extra></extra>",
    ))
    fig.add_hline(y=0, line_dash="dash", line_color="grey", opacity=0.5)
    fig.update_layout(
        title="MA Plot",
        xaxis_title="Mean Expression (log2 CPM)",
        yaxis_title="log2 Fold Change",
        template="simple_white",
        height=550,
        showlegend=False,
    )
    return fig


def plot_sample_boxplots(counts: pd.DataFrame, metadata: pd.DataFrame, group_col: str = "Condition") -> go.Figure:
    """Boxplots of log2-CPM distributions per sample."""
    from utils.analysis import log2_cpm
    log_cpm = log2_cpm(counts)
    df = log_cpm.T.reset_index().melt(id_vars=["index"], var_name="Gene", value_name="log2CPM")
    df = df.rename(columns={"index": "Sample"})
    df = df.merge(metadata[["Sample", group_col]], on="Sample")

    fig = px.box(
        df,
        x="Sample",
        y="log2CPM",
        color=group_col,
        title="Sample Expression Distribution (log2 CPM)",
        template="simple_white",
        points=False,
    )
    fig.update_layout(xaxis_tickangle=-45, height=500)
    return fig


def plot_enrichment_bar(enrich_df: pd.DataFrame) -> go.Figure:
    """Bar plot of enrichment results."""
    df = enrich_df.sort_values("P-value", ascending=True)
    fig = px.bar(
        df,
        x="-log10(P)",
        y="Term",
        orientation="h",
        color="-log10(P)",
        color_continuous_scale="Teal",
        title="Pathway Enrichment Analysis (Mock)",
        template="simple_white",
        hover_data=["Overlap", "Adjusted P-value"],
    )
    fig.update_layout(height=400, yaxis=dict(categoryorder="total ascending"))
    return fig

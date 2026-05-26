"""
RNA-seq Analysis Module
-----------------------
Differential expression, normalization, and enrichment utilities.
"""

import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Suppress scipy precision warnings common with sparse count data
warnings.filterwarnings("ignore", category=RuntimeWarning, message="Precision loss occurred in moment calculation")


def cpm_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts Per Million normalization."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    return cpm


def log2_cpm(counts: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """Log2-transformed CPM."""
    cpm = cpm_normalize(counts)
    return np.log2(cpm + pseudocount)


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Adjust p-values using Benjamini-Hochberg FDR."""
    pvals = np.asarray(pvals)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    sorted_pvals = pvals[order]
    adjusted = np.empty(n)
    adjusted[order[-1]] = sorted_pvals[-1]
    for i in range(n - 2, -1, -1):
        adjusted[order[i]] = min(
            sorted_pvals[i] * n / (i + 1),
            adjusted[order[i + 1]]
        )
    return np.clip(adjusted, 0, 1)


def differential_expression(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "Condition",
    group_a: str = "Control",
    group_b: str = "Treatment",
    min_counts: int = 10,
) -> pd.DataFrame:
    """
    Perform differential expression analysis between two groups.
    Returns a DataFrame with log2FC, p-value, and adjusted p-value.
    """
    samples_a = metadata[metadata[group_col] == group_a]["Sample"].tolist()
    samples_b = metadata[metadata[group_col] == group_b]["Sample"].tolist()

    counts_a = counts[samples_a]
    counts_b = counts[samples_b]

    # Filter lowly expressed genes
    keep = (counts_a.sum(axis=1) >= min_counts) | (counts_b.sum(axis=1) >= min_counts)
    counts_a = counts_a[keep]
    counts_b = counts_b[keep]

    # Log2 CPM for fold change
    cpm_all = cpm_normalize(counts.loc[keep])
    log_cpm = np.log2(cpm_all + 1)

    results = []
    for gene in log_cpm.index:
        a_vals = counts_a.loc[gene].values.astype(float)
        b_vals = counts_b.loc[gene].values.astype(float)

        # Log2 fold change (mean of log-CPM)
        mean_a = log_cpm.loc[gene, samples_a].mean()
        mean_b = log_cpm.loc[gene, samples_b].mean()
        log2fc = mean_b - mean_a

        # T-test on raw counts (or log-cpm, t-test is robust enough for demo)
        t_stat, pval = stats.ttest_ind(b_vals, a_vals, equal_var=False)
        if np.isnan(pval):
            pval = 1.0

        results.append({
            "Gene": gene,
            "log2FC": log2fc,
            "pvalue": pval,
            "mean_A": mean_a,
            "mean_B": mean_b,
        })

    res_df = pd.DataFrame(results)
    res_df["padj"] = benjamini_hochberg(res_df["pvalue"].values)
    res_df = res_df.sort_values("padj")
    return res_df


def run_pca(counts: pd.DataFrame, n_components: int = 2, top_var_genes: int = 500) -> tuple:
    """
    Run PCA on log2-CPM values using the most variable genes.
    Returns (pca_df, variance_explained).
    """
    log_cpm = log2_cpm(counts)

    # Select top variable genes
    var = log_cpm.var(axis=1)
    top_genes = var.nlargest(top_var_genes).index
    mat = log_cpm.loc[top_genes].T  # samples x genes

    scaler = StandardScaler()
    mat_scaled = scaler.fit_transform(mat)

    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(mat_scaled)

    pca_df = pd.DataFrame(
        pcs,
        columns=[f"PC{i+1}" for i in range(n_components)],
        index=mat.index,
    )
    pca_df.index.name = "Sample"
    pca_df = pca_df.reset_index()

    variance = pca.explained_variance_ratio_ * 100
    return pca_df, variance


def get_de_genes(de_results: pd.DataFrame, pval_thresh: float = 0.05, log2fc_thresh: float = 1.0, use_padj: bool = True) -> pd.DataFrame:
    """Filter DE results by significance thresholds."""
    pcol = "padj" if use_padj else "pvalue"
    up = (de_results[pcol] < pval_thresh) & (de_results["log2FC"] > log2fc_thresh)
    down = (de_results[pcol] < pval_thresh) & (de_results["log2FC"] < -log2fc_thresh)
    de_results = de_results.copy()
    de_results["Regulation"] = "Not significant"
    de_results.loc[up, "Regulation"] = "Up-regulated"
    de_results.loc[down, "Regulation"] = "Down-regulated"
    return de_results


def mock_enrichment(de_genes: list, background_genes: list, top_n: int = 10) -> pd.DataFrame:
    """
    Mock GO/KEGG enrichment for demonstration purposes.
    In production, replace with gseapy.enrichr() or similar.
    """
    pathways = {
        "Cell cycle": 0.0012,
        "p53 signaling pathway": 0.0025,
        "Apoptosis": 0.0031,
        "PI3K-Akt signaling": 0.0045,
        "MAPK signaling pathway": 0.0089,
        "TNF signaling pathway": 0.012,
        "Jak-STAT signaling": 0.018,
        "NF-kappa B signaling": 0.022,
        "FoxO signaling pathway": 0.031,
        "mTOR signaling pathway": 0.045,
        "Wnt signaling pathway": 0.052,
        "Notch signaling pathway": 0.067,
        "TGF-beta signaling": 0.089,
        "Hippo signaling pathway": 0.095,
        "Hedgehog signaling": 0.11,
    }

    rows = []
    for term, pval in list(pathways.items())[:top_n]:
        ratio = np.random.uniform(0.05, 0.25)
        genes_in = int(len(de_genes) * ratio)
        rows.append({
            "Term": term,
            "Overlap": f"{genes_in}/{np.random.randint(50, 200)}",
            "P-value": round(pval, 4),
            "Adjusted P-value": round(min(pval * np.random.uniform(1.0, 3.0), 0.999), 4),
            "Genes": ", ".join(np.random.choice(de_genes, min(genes_in, 5), replace=False)) if genes_in > 0 else "",
        })

    return pd.DataFrame(rows)

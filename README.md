# 🧬 RNA-seq Analysis Dashboard

An interactive, professional-grade RNA-seq analysis dashboard built with **Streamlit** and **Plotly**. Designed for exploratory analysis of gene expression data, differential expression visualization, and pathway enrichment — all in the browser.

![Dashboard Preview](https://img.shields.io/badge/Built%20with-Streamlit-FF4B4B?logo=streamlit)
![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python)

---

## ✨ Features

| Module | Description |
|--------|-------------|
| **📋 Overview** | Expression matrix preview, metadata viewer, top DE genes table |
| **📊 QC** | Per-sample boxplots, library sizes, detected genes |
| **🔬 PCA** | 2D principal component analysis with condition coloring |
| **📈 DE Analysis** | T-test based differential expression with log2FC, p-values, and BH correction. Toggle between nominal p-value and FDR |
| **🎨 Visualizations** | Interactive volcano plot, MA plot, and z-score heatmap |
| **🧪 Enrichment** | Simulated pathway enrichment (replaceable with `gseapy`) |

### Built-in datasets

| Dataset | Samples | Genes | Description |
|---------|---------|-------|-------------|
| **Demo** | 12 | 2,500 | Synthetic data with strong DE signal for quick testing |
| **Airway** | 8 | 38,694 | Real RNA-seq data from Himes et al. (2014). Human airway smooth muscle cells treated with dexamethasone vs control |

---

## 🚀 Quick Start

### 1. Clone or navigate to the project

```bash
cd rnaseq-dashboard
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the app

```bash
streamlit run app.py
```

The dashboard will open at `http://localhost:8501`.

---

## 📁 Project Structure

```
rnaseq-dashboard/
├── app.py                  # Main Streamlit application
├── requirements.txt        # Python dependencies
├── README.md               # This file
├── data/
│   ├── demo_counts.csv     # Synthetic gene count matrix (2500 genes × 12 samples)
│   └── demo_metadata.csv   # Sample metadata (condition, batch, RIN)
└── utils/
    ├── analysis.py         # DE, normalization, PCA, enrichment logic
    └── plots.py            # Plotly visualization helpers
```

---

## 📊 Input Data Format

### Counts matrix
- Rows = genes (index)
- Columns = samples
- Values = raw read counts (integers)

| Gene | Sample_01 | Sample_02 | ... |
|------|-----------|-----------|-----|
| BRCA1| 1245      | 980       | ... |
| TP53 | 560       | 612       | ... |

### Metadata
- Must contain at least `Sample` and `Condition` columns
- Additional columns (e.g., `Batch`, `RIN`) are optional and used for QC/PCA

| Sample | Condition | Batch | RIN |
|--------|-----------|-------|-----|
| Sample_01 | Control | A | 8.5 |
| Sample_02 | Treatment | A | 8.2 |

---

## 🔬 Methods

- **Normalization**: Counts Per Million (CPM) + log2 transformation
- **Differential Expression**: Welch's t-test per gene, Benjamini-Hochberg FDR correction
- **PCA**: Standardized log2-CPM of top 500 most variable genes
- **Enrichment**: Mock pathway analysis (replace `mock_enrichment()` with `gseapy.enrichr()` for real data)

---

## 🛠️ Tech Stack

- **Frontend**: Streamlit
- **Visualization**: Plotly
- **Analysis**: NumPy, Pandas, SciPy, scikit-learn

---

## 📝 Author

**Marwan El Saabi** — MSc Bioinformatics student & Bioinformatician  
[Portfolio](https://marwansaabi.github.io) · [GitHub](https://github.com/marwansaabi) · [LinkedIn](https://www.linkedin.com/in/marwansaabi/)

---

## 📄 License

MIT License — feel free to use and modify for your own research or portfolio.

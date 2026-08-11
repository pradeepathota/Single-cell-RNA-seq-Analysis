# Single-Cell RNA Analysis: Melanoma Tumor Microenvironment

A Scanpy-based pipeline for analyzing melanoma single-cell RNA-seq data (GSE72056), including quality control, dimensionality reduction, clustering, and cell type annotation.

## Overview

This pipeline processes raw expression data and performs:
- **Quality Control**: Filters low-quality cells and sparse genes
- **Preprocessing**: Normalization and selection of highly variable genes
- **Dimensionality Reduction**: PCA → kNN graph → UMAP
- **Clustering**: Leiden graph-based clustering
- **Marker Detection**: Wilcoxon rank-sum test for differentially expressed genes
- **Cell Type Annotation**: Manual annotation based on canonical markers and crosstab validation

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python scRNA.py
```

Output files will be saved to the `results/` directory:
- `*_clusters.png` - UMAP colored by Leiden clusters
- `*_malignant_vs_normal.png` - UMAP colored by malignancy status
- `*_markers.png` - Top marker genes per cluster
- `*_dotplot.png` - Canonical marker expression across clusters
- `*_final_annotated.png` - Final annotated cell types
- `cluster_vs_published_annotation.csv` - Cluster validation table

## Data

- **Source**: GSE72056 - Melanoma single-cell RNA-seq dataset
- **File**: `GSE72056_melanoma_single_cell_revised_v2.txt`
- **Format**: Gene × Cell expression matrix (genes in rows, cells in columns)
- **Metadata**: First 3 rows contain tumor ID, malignancy status, and published cell type

## Requirements

- Python 3.8+
- See `requirements.txt` for package versions

## Key Parameters

- **Highly Variable Genes**: Top 2000 genes selected
- **PCA Components**: 30 components retained
- **kNN Neighbors**: 10 neighbors used for graph construction
- **Leiden Resolution**: 0.6 (controls cluster granularity)

## Notes

- Modify `cluster_to_label` mapping in the final section based on your marker gene analysis results
- Cell type annotation may vary with different clustering parameters
- Ensure the data file is in the correct path before running the script

## Author

Pradeep Chowdary

"""
Single-cell RNA analysis pipeline for melanoma TME (tumor microenvironment).
This script processes raw expression data, performs quality control, dimensionality reduction,
clustering, and cell type annotation using the Scanpy library.

Data source: GSE72056 - Melanoma single-cell RNA-seq dataset
"""

import pandas as pd
import numpy as np
import scanpy as sc

# Load the raw expression matrix from TSV file
# Format: Genes in rows (first column as gene names), cells in columns
raw = pd.read_csv('~/Projects/scRNA/GSE72056_melanoma_single_cell_revised_v2.txt', sep="\t")

# Inspect the loaded data
print("Loaded shape:", raw.shape)
print(raw.iloc[:5, :5])

# Extract metadata from first 3 rows
# Row 1: tumor ID - identifies which tumor sample each cell comes from
# Row 2: malignant flag - 0=unresolved, 1=non-malignant, 2=malignant cancer cell
# Row 3: cell-type flag - identifies immune/stromal cell types in non-malignant cells
meta = raw.iloc[0:3, 1:]
meta.index = ["tumor", "malignant", "non_malignant_cell_type"]
meta = meta.T  # Transpose so cells are rows, metadata fields are columns

# Extract expression matrix (all rows except the 3 metadata rows)
# Set gene names as index for easier downstream access
expr = raw.iloc[3:, :].set_index(raw.columns[0])
expr = expr.apply(pd.to_numeric)  # Ensure all values are numeric

# Create AnnData object (standard format for single-cell analysis in Scanpy)
# Scanpy expects cells as rows (observations) and genes as columns (variables)
# so we transpose the expression matrix
adata = sc.AnnData(X=expr.T.values, obs=meta, var=pd.DataFrame(index=expr.index))
adata.var_names_make_unique()  # Ensure all gene names are unique (rename duplicates if any)

print(adata)

# Map numeric malignancy codes to human-readable labels
malignant_map = {"0": "Unresolved", "1": "Non-malignant", "2": "Malignant"}
adata.obs["malignant"] = adata.obs["malignant"].astype(str).map(malignant_map)

# Map numeric cell-type codes to human-readable immune/stromal cell type labels
# (only applies to non-malignant cells)
celltype_map = {
    "0": "Unresolved", "1": "T cell", "2": "B cell", "3": "Macrophage",
    "4": "Endothelial", "5": "CAF", "6": "NK",
}
adata.obs["published_cell_type"] = adata.obs["non_malignant_cell_type"].astype(str).map(celltype_map)

# Quality control: Remove low-quality cells and genes with sparse expression
# Filter cells: Keep only cells with at least 200 detected genes (minimum coverage threshold)
# This removes doublets, dead cells, and low-quality captures
sc.pp.filter_cells(adata, min_genes=200)

# Filter genes: Keep only genes expressed in at least 3 cells
# This removes noise and non-biological genes with only background signal
sc.pp.filter_genes(adata, min_cells=3)
print(f"{adata.n_obs} cells x {adata.n_vars} genes after filtering")

# Select highly variable genes for downstream analysis
# These genes have high biological signal and drive the structure in the data
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata  # Save unfiltered data for later marker gene calculations
adata = adata[:, adata.var.highly_variable]  # Keep only HVGs for analysis

# Normalize and scale expression data
# Scale: Normalize to unit variance and zero mean
sc.pp.scale(adata, max_value=10)  # Clip values to ±10 to reduce impact of outliers

# Perform PCA (Principal Component Analysis) for initial dimensionality reduction
# Reduces 2000 genes to 30 principal components capturing major variance
sc.tl.pca(adata, svd_solver="arpack")

# Compute k-nearest neighbors graph in PCA space
# Uses 30 PCs and finds 10 nearest neighbors for each cell
# This defines the local structure needed for UMAP and clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

# Compute UMAP (Uniform Manifold Approximation and Projection)
# Creates 2D visualization preserving local and global cell relationships
sc.tl.umap(adata)

# Perform Leiden graph-based clustering algorithm
# Clusters cells based on the kNN graph; resolution=0.6 controls cluster granularity
sc.tl.leiden(adata, resolution=0.6)

# Visualize clusters in UMAP space and save plot
sc.pl.umap(adata, color="leiden", save="_clusters.png")

# Visualize malignancy status in UMAP space
# Shows separation between malignant cancer cells and non-malignant cells

sc.pl.umap(adata, color="malignant", save="_malignant_vs_normal.png")

# Identify marker genes for each cluster using Wilcoxon rank-sum test
# Finds genes significantly differentially expressed between each cluster and the rest
# Statistical test: non-parametric, robust to outliers common in single-cell data
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

# Visualize top marker genes for each cluster
sc.pl.rank_genes_groups(adata, n_genes=10, save="_markers.png")

# Define canonical marker genes for known cell types from literature
# These are well-established markers to validate our unsupervised clustering
canonical_markers = {
    "T cell": ["CD3D", "CD3E", "CD2"],
    "B cell": ["CD79A", "MS4A1"],
    "Macrophage": ["CD68", "CD163"],
    "Endothelial": ["PECAM1", "VWF"],
    "CAF": ["COL1A1", "DCN"],
    "NK": ["NKG7", "GNLY"],
    "Malignant": ["MLANA", "PMEL", "MITF"],
}

# Create dot plot showing expression of canonical markers across clusters
# x-axis: Leiden clusters | y-axis: marker genes | dot size/color: expression level
sc.pl.dotplot(adata, canonical_markers, groupby="leiden", use_raw=True, save="_dotplot.png")

# Cross-tabulation: Compare Leiden clusters with published cell type annotations
# Helps validate that our unsupervised clusters align with known cell types
comparison = pd.crosstab(adata.obs["leiden"], adata.obs["published_cell_type"])
print(comparison)
comparison.to_csv("results/cluster_vs_published_annotation.csv")


# Manual cell type annotation of clusters
# Based on:
# 1. Marker gene expression patterns from rank_genes_groups analysis
# 2. Canonical marker gene validation from dot plot
# 3. Comparison with published cell type annotations
# Note: Mapping may vary with different clustering parameters (resolution, etc.)
cluster_to_label = {
    # Fill in based on YOUR marker gene + crosstab results
    # Example mappings below — these will differ based on your specific analysis
    "0": "Malignant cells",
    "1": "T cells",
    "2": "Macrophages",
    # ... Add additional clusters based on what you observe
}

# Assign human-readable cell type labels to each cell based on cluster assignment
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_label)

# Create final UMAP visualization with annotated cell types
# This is the publication-ready visualization showing biological cell types
sc.pl.umap(adata, color="cell_type", save="_final_annotated.png",
           title="Melanoma TME — Annotated Cell Types")
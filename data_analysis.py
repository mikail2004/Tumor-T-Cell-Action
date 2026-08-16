# TTCA 
# Tumor T-Cell Action
# Mikail U.

# ----------- INTRODUCTION ------------ #
# Scanpy – Single-Cell Analysis in Python (GPU version available)
# Task: Identify non-exhausted T-Cells that can be used against tumors.
# Data: Series GSE108989 (Colorectal Cancer) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989)
# Tutorial: https://www.youtube.com/watch?v=5HuOGZEu2HY

"""
Notes:
1. Z-score measures how far a specific gene's expression level in an individual 
    cell is from the average expression of that same gene across all analyzed cells.
2. Cell cycle checkpoints: stages where a cell checks its health, size, 
    and DNA before moving forward to divide (mitosis).
3. T-cells get exhausted fighting cancer: the continuous, non-stop presence of tumors (unlike other diseases)
    overstimulates them. This chronic activation effects their genes and metabolism, 
    causing them to lose their ability to multiply (checkpoints) and kill cancer cells. 
4. Tumors steal vital nutrients from the locality, starving T-cells. 
"""

# ----------- REQUIREMENTS ------------ #
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

# Loading data (No un-compression needed)
raw_data_path = "GSE108989_CRC.TCell.S11138.count.txt.gz"
raw_data = pd.read_csv(raw_data_path, sep='\t')

# Read only the first 10 rows and inspect the first 10 columns
preview = pd.read_csv(raw_data_path, sep='\t', nrows=10)

print(f"Preview DataFrame shape: {preview.shape}")

print("\n ~First 10 rows and first 10 columns~ ")
print(preview.iloc[:10, :10])

print("\n ~Column names (first 10)~ ")
print(preview.columns[:10].tolist())

print("\n ~Data types of first 5 columns~ ")
print(preview.dtypes[:5])

# ----------- PRE-PROCESSING ------------ #
gene_symbols = raw_data['symbol'].astype('str').values
drops = raw_data.drop(columns=['geneID', 'symbol']) # Remove 'geneID' and 'symbol' from dataset

print(f"Matrix shape (Gene x Cells) : {drops.shape}")

# SC Expectation: Transpose so that rows = cells (barcodes), cols = genes 
# Originally it was rows = genes, cols = cells (barcodes)
adata = sc.AnnData(X=drops.T.values, 
                   obs=pd.DataFrame(index=drops.columns), # cells
                   var=pd.DataFrame(index=gene_symbols)) # genes

print(adata)

# Annotate mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate % of Mitochondrial (mt) Genes present
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Removing low quality cells 
sc.pp.filter_cells(adata, min_genes=200) # Cells with fewer than 200 genes
sc.pp.filter_genes(adata, min_cells=3) # Genes found in fewer than 3 cells

# Storing data so far
adata.raw = adata

# Normalize and Transform Data (???)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ----------- IDENTIFY BYSTANDER T-CELLS ------------ #
# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Run PCA (Principal Component Analysis) - Reduce size, keep important info, find hidden patterns
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Plotting UMAP 
# UMAP representing each cell (dot/marker) and grouping based on transcriptome (cell types) (colors).
# Leiden = algorithm for clustering. 
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])

# Find differential expression for genes per grouping
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Finally identify top 5 genes for each grouping using <gene_symbols> in dataset
# Use this to find genes that correspond to bystander T-cells 
print(pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5))
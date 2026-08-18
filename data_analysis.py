# TTCA 
# Tumor T-Cell Action
# Mikail U. Bukhari

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
5. Matrix Dimensions: In the raw dataset, each column is an individual cell and each row is a specific gene. 
    The number at the intersection is the count of mRNA transcripts detected 
    for that specific gene in that specific cell. 
6. Normalization: It corrects for differences in sequencing depth by scaling every cell 
    to have the same total sum of counts (I.E 10,000 total counts per cell). 
    This ensures fair comparisons (I.E if one cell had 10k counts and the other 20k).
7. Log-Transformation (log1p): It applies ln(1 + x) to compress extreme differences in scale (variance stabilization) 
    and handle zeros mathematically (ln(1+0) = 0). This prevents genes with high counts from dominating PCA 
    and clustering over low-count regulatory genes.
8. Adjusted P-value: Calculates probabilistic score of expected error across multiple runs.
9. (-log(10)): Makes extreme values discernable, especially for visualizing through graphs.
"""

# ----------- REQUIREMENTS ------------ #
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import gseapy as gp
import os

output_dir = "./enrichr_results"
os.makedirs(output_dir, exist_ok=True)

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

# Remove duplicate entries
adata.var_names_make_unique()
adata.obs_names_make_unique()

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

# ----------- ASSIGN LABELS TO GROUPINGS ------------ #
# 1. Creating marker dictionary for major T-cell phenotypes
marker_database = {
    'Tregs': {'FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'BATF', 'CCR8'},
    'Terminally Exhausted CD8+': {'HAVCR2', 'ENTPD1', 'CXCL13', 'LAYN', 'ITGAE', 'TOX', 'PDCD1'},
    'Cytotoxic / Bystander': {'NKG7', 'KLRG1', 'GZMH', 'CX3CR1', 'FGFBP2', 'GNLY', 'PRF1'},
    'Naive / Stem-like Memory': {'SELL', 'TCF7', 'LEF1', 'CCR7', 'IL7R'},
    'Effector / Activated': {'GZMK', 'CCL4', 'CCL5', 'JUNB', 'FOS', 'CD69'}
}

# Repeating identification of top genes but for actual analysis now
top_n = 15  # Using top 15 genes gives better matching precision than just 5
rank_names = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(top_n)

# Assign the best-matching label based on highest gene overlap count
cluster_labels = {}
for cluster in rank_names.columns:
    cluster_top_genes = set(rank_names[cluster])

    # Calculate overlap count with each biological state
    scores = {
        cell_state: len(cluster_top_genes.intersection(genes)) 
        for cell_state, genes in marker_database.items()
    }
    
    # Pick the state with the highest overlap
    best_match = max(scores, key=scores.get)
    
    # Fallback to Unclassified if no markers match
    cluster_labels[cluster] = best_match if scores[best_match] > 0 else f'Cluster {cluster} (Unassigned)'

# Map directly back into AnnData
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_labels)

# Plot and check assignments
print("Automated Cluster Assignments:")
for cl, label in cluster_labels.items():
    print(f"Cluster {cl} -> {label}")

sc.pl.umap(adata, color=['cell_type'])

# ----------- DIFFERENTIAL EXPRESSION ------------ #
# Compare bystander T-cells to exhausted T-cells to isolate exact genes that maintain bystander behavior
sc.tl.rank_genes_groups(
    adata,
    groupby='cell_type',
    groups=['Cytotoxic / Bystander'],
    reference='Terminally Exhausted CD8+',
    method='wilcoxon')

# Extract differential expression to dataframe
diff_exp = sc.get.rank_genes_groups_df(adata, group='Cytotoxic / Bystander')

# Find significant genes from resulting differential expression
sig_genes = diff_exp[(diff_exp['pvals_adj'] < 0.05) & (diff_exp['logfoldchanges'] > 1.0)]

print('Top upregulated genes between Bystander and Exhausted T-Cells')
print(sig_genes[['names', 'logfoldchanges', 'pvals_adj']].head(10))

# ----------- PATHWAY ENRICHMENT ------------ #
# Use differential expression results to find pathways for found significant genes
gene_list = sig_genes['names'].head(150).tolist()

pathways = gp.enrichr(
    gene_list=gene_list,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'GO_Biological_Process_2023'],
    organism='human',
    outdir=output_dir,
    no_plot=True,
    verbose=True)

# Print & Plot top enriched pathways found:
print(' ~ Pathways in Bystander T-Cells ~ ')
print(pathways.results[['Gene_set', 'Term', 'Adjusted P-value']].head(10))
gp.barplot(pathways.res2d, 
           title='Pathways in Bystander T-Cells', 
           column='Adjusted P-value', 
           group='Gene_set', 
           size=10, 
           top_term=5,
           ofname="./enrichr_results/PATHWAYS_PLOT.png")

# ----------- ENGAGER TARGET WORTHINESS ------------ #
# Find which bystander T-Cells have:
# 1. high cytotoxic machinery (GZMB, NKG7) 
# 2. high CD3 receptor availability (CD3E, CD3D)
# 3. low exhaustion checkpoints (PDCD1, HAVCR2, LAG3)

# Target vs. exhaustion gene sets
target_genes = [g for g in ['CD3E', 'CD3D', 'GZMB', 'PRF1', 'NKG7'] if g in adata.var_names]
exhaustion_genes = [g for g in ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'ENTPD1'] if g in adata.var_names]

# Score cell states
sc.tl.score_genes(adata, gene_list=target_genes, score_name='target_potential', use_raw=False)
sc.tl.score_genes(adata, gene_list=exhaustion_genes, score_name='exhaustion_score', use_raw=False)

# Calculate target worthiness
adata.obs['worthiness'] = adata.obs['target_potential'] - adata.obs['exhaustion_score']

# Visualize on UMAP
#sc.pl.umap(adata, color=['cell_type', 'target_potential', 'exhaustion_score', 'worthiness'])
sc.pl.umap(adata, 
           color=['cell_type', 'worthiness'],
           cmap='coolwarm', # Color palette
           vmin='p1', # Set minimum threshold to filter extreme values
           vmax='p99', # Set max threshold to filter extreme values
           ) 
# ----------- END OF CODE ------------ #
# TTCA 
# Tumor T-Cell Action
# Mikail U.

# ----------------------- #
# Scanpy – Single-Cell Analysis in Python (GPU version available)
# Task: Identify non-exhausted T-Cells that can be used against tumors.
# Data: Series GSE108989 (Colorectal Cancer)

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
# ----------------------- #


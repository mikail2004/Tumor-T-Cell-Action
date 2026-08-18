# Tumor-T-Cell-Action
The growth of tumors attracts a response by the immune system in the form of T-Cells. T-Cells struggle to deal with tumors due to a variety of complications, resulting in exhaustion; the destruction of their offensive abilities. Yet, there is a massive presence of bystander T-Cells that are yet to engage tumors. 

This project aims to use public data regarding the composition of cells at a tumor site (Colorectal Cancer) to identify not only bystander T-Cells but also compare them to their accompanying exhausted T-Cells to mark worthy cells capable of fighting tumors. Worthiness is based on possession of requisite genes to fight tumors whilst not having genes associated with exhaustion. All code is in `data_analysis.py`.

## Gene & Cell Analysis
Grouping Cells from the Original Dataset

![alt text](Plots/Plot_1.png)

Categorizing the Groupings based on Genes (Bystander, Exhausted found)

![alt text](Plots/Plot_2.png)

Most Worthy Cells 

![alt text](Plots/Plot_3.png)

Gene Pathway Information

![alt text](Plots/Plot_4.png)

## Dataset & Dataset Structure
[Series-GSE108989-Colorectal-Cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989)

| Column Name | Type | Description |
| :--- | :--- | :--- |
| **`geneID`** | `int64` | Entrez / NCBI Gene ID |
| **`symbol`** | `object` / `str` | HGNC Gene Symbol |
| **`<Cell_Barcode_1>`** | `int64` | Raw transcript (UMI) count in Cell 1 |
| **`<Cell_Barcode_2>`** | `int64` | Raw transcript (UMI) count in Cell 2 |
| **`...`** | `int64` | Raw transcript (UMI) count in Cell *N* |

* **Rows:** Individual genes indexed by identifier and gene symbol.
* **Columns (from column 3 onward):** Individual single-cell barcodes (`NP###-YYYYMMDD`).
* **Matrix Values:** Raw non-negative integer transcript counts (UMI).

## Requirements
```
scanpy==1.9.8
pandas==2.0.3
numpy==1.24.4
scipy==1.10.1
gseapy==0.14.0
matplotlib==3.7.5
leidenalg==0.10.2
umap-learn==0.5.7
h5py==3.11.0
```

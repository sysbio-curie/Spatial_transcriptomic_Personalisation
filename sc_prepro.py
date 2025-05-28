## Single cell RNAseq preprocessing to run before Cell2location and Spatialscope
# Author Agathe Sobkowicz started : 28/05


## Imports
import numpy as np
import pandas as pd
import scanpy as sc
import mygene

import matplotlib.pyplot as plt

run = "local"

## Read data
if run == "local":
    ad_sc = sc.read("/home/agathes/work/Sc_data/Sc_cellxgene.h5ad")
    ligand_recept = pd.read_csv(
        "/home/agathes/work/SpatialScope/extdata/ligand_receptors.txt", sep="\t"
    )
    sizek_df = pd.read_csv(
        "/home/agathes/work/PhysiBoSS_Personalisation/models/sizek_dictionary.csv"
    )
elif run == "cluster":
    ad_sc = sc.read(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene.h5ad"
    )
    ligand_recept = pd.read_csv(
        "/mnt/beegfs/home/asobkow1/persistent/spatialscope/SpatialScope/extdata/ligand_receptors.txt",
        sep="\t",
    )
    sizek_df = pd.read_csv(
        "/mnt/beegfs/home/asobkow1/persistent/data/models/sizek_dictionary.csv"
    )

## Restore raw data keeping normalised data in a "normal" layer
ad_sc.layers["normalized"] = ad_sc.X.copy()
ad_sc.X = ad_sc.raw.X

# print("Max count", ad_sc.X.max())
# print("Min count", ad_sc.X.min())
# print("Mean count", ad_sc.X.mean())

## Extract data
ad_sc = ad_sc[ad_sc.obs["tissue"] == "lung"].copy()  # keep only lung tissue
ad_sc = ad_sc[
    ad_sc.obs["harm_sample.type"].isin(["normal", "tumor"])
].copy()  # keep normal and tumor
ad_sc = ad_sc[ad_sc.obs["harm_study"] == "Qian et al"].copy()  # keep one study
ad_sc = ad_sc[ad_sc.obs["assay"] == "10x 3' v2"].copy()  # keep one technique

## Filter
sc.pp.filter_cells(ad_sc, min_counts=500)
sc.pp.filter_genes(ad_sc, min_cells=3)
ad_sc = ad_sc[:, ~np.array([_.startswith("MT-") for _ in ad_sc.var.index])]
ad_sc = ad_sc[:, ~np.array([_.startswith("mt-") for _ in ad_sc.var.index])]

## Create new sample_type + cell_type and Filter according to cell type
ad_sc.obs["sample_cell_type"] = (
    ad_sc.obs["harm_sample.type"].astype(str) + "-" + ad_sc.obs["cell_type"].astype(str)
)
# Keep only sample-cell categories with more than 500 cells
selected_sample_cell_type = (
    ad_sc.obs["sample_cell_type"]
    .value_counts()[ad_sc.obs["sample_cell_type"].value_counts() > 500]
    .index.tolist()
)
ad_sc = ad_sc[ad_sc.obs["sample_cell_type"].isin(selected_sample_cell_type), :].copy()
# Avoid over representation of one sample_cell_type with no more 3000 cells per sample_cell_type
ad_sc = ad_sc[
    ad_sc.obs.index.isin(
        ad_sc.obs.groupby("sample_cell_type")
        .apply(
            lambda x: (
                x.sample(frac=3000 / x.shape[0], replace=False)
                if x.shape[0] > 3000
                else x
            )
        )
        .reset_index(level=0, drop=True)
        .index
    )
].copy()
cell_type_column = "sample_cell_type"

print("Number of cells", ad_sc.shape[0])
print("Number of genes", ad_sc.shape[1])
print("Max count", ad_sc.X.max())
print("Min count", ad_sc.X.min())
print("Mean count", ad_sc.X.mean())

## Plot preprocessing figures
if ad_sc.X.max() < 20:
    ad_sc.X = np.exp(ad_sc.X) - 1
plt.hist(ad_sc.X.sum(1), bins=100)  # Plot function changes scale
plt.title("Total counts")
plt.xlabel("Cells")
plt.ylabel("Count")
plt.tight_layout()
if run == "local":
    plt.show()
elif run == "cluster":
    plt.gcf().savefig(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/sc_distribution.png",
        bbox_inches="tight",
    )
    plt.close()

ad_sc.var.index.name = "gene_ids_2"

# Save filtered not normalised values
if run == "local":
    ad_sc.write("/home/agathes/work/Sc_data/Sc_cellxgene_filtered.h5ad")
elif run == "cluster":
    ad_sc.write(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene_filtered.h5ad"
    )


## Renormalise and log1p after filtering
sc.pp.normalize_total(ad_sc, target_sum=30000)  # According to raw counts
sc.pp.log1p(ad_sc)


## Define a class of marker genes with higly variable, cell type markers,
## ligand recepter genes, model nodes and CAF marker genes
def convert_genes_ensembles(gene_list):
    mg = mygene.MyGeneInfo()
    ensembl_dict = mg.querymany(
        gene_list, scopes="symbol", fields="ensembl.gene", species="human"
    )
    ensembl_list = []
    for gene in ensembl_dict:
        ensembl = gene.get("ensembl")
        if isinstance(ensembl, dict):
            ensembl_list.append([gene["query"], ensembl["gene"]])
        elif isinstance(ensembl, list) and len(ensembl) > 0:
            ensembl_list.append([gene["query"], ensembl[0]["gene"]])
    return np.array(ensembl_list)


# Identify highly variable and marker genes
ad_sc.raw = ad_sc.copy()
sc.pp.highly_variable_genes(ad_sc, flavor="seurat_v3", n_top_genes=1000)
sc.tl.rank_genes_groups(ad_sc, groupby=cell_type_column, method="wilcoxon")
markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
markers = list(
    set(ad_sc.var.loc[ad_sc.var["highly_variable"] == 1].index) | set(markers)
)  # highly variable genes + cell type marker genes
# Ligand recepteur genes
ligand_recept = list(set(ligand_recept.melt()["value"].values))  # keep unique values
ligand_recept = [
    _.upper() for _ in ligand_recept
]  # upper case to correspond with humans
ligand_recept = list(convert_genes_ensembles(ligand_recept)[:, 1])
# Add gene from nodes in sizek model
sizek_genes = sorted(sizek_df["Genes"][sizek_df["Genes"] != "Not_assigned"].to_list())
sizek_ensembl = list(convert_genes_ensembles(sizek_genes)[:, 1])
# Add genes from Grout et al. 2022 papper
caf_genes = [
    "NOX4",
    "PDCD1",
    "CD274",
    "MMP2",
    "CLDN5",
    "PECAM1",
    "TFF3",
    "PROX1",
    "MCAM",
    "COX4I2",
    "DES",
    "CD10",
    "VEGFD",
    "IL6",
    "PI16",
    "CLU",
    "FAP",
    "ADH1B",
    "POSTN",
    "LRRC15",
    "GREM1",
    "ACTA2",
    "MMP2",
    "MYH11",
    "ACTA2",
    "COL3A1",
    "BGN",
    "TCF21",
    "COL9A1",
    "COL27A1",
    "COL4A2",
    "COL11A1",
    "COL12A1",
]
caf_ensembl = list(convert_genes_ensembles(caf_genes)[:, 1])
other_genes = [
    "PDCD1",
    "CD274",
    "CLDN5",
    "PECAM1",
    "TFF3",
    "PROX1",
    "MCAM",
    "COX4I2",
    "DES",
]
other_ensembl = list(convert_genes_ensembles(other_genes)[:, 1])
markers = list(set(markers + ligand_recept + sizek_ensembl + caf_genes + other_genes))
print("Number of selected markers :", len(markers))

# Add marker column
ad_sc.var.loc[ad_sc.var.index.isin(markers), "Marker"] = True
ad_sc.var["Marker"] = ad_sc.var["Marker"].fillna(False)
ad_sc.var["highly_variable"] = ad_sc.var["Marker"]

# Log data
sc.pp.pca(ad_sc)
sc.pp.neighbors(ad_sc)
sc.tl.umap(ad_sc)

# Umpa plot
fig, axs = plt.subplots(1, 1, figsize=(12, 10))
sc.pl.umap(
    ad_sc,
    color="cell_type",
    size=15,
    frameon=False,
    show=False,
    ax=axs,
)

plt.tight_layout()
if run == "local":
    plt.show()
elif run == "cluster":
    plt.gcf().savefig(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/sc_umap.png",
        bbox_inches="tight",
    )
    plt.close()

## Write processed data
if run == "local":
    ad_sc.write("/home/agathes/work/Sc_data/Sc_cellxgene_normalised.h5ad")
elif run == "cluster":
    ad_sc.write(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene_normalised.h5ad"
    )

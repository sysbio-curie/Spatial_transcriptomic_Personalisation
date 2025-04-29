## Single cell RNAseq preprocessing to run Spatialscope
# Author Agathe Sobkowicz started : 29/04
# Based on the demo of Spatial scope on Human heart Visium data


## Imports
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scanpy as sc

import sys

sys.path.append("../src")
# from utils import *

import warnings

warnings.filterwarnings("ignore")

run = "local"

## Read data
if run == "local":
    ad_sc = sc.read("/home/agathes/work/Sc_data/Sc_cellxgene.h5ad")
elif run == "cluster":
    ad_sc = sc.read(
        "/mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene.h5ad"
    )

ad_sc.layers["normalized"] = ad_sc.X.copy()
ad_sc.X = ad_sc.raw.X
# Keep only lung tissue without blood samples
ad_sc = ad_sc[(ad_sc.obs["tissue"] == "lung")].copy()  # only lung
ad_sc = ad_sc[
    (ad_sc.obs["harm_sample.type"] != "blood")
].copy()  # keep cancer and normal

ad_sc.obs["sample_cell_type"] = (
    ad_sc.obs["harm_sample.type"].astype(str) + "-" + ad_sc.obs["cell_type"].astype(str)
)
# Keep only sample-cell categories with more than 500 cells
selected_sample_cell_type = (
    ad_sc.obs["sample_cell_type"]
    .value_counts()[ad_sc.obs["sample_cell_type"].value_counts() > 500]
    .index.tolist()
)
adata_ref = ad_sc[
    ad_sc.obs["sample_cell_type"].isin(selected_sample_cell_type), :
].copy()

cell_type_column = "sample_cell_type"

sc.pp.filter_cells(ad_sc, min_counts=500)
sc.pp.filter_cells(ad_sc, max_counts=20000)
sc.pp.filter_genes(ad_sc, min_cells=5)
ad_sc = ad_sc[:, ~np.array([_.startswith("MT-") for _ in ad_sc.var.index])]
ad_sc = ad_sc[:, ~np.array([_.startswith("mt-") for _ in ad_sc.var.index])]

ad_sc = ad_sc[
    ad_sc.obs["harm_study"] == "Qian et al"
].copy()  # reduce batch effect due to different studies
ad_sc = ad_sc[
    ad_sc.obs["assay"] == "10x 3' v2"
].copy()  # reduce batche effect between experiments
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
].copy()  # avoid over representation of one sample_cell_type

print(ad_sc.X.max(), ad_sc.shape)

# Plot preprocessing figures
sns.set_context("paper", font_scale=1.6)
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
sc.pl.umap(
    ad_sc,
    color=cell_type_column,
    size=15,
    frameon=False,
    show=False,
    ax=axs,
    legend_loc="on data",
)
plt.tight_layout()

if ad_sc.X.max() < 20:
    ad_sc.X = np.exp(ad_sc.X) - 1
plt.hist(ad_sc.X.sum(1), bins=100)
plt.show()

# Normalise
ad_sc.raw = ad_sc.copy()
sc.pp.normalize_total(ad_sc, target_sum=2000)

# Identify import genes
sc.pp.highly_variable_genes(ad_sc, flavor="seurat_v3", n_top_genes=1000)
sc.tl.rank_genes_groups(ad_sc, groupby=cell_type_column, method="wilcoxon")
markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
markers = list(
    set(ad_sc.var.loc[ad_sc.var["highly_variable"] == 1].index) | set(markers)
)  # highly variable genes + cell type marker genes

ligand_recept = list(
    set(
        pd.read_csv(
            "/home/agathes/work/SpatialScope/extdata/ligand_receptors.txt", sep="\t"
        )
        .melt()["value"]
        .values
    )
)
# if scRNA-seq reference is from human tissue, run following code to make gene name consistent
ligand_recept = [_.upper() for _ in ligand_recept]
add_genes = (
    "DCN GSN PDGFRA\
RGS5 ABCC9 KCNJ8\
MYH11 TAGLN ACTA2\
GPAM FASN LEP\
MSLN WT1 BNC1\
VWF PECAM1 CDH5\
CD14 C1QA CD68\
CD8A IL7R CD40LG\
NPPA MYL7 MYL4\
MYH7 MYL2 FHL2\
DLC1 EBF1 SOX5\
FHL1 CNN1 MYH9\
CRYAB NDUFA4 COX7C\
PCDH7 FHL2 MYH7\
PRELID2 GRXCR2 AC107068.2\
MYH6 NPPA MYL4\
CNN1 MYH9 DUSP27\
CKM COX41L NDUFA4\
DLC1 PLA2GS MAML2\
HAMP SLIT3 ALDH1A2\
POSTN TNC FAP\
SCN7A BMPER ACSM1\
FBLN2 PCOLCE2 LINC01133\
CD36 EGFLAM FTL1\
CFH ID4 KCNT2\
PTX3 OSMR IL6ST\
DCN PTX3 C1QA\
NOTCH1 NOTCH2 NOTCH3 NOTCH4 DLL1 DLL4 JAG1 JAG2\
CDH5 SEMA3G ACKR1 MYH11".split()
    + ["TNNT2", "PLN", "MYH7", "MYL2", "IRX3", "IRX5", "MASP1"]
)  # some important genes that we interested

markers = markers + add_genes + ligand_recept
len(markers)

# Add marker column
ad_sc.var.loc[ad_sc.var.index.isin(markers), "Marker"] = True
ad_sc.var["Marker"] = ad_sc.var["Marker"].fillna(False)
ad_sc.var["highly_variable"] = ad_sc.var["Marker"]

# Log data
sc.pp.log1p(ad_sc)
sc.pp.pca(ad_sc, svd_solver="arpack", n_comps=30, use_highly_variable=True)

ad_sc.X.max()

sc.pp.neighbors(ad_sc, metric="cosine", n_neighbors=30, n_pcs=30)
sc.tl.umap(ad_sc, min_dist=0.5, spread=1, maxiter=60)

# preprocessing plots
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
sc.pl.umap(
    ad_sc,
    color=cell_type_column,
    size=15,
    frameon=False,
    show=False,
    ax=axs,
    legend_loc="on data",
)
plt.tight_layout()

## Write processed data
if run == "local":
    ad_sc.write("/home/agathes/work/Sc_data//Sc_cellxgene_processed_lung_batch.h5ad")

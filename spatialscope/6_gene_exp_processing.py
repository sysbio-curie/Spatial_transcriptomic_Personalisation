# Transforme Spatialscope cell type identification
# and gene expression decomposition results
# started 3/06/2025 author: Agathe Sobkowicz

# Import
import os
import scanpy as sc
from src_custom_image.utils import *


## Define paths
results_dir = "/home/agathes/work/results"
deconv_results_dir = os.path.join(
    results_dir,
    "spatialscope_deconv/third_run/deconv_normalised",
)
sc_data_dir = os.path.join(
    results_dir,
    "spatialscope_deconv/third_run/Sc_cellxgene_normalised.h5ad",
)
slides = [
    "17P04394",
    "18P06762",
    # "18P08140",
    # "18P06593",
    # "18P03122",
    # "18P02831",
    # "17P02529",
]

for slide in slides:
    slide_deconv_dir = os.path.join(deconv_results_dir, slide)
    sp_adata = sc.read(os.path.join(slide_deconv_dir, "sp_adata.h5ad"))

    # Merge gene expression decomposition results
    END_INDEX = sp_adata.obs.shape[0] - (sp_adata.obs["cell_count"] == 0).sum()
    CHUNK_SIZE = 200
    spot_ranges = []
    start = 0
    while start < END_INDEX:
        end = start + CHUNK_SIZE
        if end > END_INDEX:
            end = END_INDEX
        spot_ranges.append((start, end))
        start = end

    print("Found spot ranges:", spot_ranges)

    gene_adata = [
        sc.read(f"{slide_deconv_dir}/generated_cells_spot{start}_{end}.h5ad")
        for start, end in spot_ranges
    ]
    gene_adata = gene_adata[0].concatenate(
        gene_adata[1:], batch_key="_", uns_merge="unique", index_unique=None
    )
    gene_adata.uns = sp_adata.uns
    sc_adata = sc.read_h5ad(sc_data_dir)
    shared_vars = gene_adata.var_names.intersection(sc_adata.var_names)
    gene_adata.var.loc[shared_vars, "feature_name"] = sc_adata.var.loc[
        shared_vars, "feature_name"
    ]
    gene_adata.write(os.path.join(slide_deconv_dir, f"decomposed_gene_exp.h5ad"))

    # Plot identified cell types
    # fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=100)
    # fig.patch.set_alpha(0.1)decomposed_
    # PlotVisiumCells(
    #     sp_adata, "discrete_label_ct", size=1.5, alpha_img=0.4, lw=0.3, ax=ax
    # )

    # fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=100)
    # fig.patch.set_alpha(0.1)
    # PlotVisiumCells(
    #     sp_adata, "discrete_label_ct", size=0.4, alpha_img=0.3, lw=0.3, ax=ax
    # )
    # ax.set_xlim(9500, 11500)
    # ax.set_ylim(10000, 12000)
    # ax.axis("off")
    # plt.tight_layout()
    # plt.show()

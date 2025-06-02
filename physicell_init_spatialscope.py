# Transform Celloc deconvolution results to PhysiCell initilisation
# started 09/04/2025 author: Agathe Sobkowicz

# Import
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from spatialscope.src_custom_image.utils import *


## Define paths to deconvolution results
results_dir = "/home/agathes/work/results"
deconv_results_dir = os.path.join(
    results_dir,
    "spatialscope_deconv/third_run/deconv_normalised",
)
slides = [
    "17P04394",
    "18P06762",
    "18P08140",
    "18P06593",
    "18P03122",
    "18P02831",
    "17P02529",
]
model_init_results_dir = os.path.join(
    results_dir, "model_initialisation/second_run/spatialscope_deconv"
)
os.makedirs(model_init_results_dir, exist_ok=True)

## Import data
for slide in slides:
    slide_deconv_dir = os.path.join(deconv_results_dir, slide)
    sp_adata = sc.read(os.path.join(slide_deconv_dir, "sp_adata.h5ad"))
    cell_sp_data = sp_adata.uns["cell_locations"]
    slide_model_init = os.path.join(model_init_results_dir, slide)
    os.makedirs(slide_model_init, exist_ok=True)
    cell_sp_data.to_csv(
        os.path.join(slide_model_init, f"cell_init_spatialscope_{slide}.csv")
    )

    fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=100)
    fig.patch.set_alpha(0.1)
    PlotVisiumCells(
        sp_adata, "discrete_label_ct", size=1.5, alpha_img=0.4, lw=0.3, ax=ax
    )

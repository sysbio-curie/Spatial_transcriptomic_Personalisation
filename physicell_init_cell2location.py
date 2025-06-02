# Transform Celloc deconvolution results to PhysiCell initilisation
# started 09/04/2025 author: Agathe Sobkowicz

# Import
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad


## Define paths to deconvolution results
results_dir = "/home/agathes/work/results"
deconv_results_dir = os.path.join(
    results_dir,
    "Cell2loc_deconv/Sencond_run_changed_filter_sc_sp",
)
slides = [
    "17P02529",
    "18P06762",
    "18P08140",
    "17P04394",
    "18P06593",
    "18P03122",
    "18P02831",
]
model_init_results_dir = os.path.join(
    results_dir, "model_initialisation/second_run/cell2loc_deconv"
)
os.makedirs(model_init_results_dir, exist_ok=True)

## Import data
for slide in slides:
    slide_deconv_dir = os.path.join(deconv_results_dir, slide)
    adata_map = sc.read_h5ad(
        os.path.join(slide_deconv_dir, "cell2location_map/st_map.h5ad")
    )

    slide_model_init_results = os.path.join(model_init_results_dir, slide)
    os.makedirs(slide_model_init_results, exist_ok=True)

    ## Extract coordinates and dominant cell type of each spot
    coords_df = pd.DataFrame(
        adata_map.obsm["spatial"], columns=["x", "y"], index=adata_map.obs_names
    )
    cell_types = [
        name.split("sf_")[-1]
        for name in adata_map.obsm["means_cell_abundance_w_sf"].columns
    ]
    cell_abundance_df = pd.DataFrame(adata_map.obsm["means_cell_abundance_w_sf"])
    cell_abundance_df.columns = cell_types
    dominant_cell_df = cell_abundance_df.idxmax(axis=1).astype(str)

    model_init_df = pd.concat([coords_df, dominant_cell_df], axis=1)
    model_init_df.columns = ["x", "y", "dominant_cell_type"]
    model_init_df.to_csv(
        os.path.join(slide_model_init_results, f"cell_init_{slide}.csv")
    )
    ## Visualise dominant cell types in all spot of the slide
    sns.countplot(data=model_init_df, x="dominant_cell_type")
    plt.title("Count of Different Cell Types")
    plt.xlabel("Cell Type")
    plt.ylabel("Count")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            slide_model_init_results, f"Plot_dominant_cells_slide_{slide}.png"
        ),
    )
    plt.show()

    adata_init = ad.AnnData(X=model_init_df[["x", "y"]].values)
    adata_init.obs["dominant_cell_type"] = dominant_cell_df.to_string()
    adata_init.write(os.path.join(slide_model_init_results, f"cell_init_{slide}.h5ad"))

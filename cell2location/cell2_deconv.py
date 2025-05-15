# Trainning Cell2location model for lusc v2 deconvolution with cell2ocation - cluster implementation
# started 02/04/2025 author: Agathe Sobkowicz
# st : spatial transcriptomics, sc: single-cell

# Import
import os
import numpy as np
import scanpy as sc
import cell2location as cell2loc
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
import celltypist as ct
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description="Process a parameter from SLURM.")
parser.add_argument(
    "--param", type=str, required=True, help="Parameter from SLURM job array"
)
args = parser.parse_args()
slide = args.param
print(f"Running script on slide: {slide}")

# Define paths to data, models an results
current_dir = os.getcwd()
vis_data_dir = os.path.join(current_dir, "data/LUSC_v2")
vis_slide_dir = os.path.join(vis_data_dir, slide)

results_dir = os.path.join(current_dir, "results")
os.makedirs(results_dir, exist_ok=True)

ref_results_dir = os.path.join(results_dir, "reference_signatures")
os.makedirs(ref_results_dir, exist_ok=True)

slide_results_dir = os.path.join(results_dir, slide)
os.makedirs(slide_results_dir, exist_ok=True)

st_prepro_results_dir = os.path.join(slide_results_dir, "spatial_preprocessing")
os.makedirs(st_prepro_results_dir, exist_ok=True)

map_results_dir = os.path.join(slide_results_dir, "cell2location_map")
os.makedirs(map_results_dir, exist_ok=True)


## Preprocessing functions
def tot_filter_genes(adata, min_counts=200, pct_min_cells=0.05):
    """All filters applied on genes of st data"""
    adata.var["Gene outlier"] = (
        ~sc.pp.filter_genes(adata, min_counts=min_counts, inplace=False)[0]
        | ~sc.pp.filter_genes(
            adata, min_cells=int(pct_min_cells * adata.n_obs), inplace=False
        )[0]
    )  # remove genes expressed weakly or in less than 5% of cells
    # returns tuple, first element is a True False array where
    # True are spots above min count to keep,
    # INVERSED array so that True is an outlier

    print(adata.var["Gene outlier"].value_counts())

    return adata


def tot_filter_cells(adata, min_counts=200):
    """All filters applied on cells of st data"""
    adata.obs["Cell outlier"] = ~sc.pp.filter_cells(
        adata, min_counts=min_counts, inplace=False
    )[0]

    print(adata.obs["Cell outlier"].value_counts())

    return adata


def st_processing(adata, slide):
    """Preprocessing steps: QC netrics, filter, normalise, dimension reduction and cluster"""

    adata.var["mt"] = [gene.startswith("MT-") for gene in adata.var_names]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=[20])
    adata.obs["mt_frac"] = (
        adata[:, adata.var["mt"]].X.sum(1).A.squeeze() / adata.obs["total_counts"]
    )

    # Filter
    adata = tot_filter_genes(adata)
    adata = adata[:, (~adata.var["Gene outlier"])].copy()

    adata = tot_filter_cells(adata)
    sc.pl.spatial(
        adata,
        img_key="hires",
        color="Cell outlier",
        title="Outlier",
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(st_prepro_results_dir, f"outlier_spatial_plot_{slide}.png")
    )
    plt.close()
    adata = adata[(~adata.obs["Cell outlier"])].copy()

    # Normalising on a copy of data object as value can not be normalised for deconvolution
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, inplace=True)  # Shifted algorithm
    sc.pp.log1p(adata_norm)  # log1p transform

    # Dimesion reduction
    sc.pp.highly_variable_genes(
        adata_norm, flavor="seurat", n_top_genes=2000, inplace=True
    )
    sc.pp.pca(adata_norm, n_comps=50, mask_var="highly_variable", svd_solver="arpack")
    sc.pp.neighbors(adata_norm)
    sc.tl.umap(adata_norm)
    sc.tl.leiden(adata_norm, key_added="leiden_res0_5", resolution=0.5)

    # Visualise
    sc.pl.spatial(
        adata_norm,
        img_key="hires",
        show=False,
    )
    plt.gcf().savefig(os.path.join(st_prepro_results_dir, f"histo_{slide}.png"))
    plt.close()
    sc.pl.spatial(
        adata_norm,
        img_key="hires",
        color=["total_counts", "n_genes_by_counts"],
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(st_prepro_results_dir, f"count_spatial_plot_{slide}.png")
    )
    plt.close()
    sc.pl.umap(
        adata_norm,
        color=["total_counts", "n_genes_by_counts"],
        legend_loc="on data",
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(st_prepro_results_dir, f"count_umap_plot_{slide}.png")
    )
    plt.close()
    sc.pl.spatial(
        adata_norm,
        img_key="hires",
        color="leiden_res0_5",
        alpha=0.7,
        size=1.3,
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(st_prepro_results_dir, f"cluster_umap_plot_{slide}.png")
    )
    plt.close()

    ct.models.download_models(model=["Immune_All_High.pkl"], force_update=True)
    model = ct.models.Model.load(model="Immune_All_High.pkl")
    predictions = ct.annotate(
        adata_norm,
        model="Immune_All_High.pkl",
        majority_voting=True,
        over_clustering="leiden_res0_5",
    )
    adata_norm = predictions.to_adata()

    sc.pl.umap(adata_norm, color=["predicted_labels"], show=False)
    plt.gcf().savefig(
        os.path.join(
            st_prepro_results_dir,
            f"predicted_annotation_cluster_umap_{slide}.png",  # predicted_
        ),
        bbox_inches="tight",
    )
    plt.close()
    sc.pl.umap(adata_norm, color=["majority_voting"], show=False)
    plt.gcf().savefig(
        os.path.join(
            st_prepro_results_dir, f"voted_annotation_cluster_umap_{slide}.png"
        ),
        bbox_inches="tight",
    )
    plt.close()
    sc.pl.spatial(
        adata_norm,
        img_key="hires",
        color=["predicted_labels"],
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(
            st_prepro_results_dir,
            f"predicted_annotation_clusters_spatial_plot_{slide}.png",
        ),
        bbox_inches="tight",
    )
    plt.close()
    sc.pl.spatial(
        adata_norm,
        img_key="hires",
        color=["majority_voting"],
        show=False,
    )
    plt.gcf().savefig(
        os.path.join(
            st_prepro_results_dir, f"voted_annotation_clusters_spatial_plot_{slide}.png"
        ),
        bbox_inches="tight",
    )
    plt.close()
    return adata


## Regression model, load saved model and output h5ad
print("Importing Regression model")

adata_file = os.path.join(ref_results_dir, "sc_ref.h5ad")
adata_ref = sc.read_h5ad(adata_file)
mod = cell2loc.models.RegressionModel.load(
    os.path.join(ref_results_dir, "reg_model"), adata_ref
)

# Export estimated expression in each cluster
if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
inf_aver.iloc[0:5, 0:5]


## Import and Preprocess St data
print(f"Started analysing slide: {slide}")
adata_vis = sc.read_visium(vis_slide_dir)

# Add column to obs with sample name
adata_vis.obs["sample"] = list(adata_vis.uns["spatial"].keys())[0]
# For deconvolution gene/var names are replaced by ENSEMBL ID to match sc data
adata_vis.var["SYMBOL"] = adata_vis.var_names
adata_vis.var_names_make_unique()
if adata_vis.obsm["spatial"].dtype == "object":
    adata_vis.obsm["spatial"] = np.array(adata_vis.obsm["spatial"], dtype=float)

# Pre processing
adata_vis = st_processing(adata_vis, slide)

# Rename genes to ENSEMBL ID for matching between single cell and spatial data
adata_vis.var.set_index("gene_ids", drop=True, inplace=True)

# Remove MT genes in ST data because it represents aretefacts in sc data
adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
adata_vis.obsm["MT"] = adata_vis[
    :, adata_vis.var["MT_gene"].values
].X.toarray()  # Keep MT counts in obsm
adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]

# Find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver_slide = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2loc.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
Cell2locMod = cell2loc.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver_slide,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20,
)
Cell2locMod.view_anndata_setup()

Cell2locMod.train(
    max_epochs=30000,
    # train using full data (batch_size=None)
    batch_size=None,
    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1,
)

# plot ELBO loss history during training, # removing first 100 epochs from the plot
Cell2locMod.plot_history(1)
plt.legend(labels=["full data training"])
plt.gcf().savefig(os.path.join(map_results_dir, f"cell2loc_train_loss_{slide}.png"))
plt.close()

# Export the estimated cell abundance (summary of the posterior distribution).
adata_vis = Cell2locMod.export_posterior(
    adata_vis,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": Cell2locMod.adata.n_obs,
    },
)

# Save model
Cell2locMod.save(os.path.join(map_results_dir, "cell2loc_mod"), overwrite=True)

# Save anndata object with results
adata_file = os.path.join(map_results_dir, "st_map.h5ad")
adata_vis.write(adata_file)
adata_file

# # Load saved model and output h5ad
# adata_file = os.path.join(map_results_dir,"st_map.h5ad")
# adata_vis = sc.read_h5ad(adata_file)
# Cell2locMod = cell2loc.models.Cell2location.load(os.path.join(map_results_dir, "cell2loc_mod"), adata_vis)

Cell2locMod.plot_QC()
plt.gcf().savefig(os.path.join(map_results_dir, f"cell2loc_mod_qc_{slide}.png"))
plt.close()

fig = Cell2locMod.plot_spatial_QC_across_batches()
plt.gcf().savefig(os.path.join(map_results_dir, f"cell2loc_mod_qc_mapped_{slide}.png"))
plt.close()

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm[
    "q05_cell_abundance_w_sf"
]

# plot in spatial coordinates
with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [4.5, 5]}):
    sc.pl.spatial(
        adata_vis,
        cmap="magma",
        # show first 8 cell types
        color=[
            "normal-T cell",
            "normal-epithelial cell",
            "normal-fibroblast",
            "tumor-T cell",
            "tumor-epithelial cell",
            "tumor-fibroblast",
            "tumor-malignant cell",
        ],
        ncols=4,
        size=1.3,
        img_key="hires",
        # limit color scale at 99.2% quantile of cell abundance
        vmin=0,
        vmax="p99.2",
        show=False,
    )
    plt.gcf().savefig(os.path.join(map_results_dir, f"cell_abundance_{slide}.png"))
    plt.close()

# up to 6 clusters
clust_labels = [
    "normal-T cell",
    "normal-fibroblast",
    "tumor-T cell",
    "tumor-fibroblast",
    "tumor-malignant cell",
]
clust_col = [
    "" + str(i) for i in clust_labels
]  # in case column names differ from labels

with mpl.rc_context({"figure.figsize": (15, 15)}):
    fig = cell2loc.plt.plot_spatial(
        adata=adata_vis,
        # labels to show on a plot
        color=clust_col,
        labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style="fast",
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position="right",
    )
    plt.gcf().savefig(
        os.path.join(map_results_dir, f"cell_abundance_multiple_{slide}.png")
    )
    plt.close()

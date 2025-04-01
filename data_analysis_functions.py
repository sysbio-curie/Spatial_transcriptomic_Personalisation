# Data analysis functions
# started: 21/02/2025 author: Agathe Sobkowicz


import os
import seaborn as sb
import matplotlib.pyplot as plt
import scanpy as sc
import celltypist as ct
import cell2location as cell2loc


def print_tree(startpath, indent=""):
    """Prints a tree-like folder structure, ignoring hidden files and folders."""
    for root, _, files in os.walk(startpath):
        level = root.replace(startpath, "").count(os.sep)
        indent = "│   " * level + "├── "
        print(f"{indent}{os.path.basename(root)}/")
        subindent = "│   " * (level + 1) + "├── "
        for f in files:
            if not f.startswith("."):
                print(f"{subindent}{f}")


def plot_count_gene_dist(adata):
    """PLot de the total count and gene per count distribution for given data"""
    fig, axs = plt.subplots(1, 5, figsize=(25, 4))
    fig.suptitle("QC")
    sb.histplot(
        adata.obs["total_counts"],
        kde=False,
        bins=160,
        ax=axs[0],
    )  # Number of obs/spot with given amount of counts
    axs[0].set_ylabel("Number of Observations")  # Change y-axis label
    sb.histplot(
        adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
        kde=False,
        bins=20,
        ax=axs[1],
    )
    axs[1].set_ylabel("Number of Observations")  # Change y-axis label
    sb.histplot(
        adata.obs["n_genes_by_counts"], kde=False, bins=55, ax=axs[2]
    )  # Number of obs/spot with given amount of unique genes per spot/obs
    axs[2].set_ylabel("Number of genes per observation")  # Change y-axis label
    sb.histplot(
        adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
        kde=False,
        bins=20,
        ax=axs[3],
    )
    axs[3].set_ylabel("Number of genes per observation")  # Change y-axis label

    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        ax=axs[4],
        show=False,
    )
    # Remove the extra color bar (if present)
    for cbar in axs[4].figure.axes:
        if hasattr(cbar, "collections") and len(cbar.collections) > 0:
            for collection in cbar.collections:
                if collection.colorbar is not None:
                    collection.colorbar.remove()  # Remove extra color bar
    # Add color bar separately and reposition it
    cbar = fig.colorbar(
        axs[4].collections[0],
        ax=axs[4],
        orientation="vertical",
        fraction=0.05,
        pad=0.05,
    )
    cbar.set_label("pct_counts_mt")

    # Adjust layout to prevent overlap
    plt.tight_layout()


def plot_count_gene_dist_normal(adata):
    _, axes = plt.subplots(1, 4, figsize=(18, 4))
    sb.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    sb.histplot(adata.X.sum(1), bins=25, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm - Total counts")
    sb.histplot(adata.obs["n_genes_by_counts"], bins=55, kde=False, ax=axes[2])
    axes[2].set_title("Genes per count")
    sb.histplot((adata.X > 0).sum(1), bins=20, kde=False, ax=axes[3])
    axes[3].set_title("Shifted logarithm - Genes per count")


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


def st_processing(adata, slide, st_prepro_results_dir):
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
    sc.pp.normalize_total(adata_norm, inplace=True, target_sum=1e4)  # Shifted algorithm
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


def sc_preprocessing(adata_ref):
    """Preprocessing of single cell data used to define cell type references"""
    # Convert downloaded preprocessed data back to raw data saving normalised values
    adata_ref.layers["normalized"] = adata_ref.X.copy()
    adata_ref.X = adata_ref.raw.X
    # Keep only lung tissue without blood samples
    adata_ref = adata_ref[(adata_ref.obs["tissue"] == "lung")].copy()  # only lung
    adata_ref = adata_ref[
        (adata_ref.obs["harm_sample.type"] != "blood")
    ].copy()  # keep cancer and normal

    # Create category storing both sample and cell types
    adata_ref.obs["sample_cell_type"] = (
        adata_ref.obs["harm_sample.type"].astype(str)
        + "-"
        + adata_ref.obs["cell_type"].astype(str)
    )
    # Keep only sample-cell categories with more than 500 cells
    selected_sample_cell_type = (
        adata_ref.obs["sample_cell_type"]
        .value_counts()[adata_ref.obs["sample_cell_type"].value_counts() > 500]
        .index.tolist()
    )
    adata_ref = adata_ref[
        adata_ref.obs["sample_cell_type"].isin(selected_sample_cell_type), :
    ].copy()

    # Genes/var are already names by ENSEBL ID, which is necessary for matching between single cell and spatial data
    # Change column name to match with the spatial dataset
    adata_ref.var.index.name = "gene_ids_2"

    # Mitochondrial genes
    adata_ref.var["mt"] = adata_ref.var["feature_name"].str.startswith("MT-")
    # Ribosomal genes
    adata_ref.var["ribo"] = adata_ref.var["feature_name"].str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata_ref, qc_vars=["mt", "ribo"], inplace=True, percent_top=[20], log1p=True
    )

    # Filter out none relevant genes
    selected_genes = cell2loc.utils.filtering.filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.2,
    )
    adata_ref = adata_ref[:, selected_genes].copy()

    return adata_ref

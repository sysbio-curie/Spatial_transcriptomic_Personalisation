# Trainning Regression model for lusc v2 deconvolution with cell2ocation - cluster implementation
# started 02/04/2025 author: Agathe Sobkowicz
# st : spatial transcriptomics, sc: single-cell
# depredicated on 28/05/2025 to match sc preprocessing between spatialscope and cell2location

# Import
import os
import scanpy as sc
import cell2location as cell2loc
import matplotlib.pyplot as plt


## Define paths to data and results
current_dir = os.getcwd()

sc_data_dir = os.path.join(current_dir, "data/Sc_data/Sc_cellxgene.h5ad")
results_dir = os.path.join(current_dir, "results")
os.makedirs(results_dir, exist_ok=True)

ref_results_dir = os.path.join(results_dir, "reference_signatures")
os.makedirs(ref_results_dir, exist_ok=True)


## Preprocessing functions
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


## Import and Preprocess single-cell data
print("Importing and preprocessing single cell reference dataset :")
adata_ref = sc.read_h5ad(sc_data_dir)
adata_ref = sc_preprocessing(adata_ref)


## Regression model
print("Setting up data for Regression model")

# Prepare anndata for the regression model
cell2loc.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key="donor_id",  # 10X reaction / sample / batch
    labels_key="sample_cell_type",  # sample-cell type, covariate used for constructing signatures
    categorical_covariate_keys=[
        "assay",
        "harm_study",
    ],  # multiplicative technical effects (platform, 3' vs 5')
)

# Create the regression model
RegModel = cell2loc.models.RegressionModel(adata_ref)

# View anndata_setup as a sanity check
RegModel.view_anndata_setup()

print("Started training Regression model")
# Train model
RegModel.train(max_epochs=250)

# Check training
RegModel.plot_history()
plt.gcf().savefig(os.path.join(ref_results_dir, "reg_model_loss_plot.png"))
plt.close()

# Export estimated cell abundance - summary of posterior distribution
adata_ref = RegModel.export_posterior(
    adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
)

# Save model
RegModel.save(os.path.join(ref_results_dir, "reg_model"), overwrite=True)
# Save anndata object with results
adata_file = os.path.join(ref_results_dir, "sc_ref.h5ad")
adata_ref.write(adata_file)
adata_file

adata_ref = RegModel.export_posterior(
    adata_ref,
    use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={"batch_size": 2500},
)

# QC of prior distribution
RegModel.plot_QC(summary_name="q50")

# # Load saved model and output h5ad
# adata_file = os.path.join(ref_results_dir,"sc_ref.h5ad")
# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2loc.models.RegressionModel.load(os.path.join(ref_results_dir,"reg_model"), adata_ref)

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

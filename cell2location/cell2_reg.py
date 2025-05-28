# Trainning Regression model for lusc v2 deconvolution with cell2ocation - cluster implementation
# started 02/04/2025 author: Agathe Sobkowicz
# st : spatial transcriptomics, sc: single-cell

# Import
import os
import scanpy as sc
import cell2location as cell2loc
import matplotlib.pyplot as plt

current_dir = os.getcwd()

## Import data
sc_data_dir = os.path.join(current_dir, "data/Sc_data/Sc_cellxgene_filtered.h5ad")
adata_ref = sc.read_h5ad(sc_data_dir)

## Paths to results
results_dir = os.path.join(current_dir, "results")
os.makedirs(results_dir, exist_ok=True)

ref_results_dir = os.path.join(results_dir, "reference_signatures")
os.makedirs(ref_results_dir, exist_ok=True)

## Regression model
print("Setting up data for Regression model")

# Prepare anndata for the regression model
cell2loc.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key="donor_id",  # 10X reaction / sample / batch
    labels_key="sample_cell_type",  # sample-cell type, covariate used for constructing signatures
    categorical_covariate_keys=[],  # ["assay","harm_study"]
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
# mod = cell2loc.models.RegressionModel.load(os.path.join(ref_results_dir,"reg_model"), adata_ref

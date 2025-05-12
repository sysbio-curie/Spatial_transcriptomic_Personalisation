## Imports
import argparse
import scanpy as sc
import os

parser = argparse.ArgumentParser(description='simulation sour_sep')
parser.add_argument('--ST_Data', type=str, help='ST data path', default=None)
args = parser.parse_args()

## Read Data
par_dir = os.path.dirname(args.ST_Data)
adata_sp = sc.read_visium(par_dir)

## Change names of vars and column name while saving the original var column to have enselbe_ids as var names
adata_sp.var["feature_name"] = adata_sp.var_names
adata_sp.var.rename(columns={"gene_ids": "ensembl_id"}, inplace=True)
adata_sp.var.set_index("ensembl_id", drop=True, inplace=True)

## Save file
adata_sp.write(os.path.join(par_dir,"filtered_feature_bc_matrix_prepro.h5ad"))

## Spatial transcriptomic data preprocessing to run Spatialscope
# Author Agathe Sobkowicz started : ?/05/2025

## Imports
import argparse
import scanpy as sc
import os
import glob
import openslide

from image_alignment import *

use = "local"

## Read spatial transcripotmic data and the two image data files
if use == "cluster":
    parser = argparse.ArgumentParser(description="simulation sour_sep")
    parser.add_argument(
        "--Slide_Path",
        type=str,
        help="Path to spatial and image data of one slide",
        default=None,
    )
    args = parser.parse_args()

    adata_sp = sc.read_visium(args.Slide_Path)
    ndpi_file = glob.glob(os.path.join(args.Slide_Path, "*.ndpi"))[0]  # wsi
    tiff_file = glob.glob(os.path.join(args.Slide_Path, "*.tiff"))[0]  # template
elif use == "local":
    path_to_slide = "/home/agathes/work/LUSC_v2/18P06762_test"
    adata_sp = sc.read_visium(path_to_slide)
    ndpi_file = glob.glob(os.path.join(path_to_slide, "*.ndpi"))[0]  # wsi
    tiff_file = glob.glob(os.path.join(path_to_slide, "*.tiff"))[0]  # template


## Change names of vars and column name while saving the original var column to have enselbe_ids as var names
adata_sp.var["feature_name"] = adata_sp.var_names
adata_sp.var.rename(columns={"gene_ids": "ensembl_id"}, inplace=True)
adata_sp.var.set_index("ensembl_id", drop=True, inplace=True)

## Filter cells and genes
sc.pp.filter_cells(adata_sp, min_counts=100)
sc.pp.filter_genes(adata_sp, min_cells=3)
#  Remove MT genes in ST data because it represents aretefacts in sc data
adata_sp.var["MT_gene"] = adata_sp.var["feature_name"].str.upper().str.startswith("MT-")
adata_sp.obsm["MT"] = adata_sp[
    :, adata_sp.var["MT_gene"].values
].X.toarray()  # Keep MT counts in obsm
adata_sp = adata_sp[:, ~adata_sp.var["MT_gene"].values]


def correspondance_sp_image(adata_sp, ndpi_file, tiff_file, analysis_level: int = 1):
    ## Align ndpi and tiff images
    # Get coresponding transformation matrix and
    # coordinate to crop ndpi image at different levels to get tiff image
    transformation_matrix, cropped_coordinates_matrix = template_matching_workflow(
        ndpi_file, tiff_file
    )

    ## Change cooridnates of spatial data to match cropped ndpi image
    spatial_coordinates = np.array(adata_sp.obsm["spatial"], dtype=int)
    # Saved transformation matrix is defined to get coordinates in level 0 of ndpi image
    new_spatial_coordinates = transform_spatial_coordinates(
        spatial_coordinates, (transformation_matrix / (2**analysis_level))
    )

    ## Crop ndpi image, keep only level of analysis and save as tiff image for futur use
    xmin_level_0, ymin_level_0, _, _ = (
        cropped_coordinates_matrix[cropped_coordinates_matrix["level"] == 0]
        .iloc[:, 1:]
        .astype(int)
        .values[0]
    )
    xmin, ymin, w, h = (
        cropped_coordinates_matrix[
            cropped_coordinates_matrix["level"] == analysis_level
        ]
        .iloc[:, 1:]
        .astype(int)
        .values[0]
    )

    # Adapt spatial coordinates
    # new_spatial_coordinates[:, 0] = new_spatial_coordinates[:, 0] - xmin
    # new_spatial_coordinates[:, 1] = new_spatial_coordinates[:, 1] - ymin
    adata_sp.obsm["spatial"] = new_spatial_coordinates

    slide = openslide.OpenSlide(ndpi_file)
    region = slide.read_region(
        (xmin_level_0, ymin_level_0), level=analysis_level, size=(w, h)
    )
    #    region = slide.read_region((0, 0), level=analysis_level, size=slide.level_dimensions[analysis_level])
    region_rgb = region.convert("RGB")

    # Add custom high resolution image corresponding to new spatial coordinates
    slide_id = list(adata_sp.uns["spatial"].keys())[0]
    adata_sp.uns["spatial"][slide_id]["images"]["custom"] = np.array(region_rgb)
    adata_sp.uns["spatial"][slide_id]["scalefactors"]["tissue_custom_scalef"] = 1

    return adata_sp, region_rgb


adata_sp, region_rgb = correspondance_sp_image(adata_sp, ndpi_file, tiff_file)
# adata_sp = correspondance_sp_image(adata_sp, ndpi_file, tiff_file)

## Save preprocessed image and spatial transcriptomic data
if use == "cluster":
    region_rgb.save(os.path.join(args.Slide_Path, "prepro_image.tiff"), format="TIFF")
    adata_sp.write(
        os.path.join(args.Slide_Path, "filtered_feature_bc_matrix_prepro.h5ad")
    )
elif use == "local":
    region_rgb.save(os.path.join(path_to_slide, "prepro_image.tiff"), format="TIFF")
    adata_sp.write(
        os.path.join(path_to_slide, "filtered_feature_bc_matrix_prepro.h5ad")
    )

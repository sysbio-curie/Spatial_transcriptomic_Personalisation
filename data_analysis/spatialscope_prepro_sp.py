## Imports
import argparse
import scanpy as sc
import os
import glob
import openslide

from image_alignement import *


parser = argparse.ArgumentParser(description="simulation sour_sep")
parser.add_argument(
    "--Slide_Path",
    type=str,
    help="Path to spatial and image data of one slide",
    default=None,
)
args = parser.parse_args()

## Read spatial transcripotmic data and the two image data files
adata_sp = sc.read_visium(args.Slide_Path)
ndpi_file = glob.glob(os.path.join(args.Slide_Path, "*.ndpi"))  # wsi
tiff_file = glob.glob(
    os.path.join(args.Slide_Path, "*.tiff")
)  # cropped template for visium

## Change names of vars and column name while saving the original var column to have enselbe_ids as var names
adata_sp.var["feature_name"] = adata_sp.var_names
adata_sp.var.rename(columns={"gene_ids": "ensembl_id"}, inplace=True)
adata_sp.var.set_index("ensembl_id", drop=True, inplace=True)


def correspondance_sp_image(adata_sp, ndpi_file, tiff_file, analysis_level: int = 1):
    ## Align ndpi and tiff images
    # Get coresponding transformation matrix and
    # coordinate to crop ndpi image at different levels to get tiff image
    transformation_matrix, cropped_coordinates_matrix = template_matching_workflow(
        ndpi_file, tiff_file
    )

    ## Change cooridnates of spatial data to match cropped ndpi image
    spatial_coordinates = adata_sp.obsm["spatial"].values
    # Saved transformation matrix is defined to get coordinates in level 0 of ndpi image
    new_spatial_coordinates = transform_spatial_coordinates(
        spatial_coordinates, (transformation_matrix / (2**analysis_level))
    )
    adata_sp.obsm["spatial"] = new_spatial_coordinates

    ## Crop ndpi image, keep only lelvel of analysis and save as tiff image for futur use
    xmin, ymin, w, h = cropped_coordinates_matrix[
        cropped_coordinates_matrix["level"] == analysis_level
    ].iloc[:, 1:]
    slide = openslide.OpenSlide(ndpi_file)
    region = slide.read_region((xmin, ymin), level=analysis_level, size=(w, h))
    region_rgb = region.convert("RGB")
    return adata_sp, region_rgb


adata_sp, region_rgb = correspondance_sp_image(adata_sp, ndpi_file, tiff_file)

## Save image
region_rgb.save(os.path.join(args.Slide_Path, "prepro_image.tiff"), format="TIFF")
## Save file
adata_sp.write(os.path.join(args.Slide_Path, "filtered_feature_bc_matrix_prepro.h5ad"))

# df = pd.read_csv(tissue_position_path)
# spatial_coordinates = df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
# new_spatial_coordinates = transform_spatial_coordinates(
#     spatial_coordinates, transformation_matrix
# )
# df[["pxl_col_in_fullres", "pxl_row_in_fullres"]] = new_spatial_coordinates
# df.to_csv(os.path.join(sample_folder, "tissue_positions.csv"), index=False)

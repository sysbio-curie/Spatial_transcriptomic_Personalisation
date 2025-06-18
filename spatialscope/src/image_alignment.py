# Functions to find the correspondance between the tiff and ndpi image in the LUSC v2 data
# Extracted from image_alignement.py written by Loic
# Agathe started : 09/05/2025


from matplotlib import pyplot as plt
from openslide import OpenSlide
from tifffile import imread
from PIL import ImageOps, Image
import cv2
import numpy as np
import pandas as pd

from numpy import ndarray


def multiscale_template_matching(
    reference_image: ndarray,
    template_image: ndarray,
    min_scale: float = 0.1,
    max_scale: float = 0.3,
    scale_step: float = 0.01,
) -> dict:
    best_result = None
    best_score = float("-inf")

    for flip_x, flip_y in [(False, False), (True, False), (False, True), (True, True)]:
        flipped_template = template_image.copy()

        if flip_x:
            flipped_template = cv2.flip(flipped_template, 1)
        if flip_y:
            flipped_template = cv2.flip(flipped_template, 0)

        for scale_factor in np.arange(min_scale, max_scale + scale_step, scale_step):
            scaled_template = cv2.resize(
                flipped_template, (0, 0), fx=scale_factor, fy=scale_factor
            )

            result = cv2.matchTemplate(
                reference_image, scaled_template, cv2.TM_CCOEFF_NORMED
            )
            _, max_val, _, max_loc = cv2.minMaxLoc(result)

            if max_val > best_score:
                best_score = max_val
                best_result = {
                    "location": max_loc,
                    "scale_factor": scale_factor,
                    "confidence_score": max_val,
                    "flip_x": flip_x,
                    "flip_y": flip_y,
                    "scaled_template_size": scaled_template.shape,
                    "template_size": template_image.shape,
                }

    return best_result


# Assuming best_result is the output of multiscale_template_matching
def get_transformation_matrix_from_template_matching(best_result: dict):
    # Extract best match properties
    scale_factor = best_result["scale_factor"]
    flip_x = best_result["flip_x"]
    flip_y = best_result["flip_y"]
    match_location = best_result["location"]

    # Compute scaling factors, including mirroring if applicable
    scale_x = scale_factor * (-1 if flip_x else 1)
    scale_y = scale_factor * (-1 if flip_y else 1)
    scaling_mirror_matrix = np.array([[scale_x, 0, 0], [0, scale_y, 0], [0, 0, 1]])

    # Adjust translation to account for flipping
    adjusted_location_x = match_location[0]
    adjusted_location_y = match_location[1]

    if flip_x:
        adjusted_location_x += int(best_result["template_size"][1] * scale_factor)
    if flip_y:
        adjusted_location_y += int(best_result["template_size"][0] * scale_factor)

    # Translation matrix using the adjusted location
    translation_matrix = np.array(
        [[1, 0, adjusted_location_x], [0, 1, adjusted_location_y], [0, 0, 1]]
    )

    # Final transformation matrix (translation * scaling/mirroring)
    transformation_matrix = translation_matrix @ scaling_mirror_matrix
    return transformation_matrix


def transform_point(point: ndarray, matrix: ndarray) -> ndarray:
    point_homogeneous = np.append(point, 1)  # Convert to homogeneous coordinates
    transformed_point = matrix @ point_homogeneous
    return transformed_point[:2].astype(int)


def transform_spatial_coordinates(points: ndarray, matrix: ndarray) -> ndarray:
    homogeneous_coordinate = np.hstack([points, np.ones((points.shape[0], 1))])
    transformed_points = np.dot(homogeneous_coordinate, matrix.T)
    return transformed_points[:, :2].astype(int)


def coordinate_correspondence(match_result: ndarray, transformation_matrix: ndarray):
    # Get cropped coordinates on full resolution image
    top_left = transform_point(np.asarray([0, 0]), transformation_matrix)
    bottom_right = transform_point(
        np.asarray(match_result["template_size"][::-1]), transformation_matrix
    )
    x_min, y_min = np.min([bottom_right, top_left], axis=0)
    w, h = np.abs(bottom_right - top_left)
    return x_min, y_min, w, h


def plot_template_location(
    reference_image: ndarray, template_image: ndarray, match_result: dict
) -> None:
    transformation_matrix = get_transformation_matrix_from_template_matching(
        match_result
    )
    x_min, y_min, w, h = coordinate_correspondence(match_result, transformation_matrix)
    # print([x_min, y_min, w, h])

    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    axs[0].imshow(template_image, cmap="gray")
    axs[0].set_title("Template (non oriented)")
    axs[1].imshow(reference_image, cmap="gray")
    axs[1].add_patch(
        plt.Rectangle((x_min, y_min), w, h, edgecolor="r", facecolor="none")
    )
    axs[1].set_title("Matched location")
    plt.show()


def template_matching_workflow(wsi_path: str, template_path: str, level: int = 4):
    # Load the WSI image
    wsi = OpenSlide(wsi_path)
    raw_image = wsi.read_region((0, 0), level=level, size=wsi.level_dimensions[level])
    gray_image = np.asarray(ImageOps.grayscale(raw_image))

    # Load the template image
    template = imread(template_path)
    gray_template = np.asarray(ImageOps.grayscale(Image.fromarray(template)))

    # Perform multi-scale template matching
    match_result = multiscale_template_matching(gray_image, gray_template)

    # Plot template location
    plot_template_location(gray_image, gray_template, match_result)

    # Compute transformation matrix
    transformation_matrix = get_transformation_matrix_from_template_matching(
        match_result
    )

    cropped_coordinates_matrix = pd.DataFrame(
        columns=["level", "xmin", "ymin", "w", "h"]
    )

    x_min, y_min, w, h = coordinate_correspondence(match_result, transformation_matrix)
    cropped_coordinates_matrix.loc[0] = [level, int(x_min), int(y_min), int(w), int(h)]

    transformation_matrix = transformation_matrix * (2 ** (level - 2))

    x_min, y_min, w, h = coordinate_correspondence(match_result, transformation_matrix)
    cropped_coordinates_matrix.loc[1] = [2, int(x_min), int(y_min), int(w), int(h)]

    transformation_matrix = transformation_matrix * 2

    x_min, y_min, w, h = coordinate_correspondence(match_result, transformation_matrix)
    cropped_coordinates_matrix.loc[2] = [1, int(x_min), int(y_min), int(w), int(h)]

    transformation_matrix = transformation_matrix * 2

    x_min, y_min, w, h = coordinate_correspondence(match_result, transformation_matrix)
    cropped_coordinates_matrix.loc[3] = [0, int(x_min), int(y_min), int(w), int(h)]

    # Update transformation matrix to reflect the cropping
    new_transformation_matrix = transformation_matrix.copy()
    new_transformation_matrix[0, 2] -= x_min
    new_transformation_matrix[1, 2] -= y_min

    return new_transformation_matrix, cropped_coordinates_matrix

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import numpy as np
from skimage import draw
import math
from numpy import matlib
from loguru import logger

import constants


def compute_nucleus_and_cytoplasm_line_segments(nucleus_mask: np.ndarray, cytoplasm_mask: np.ndarray,
                                                nucleus_centroid: np.ndarray, x_slope: float, y_slope: float,
                                                max_cell_radius=None,
                                                image_width=None, image_height=None) -> (np.ndarray, np.ndarray):
    """
    Given the nucleus and cytoplasm mask and x, y slopes of a line, compute the segment
    of the line that falls withing the nucleus and the segment that falls within the cytoplasm
    """
    max_cell_radius = max_cell_radius or constants.analysis_config['MAX_CELL_RADIUS']
    image_width = image_width or constants.dataset_config['IMAGE_WIDTH']
    image_height = image_height or constants.dataset_config['IMAGE_HEIGHT']
    cytoplasm_segment = np.zeros(max_cell_radius)
    nucleus_segment = np.zeros(max_cell_radius)

    # compute for each point of the line (length MAX_CELL_RADIUS) if it falls within the nucleus or the cytoplasm
    for point in range(max_cell_radius):
        x = int(round(nucleus_centroid[0] + point * x_slope))
        y = int(round(nucleus_centroid[1] + point * y_slope))
        if not (x < 0 or x > (image_width - 1) or
                y < 0 or y > (image_height - 1)):
            cytoplasm_segment[point] = cytoplasm_mask[y, x]
            nucleus_segment[point] = nucleus_mask[y, x]
        else:
            # TODO use a global counter to report this kind of error
            pass
    return nucleus_segment, cytoplasm_segment


def compute_edge_points(nucleus_segment: np.ndarray, cytoplasm_segment: np.ndarray) -> (int, int):
    """
    Given line segments, return points that fall on the nucleus and on the cytoplasm edges
    :param nucleus_segment:
    :param cytoplasm_segment:
    :return:
    """
    nucleus_points = np.where(nucleus_segment == 1)
    if len(nucleus_points[0]) > 0:
        nucleus_edge_point = nucleus_points[0][-1]  # TODO : should be = nucleus_points[0][-1] + 1
    else:
        nucleus_edge_point = 1
    cytoplasm_points = np.where(cytoplasm_segment == 1)
    if len(cytoplasm_points[0]) > 0:
        cytoplasm_edge_point = cytoplasm_points[0][-1]  # TODO : should be = cytoplasm_points[0][-1] + 1
    else:
        cytoplasm_edge_point = 1
    return nucleus_edge_point, cytoplasm_edge_point


def compute_contour_points(nucleus_mask, nucleus_centroid, cytoplasm_mask, num_contours=None,
                           max_cell_radius=None, image_width=None, image_height=None) -> np.ndarray:
    """
    Computes contours within the cytoplasm that form concentric isolines between the nucleus
    and the cytoplasm periphery. Each contour is defined as (x,y) coordinates of 360 points
    Returns an array with coordinates of points for each of num_contours contours
    """
    num_contours = num_contours or constants.analysis_config['NUM_CONTOURS']
    contour_points = np.zeros((360, num_contours, 2))
    x_max = np.max(np.where(cytoplasm_mask == 1)[0])  # slow ?
    y_max = np.max(np.where(cytoplasm_mask == 1)[1])

    for degree in range(360):
        angle = degree * 2 * math.pi / 360
        x_slope, y_slope = math.sin(angle), math.cos(angle)
        segments = compute_nucleus_and_cytoplasm_line_segments(nucleus_mask, cytoplasm_mask, nucleus_centroid, x_slope,
                                                               y_slope, max_cell_radius, image_width, image_height)
        nucleus_edge_point, cytoplasm_edge_point = compute_edge_points(segments[0], segments[1])
        segment_length = cytoplasm_edge_point - nucleus_edge_point

        for index in range(num_contours):
            point = nucleus_edge_point + segment_length * index // num_contours  # TODO : (i) int division, (ii) should be (num_contours - 1)
            x = int(round(nucleus_centroid[0] + point * x_slope).astype(int))
            y = int(round(nucleus_centroid[1] + point * y_slope).astype(int))
            contour_points[degree, index, :] = [x, y]
            if x > x_max or y > y_max:
                #     TODO use a global counter to report this kind of error
                pass
    return contour_points


def create_mask(row_coordinates, column_coordinates, image_shape) -> np.ndarray:
    rr, cc = draw.polygon(row_coordinates, column_coordinates, shape=image_shape)
    mask = np.zeros(image_shape, dtype=np.bool)  # TODO should be int ?
    mask[rr, cc] = True  # TODO should be 1 ?
    return mask


def compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points, num_contours=None) -> np.ndarray:
    assert (cytoplasm_mask.shape == nucleus_mask.shape), 'cytoplasm mask and nucleus mask do not have the same size'
    num_contours = num_contours or constants.analysis_config['NUM_CONTOURS']
    cytoplasm_truth_mask = cytoplasm_mask.astype(bool)
    cell_mask_distance_map = np.zeros((cytoplasm_mask.shape[0], cytoplasm_mask.shape[1]), dtype=np.int)
    for index in range(num_contours):
        if index == 0:
            peripheral_mask = cytoplasm_mask
        else:
            contour_num = num_contours - index
            peripheral_mask = create_mask(contour_points[:, contour_num, 1], contour_points[:, contour_num, 0],
                                          (cytoplasm_mask.shape[0], cytoplasm_mask.shape[1]))
            peripheral_mask &= cytoplasm_truth_mask
        cell_mask_distance_map[(peripheral_mask == 1)] = index + 1

    return cell_mask_distance_map


def compute_all_distances_to_nucleus_centroid(nucleus_centroid: np.ndarray, image_width=None,
                                              image_height=None) -> np.ndarray:
    """
    Compute distances between all points and nucleus_centroid in a IMAGE_WIDTH x IMAGE_HEIGHT matrix
    """
    image_width = image_width or constants.dataset_config['IMAGE_WIDTH']
    image_height = image_height or constants.dataset_config['IMAGE_HEIGHT']
    if image_width != image_height:
        raise IndexError("Implemented only for images with IMAGE_WIDTH == IMAGE_HEIGHT, {} != {}", image_width,
                         image_height)

    i, j = np.meshgrid(np.arange(image_height), np.arange(image_width))
    dist = np.sqrt((i - nucleus_centroid[0]) ** 2 + (j - nucleus_centroid[1]) ** 2)

    return dist

def compute_all_distances_to_nucleus_centroid3d(heightmap: np.ndarray, nucleus_centroid: np.ndarray,
                                                image_width=None, image_height=None) -> np.ndarray:
    """
    Compute distances within the cytoplasm between all points and nucleus_centroid in a
    IMAGE_WIDTH x IMAGE_HEIGHT x cytoplasm_height matrix (max height of the cytoplasm)
    """
    image_width = image_width or constants.dataset_config['IMAGE_WIDTH']
    image_height = image_height or constants.dataset_config['IMAGE_HEIGHT']
    cytoplsam_height = np.max(heightmap)
    nucleus_centroid_z = heightmap[nucleus_centroid[0], nucleus_centroid[1]] // 2
    if image_width != image_height:
        raise IndexError("Implemented only for images with IMAGE_WIDTH == IMAGE_HEIGHT, {} != {}",
                         image_width, image_height)

    i, j, k = np.meshgrid(np.arange(image_height), np.arange(image_width), np.arange(cytoplsam_height))
    dist = np.sqrt((j - nucleus_centroid[0]) ** 2 +
                   (i - nucleus_centroid[1]) ** 2 +
                   (k - nucleus_centroid_z) ** 2)

    return dist


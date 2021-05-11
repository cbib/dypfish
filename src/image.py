#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math
from typing import List

import numpy as np
from loguru import logger
from scipy import spatial
from sklearn.metrics.pairwise import pairwise_distances

import constants
import helpers
import image_processing as ip
from constants import CELL_MASK_DISTANCE_PATH_SUFFIX
from constants import CELL_MASK_PATH_SUFFIX
from constants import CYTOPLASM_MASK_PATH_SUFFIX
from constants import MTOC_LEADING_EDGE_SUFFIX
from constants import MTOC_POSITION_PATH_SUFFIX
from constants import NUCLEUS_CENTROID_PATH_SUFFIX
from constants import NUCLEUS_MASK_PATH_SUFFIX
from constants import PERIPHERAL_QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX
from constants import PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX
from constants import QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX
from constants import QUADRANT_DENSITIES_PATH_SUFFIX
from repository import Repository


class Image(object):
    """ Represents a generic image. It has at least a cell_mask, a nucleus_mask and a nucleus_centroid """

    def __init__(self, repository: Repository, image_path: str):
        """
        Sets instance variable path to the image location (image_path) in an HDF5 file
        """
        self._repository = repository
        if not self._repository.is_present(image_path):
            raise LookupError("Cannot init, image %s is absent from repository" % image_path)
        if not self._repository.is_present(image_path + CELL_MASK_PATH_SUFFIX) or \
                not self._repository.is_present(image_path + NUCLEUS_MASK_PATH_SUFFIX) or \
                not self._repository.is_include(image_path, image_path + NUCLEUS_CENTROID_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)
        self._path = image_path

    def __eq__(self, img2) -> bool:
        return self._repository == img2._repository and self._path == img2._path

    def get_cell_mask(self) -> np.ndarray:
        return np.array(self._repository.get(self._path + CELL_MASK_PATH_SUFFIX)).astype(int)

    def get_nucleus_mask(self) -> np.ndarray:
        return np.array(self._repository.get(self._path + NUCLEUS_MASK_PATH_SUFFIX)).astype(int)

    def get_nucleus_centroid(self) -> np.ndarray:
        return np.around(
            self._repository.get_by_regex(self._path, self._path + NUCLEUS_CENTROID_PATH_SUFFIX)).flatten().astype(
            np.int)

    def compute_nucleus_area(self):
        """compute nucleus surface in pixel using nucleus mask"""
        nucleus_mask = self.get_nucleus_mask()
        area = nucleus_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        return area

    def compute_cell_area(self):
        """compute cell surface in real pixel size using cell mask"""
        cell_mask = self.get_cell_mask()
        area = cell_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        return area

    def compute_areas_from_periphery(self):
        # compute cell area per isoline in cytoplasm
        cytoplasm_mask = self.get_cytoplasm_mask()
        distance_map = self.get_cell_mask_distance_map()
        areas = np.zeros(constants.analysis_config['NUM_CONTOURS'])
        for i in range(constants.analysis_config['NUM_CONTOURS']):
            areas[i] = cytoplasm_mask[(distance_map <= i+1)].sum() * helpers.surface_coeff()
        return areas

    @helpers.checkpoint_decorator(CYTOPLASM_MASK_PATH_SUFFIX, float)
    def get_cytoplasm_mask(self) -> np.ndarray:
        return self.compute_cytoplasm_mask()

    def compute_cytoplasm_mask(self) -> np.ndarray:
        cell_mask = self.get_cell_mask()
        nucleus_mask = 1 - self.get_nucleus_mask()
        return np.multiply(cell_mask, nucleus_mask)

    def compute_peripheral_mask(self) -> np.ndarray:
        cytoplasm_mask = self.get_cytoplasm_mask()
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = ((cell_mask_dist_map > 0) &
                                  (cell_mask_dist_map <= peripheral_fraction_threshold)).astype(int)
        return np.multiply(peripheral_binary_mask, cytoplasm_mask)

    @helpers.checkpoint_decorator(CELL_MASK_DISTANCE_PATH_SUFFIX, dtype=np.int)
    def get_cell_mask_distance_map(self):
        return self.compute_cell_mask_distance_map()

    def compute_cell_mask_distance_map(self):
        cytoplasm_mask = self.get_cytoplasm_mask()
        nucleus_mask = self.get_nucleus_mask()
        nucleus_centroid = self.get_nucleus_centroid().transpose()  # TODO why the transpose ??

        # for each degree, analyse the line segment between the nucleus and the periphery
        contour_points = ip.compute_contour_points(nucleus_mask, nucleus_centroid, cytoplasm_mask,
                                                   num_contours=constants.analysis_config['NUM_CONTOURS'])
        cell_mask_distance_map = ip.compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points)
        return cell_mask_distance_map

    def is_in_cytoplasm(self, position: List) -> bool:
        cytoplasm_mask = self.get_cytoplasm_mask()
        return cytoplasm_mask[position[0], position[1]] == 1

    def compute_cell_diameter(self) -> float:
        cell_mask = self.get_cell_mask()
        cell_mask_points = np.argwhere(cell_mask == 1)
        convex_hull = spatial.ConvexHull(cell_mask_points)
        boundary_points = cell_mask_points[convex_hull.vertices]
        d = np.max(pairwise_distances(boundary_points))
        return d

    def compute_cytoplasmic_coordinates_peripheral_distance(self, coordinates) -> np.ndarray:
        """
        Return an array of integer distances to the periphery for coordinates
        Returns np.nan for those coordinates that ate outside of the cytoplasm
        """
        logger.info("Computing {} coordinates 2D peripheral distances for {} coordinates", len(coordinates), self._path)
        assert coordinates.shape[1] == 2, "2D coordinates needed for distance to the periphery"

        peripheral_distance_map = self.get_cell_mask_distance_map()
        distances = np.array([])
        for c in coordinates:
            if self.is_in_cytoplasm(c[::-1]):
                distances = np.append(distances, peripheral_distance_map[c[1], c[0]])
            else:
                distances = np.append(distances, np.nan)
        logger.info("  found {} coordinates outside of cytoplasm", len(distances[np.isnan(distances)]))
        return distances

    def compute_signal_from_periphery(self) -> np.ndarray:
        raise NotImplementedError


class ImageWithMTOC(Image):
    """ Represents an image with identified MTOC coordinates and the MTOC quadrant """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_present(
            path + MTOC_POSITION_PATH_SUFFIX)  # TODO : check if this is needed : and repo.is_present(path + MTOC_LEADING_EDGE_SUFFIX)

    def mtoc_is_in_leading_edge(self) -> bool:
        return np.array(self._repository.get(self._path + MTOC_LEADING_EDGE_SUFFIX)) == 1

    def get_mtoc_position(self) -> np.ndarray:
        return np.around(self._repository.get(self._path + MTOC_POSITION_PATH_SUFFIX)).flatten().astype(np.int)

    def rotate_quadrant_mask(self, degree, slices_num=4, nucleus_centroid=None, mtoc_position=None,
                             cell_mask=None, image_width=None, image_height=None):
        """
        Computes the quadrant mask (slices are quadrants if slices_num == 4) of a cell anchored in
        the MTOC position, rotated by a degree. The original quadrant of the MTOC (before rotation)
        is defined by two lines anchored at the nucleus_centroid and at 45 degrees relative to
        the nucleus_centroid-MTOC line
        Returns a mask where each pixel of the cell has it's quadrant number
        """
        # TODO problem with quadrant num when mtoc and nucleus centroid are in the top right quadrant of the whole image
        # TODO problem with quadrant num when mtoc and nucleus centroid are in the bottom left quadrant of the whole image
        nucleus_centroid = nucleus_centroid or self.get_nucleus_centroid()
        mtoc_position = mtoc_position or self.get_mtoc_position()
        if cell_mask is None: cell_mask = self.get_cell_mask()
        image_width = image_width or constants.dataset_config['IMAGE_WIDTH']
        image_height = image_height or constants.dataset_config['IMAGE_HEIGHT']
        right_point = helpers.rotate_point(self, nucleus_centroid, mtoc_position, degree)

        s = helpers.slope_from_points(self, nucleus_centroid, right_point)
        radians = np.arctan(s)  # angle wrt to x axis
        xx, yy = np.meshgrid(np.array(range(0, image_width)) - nucleus_centroid[0],
                             np.array(range(0, image_height)) - nucleus_centroid[1])
        rotated_xx, rotated_yy = helpers.rotate_meshgrid(xx, yy, -radians)

        # Arbitrarily assign number to each slice
        # This trunc of pi is used for compatibility between macOS and Linux vs how they deal
        # with infinite decimals (0.9999999999999...)
        pi = format(math.pi, '.10f')
        pi = float(pi)

        sliceno = ((pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / ((8 / slices_num) * pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == slices_num + 1] = slices_num  # int conversion sometimes rounds the value
        quadrant_mask[cell_mask == 0] = 0
        return quadrant_mask.astype(int)

    @helpers.checkpoint_decorator(PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_quadrants_densities(self, quadrants_num=4):
        return self.compute_quadrant_densities(quadrants_num=quadrants_num, peripheral_flag=True)

    @helpers.checkpoint_decorator(QUADRANT_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_quadrants_densities(self, quadrants_num=4):
        return self.compute_quadrant_densities(quadrants_num=quadrants_num)

    @helpers.checkpoint_decorator(PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_quadrants_densities(self, quadrants_num=4, peripheral_flag=True):
        return self.compute_quadrant_densities(quadrants_num=quadrants_num, peripheral_flag=peripheral_flag)

    @helpers.checkpoint_decorator(QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_quadrants_and_slices_densities(self, quadrants_num=4, stripes=3, stripes_flag=True):
        return self.compute_quadrant_densities(quadrants_num=quadrants_num, peripheral_flag=False,
                                               stripes=stripes, stripes_flag=stripes_flag)

    @helpers.checkpoint_decorator(PERIPHERAL_QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_quadrants_and_slices_densities(self, quadrants_num=4, peripheral_flag=True,
                                                      stripes=3, stripes_flag=True):
        return self.compute_quadrant_densities(quadrants_num=quadrants_num, peripheral_flag=peripheral_flag,
                                               stripes=stripes, stripes_flag=stripes_flag)

    def get_or_compute_quadrant_densities(self, quadrants_num=4, peripheral_flag=False, stripes=3, stripes_flag=False):
        if (not peripheral_flag) and (not stripes_flag):
            density_per_quadrant = self.get_quadrants_densities(quadrants_num)
        if peripheral_flag and (not stripes_flag):
            density_per_quadrant = self.get_peripheral_quadrants_densities(quadrants_num, peripheral_flag)
        if (not peripheral_flag) and stripes_flag:
            density_per_quadrant = self.get_quadrants_and_slices_densities(quadrants_num, stripes, stripes_flag)
        if peripheral_flag and stripes_flag:
            density_per_quadrant = self.get_peripheral_quadrants_and_slices_densities(quadrants_num, peripheral_flag,
                                                                                      stripes, stripes_flag)
        return density_per_quadrant

    def compute_quadrant_densities(self, quadrants_num=4, peripheral_flag=False, stripes=1, stripes_flag=False) -> np.ndarray:
        """
        For all possible subdivisions of the cell in quadrants (90 possible) and slice (if relevant)
        computes the normalized density (vs whole cytoplasm) per quadrant
        and keeps the subdivision such that the MTOC containing quadrant is the densiest.
        The anchor for the computation is the MTOC containing quadrant.
        Returns an array with the density values per quadrant (and slice) and associated MTOC flags
        """
        if not quadrants_num in [2, 3, 4, 5, 6, 8, 9]:  # just in case
            raise (RuntimeError, "Unexpected number of quadrants %i" % quadrants_num)

        max_density = 0.0
        quadrants_max_MTOC_density = np.zeros((quadrants_num*stripes, 2))
        mtoc_position = self.get_mtoc_position()

        degree_span = 360 // quadrants_num
        for degree in range(degree_span):
            quadrant_mask = self.rotate_quadrant_mask(degree, quadrants_num)
            mtoc_quad_num = quadrant_mask[mtoc_position[1], mtoc_position[0]]
            # assign each spot to the corresponding quadrant excluding those in the nucleus
            if (not peripheral_flag) and (not stripes_flag):
                density_per_quadrant = self.compute_density_per_quadrant(mtoc_quad_num, quadrant_mask, quadrants_num)
            if peripheral_flag and (not stripes_flag):
                density_per_quadrant = self.compute_peripheral_density_per_quadrant(mtoc_quad_num, quadrant_mask,
                                                                                    quadrants_num)
            if (not peripheral_flag) and stripes_flag:
                density_per_quadrant = self.compute_density_per_quadrant_and_slices(mtoc_quad_num, quadrant_mask,
                                                                                    stripes, quadrants_num)
            if peripheral_flag and stripes_flag:
                density_per_quadrant = self.compute_peripheral_density_per_quadrant_and_slices(mtoc_quad_num, quadrant_mask,
                                                                                               stripes, quadrants_num)
            mtoc_density = density_per_quadrant[density_per_quadrant[:,1] == 1][:,0].sum()
            if mtoc_density > max_density:
                max_density = mtoc_density
                quadrants_max_MTOC_density = density_per_quadrant

        assert quadrants_max_MTOC_density.shape[0] == quadrants_num * stripes, "Density array of wrong shape"
        return quadrants_max_MTOC_density

    def compute_density_per_quadrant(self, mtoc_quad_num, quadrant_mask, quadrants_num):
        raise NotImplementedError

    def compute_peripheral_density_per_quadrant(self, mtoc_quad_num, quadrant_mask, quadrants_num):
        raise NotImplementedError

    def compute_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes, quadrants_num):
        raise NotImplementedError

    def compute_peripheral_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes, quadrants_num):
        raise NotImplementedError

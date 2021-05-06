#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math
import numpy as np
import tqdm

import helpers
# secondary basic image descriptors
from constants import NUCLEUS_CENTROID_PATH_SUFFIX
from constants import ZLINES_PATH_SUFFIX
from constants import Z_LINE_DISTANCE_PATH_SUFFIX
from repository import Repository

from image import Image
from imageWithSpots import ImageWithSpots, ImageWithSpotsAndMTOC
from imageWithIntensities import ImageWithIntensities, ImageWithIntensitiesAndMTOC

class ImageWithSpotsAndIntensities(ImageWithSpots, ImageWithIntensities):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithSpots.is_a(repo, path) and ImageWithIntensities.is_a(repo, path)


class ImageWithSpotsAndIntensitiesAndMTOC(ImageWithSpotsAndMTOC, ImageWithIntensitiesAndMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithSpotsAndMTOC.is_a(repo, path) and ImageWithIntensitiesAndMTOC.is_a(repo, path)


class imageWithSpotsAndZlines(ImageWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_include(path, path + ZLINES_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(imageWithSpotsAndZlines, self).__init__(repository, image_path)
        if not self._repository.is_include(image_path, image_path + ZLINES_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_z_lines_masks(self):
        descriptor = self._path + ZLINES_PATH_SUFFIX
        if not self._repository.is_include(self._path, descriptor):
            raise LookupError("No zlines for image %s" % self._path)
        raw_value = self._repository.get_multiple(self._path, descriptor)
        z_lines = []
        for zmask in raw_value:
            z_linemask = np.array(zmask)
            z_lines.append(z_linemask)
        return z_lines

    @helpers.checkpoint_decorator(Z_LINE_DISTANCE_PATH_SUFFIX, float)
    def get_minimal_z_line_distance(self, z_line_spacing):
        return self.compute_minimal_z_line_distance(z_line_spacing)

    def compute_minimal_z_line_distance(self, z_line_spacing):
        spots = self.get_spots()
        z_lines = self.get_z_lines_masks()
        image_spots_min_distance = np.zeros(len(spots))
        z_lines_idx = helpers.reduce_z_line_mask(z_lines, spots)
        spots_reduced = spots[z_lines_idx[0] <= spots[:, 2]]
        spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
        spots_count = 0
        z_line_distance_profile = np.zeros(z_line_spacing)
        for spot in tqdm.tqdm(spots_reduced, desc="Spots"):
            z_line_mask = z_lines[int(spot[2])]
            if z_line_mask[spot[1], spot[0]] == 0:
                total_segment = np.zeros((360, z_line_spacing))
                for degree in range(360):
                    z_line_segment = np.zeros(z_line_spacing)
                    line = np.zeros((z_line_spacing, 2))
                    angle = degree * 2 * math.pi / 360
                    x_slope, y_slope = math.sin(angle), math.cos(angle)
                    for point in range(z_line_spacing):
                        x = int(round(spot[0] + point * x_slope))
                        y = int(round(spot[1] + point * y_slope))
                        line[point, 0] = x
                        line[point, 1] = y
                        if (x >= 0 and x < z_line_mask.shape[1] and y >= 0 and y < z_line_mask.shape[0]):
                            z_line_segment[point] = z_line_mask[y, x]
                            total_segment[degree, point] = z_line_mask[y, x]
                        else:
                            z_line_segment[point] = 0
                            total_segment[degree, point] = 0

                distance = helpers.compute_minimal_distance(np.sum(total_segment, axis=0))
                image_spots_min_distance[spots_count] = distance
            else:
                image_spots_min_distance[spots_count] = 0

            spots_count += 1

        for i in range(z_line_spacing):
            z_line_distance_profile[i] = float(len(np.where(image_spots_min_distance == i)[0])) / float(
                len(image_spots_min_distance))

        return z_line_distance_profile


class ImageMultiNucleus(Image):
    """ Represents a generic image with one cell mask and multiple nucleus. It has at least a cell_mask, a nucleus_mask and one or more nucleus_centroid """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_include(path, path + NUCLEUS_CENTROID_PATH_SUFFIX) and repo.is_multiple(path,
                                                                                               path + NUCLEUS_CENTROID_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        """
        Sets instance variable path to the image location (image_path) in an HDF5 file
        """
        self._repository = repository

        if not self._repository.is_present(image_path):
            raise LookupError("Cannot init, image %s is absent from repository" % image_path)
        if not self._repository.is_include(image_path, image_path + NUCLEUS_CENTROID_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image - no nucleus %s" % image_path)
        if not self._repository.is_multiple(image_path, image_path + NUCLEUS_CENTROID_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image - no multiple nucleus %s" % image_path)
        self._path = image_path

    def get_nucleus_centroid(self) -> np.ndarray:
        return np.around(self._repository.get(self._path + NUCLEUS_CENTROID_PATH_SUFFIX)).flatten().astype(np.int)

    def get_multiple_nucleus_centroid(self) -> np.ndarray:
        return [np.around(val).flatten().astype(np.int) for val in
                self._repository.get_multiple(self._path, self._path + NUCLEUS_CENTROID_PATH_SUFFIX)]

    def compute_nucleus_area(self):
        """compute nucleus surface in real pixel size using nucleus mask"""
        nucleus_mask = self.get_nucleus_mask()
        area = nucleus_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return area / len(self.get_multiple_nucleus_centroid())
        else:
            return area

    def compute_cell_area(self):
        """compute cell surface in real pixel size using cell mask"""
        cell_mask = self.get_cell_mask()
        area = cell_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return area / len(self.get_multiple_nucleus_centroid())
        else:
            return area


class ImageMultiNucleusWithSpots(ImageMultiNucleus, ImageWithSpots):
    """ Represents an image with identified spots (e.g. from FISH), has to have spots descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageMultiNucleus.is_a(repo, path) and ImageWithSpots.is_a(repo, path)


class imageMultiNucleusWithSpotsAndZlines(ImageMultiNucleusWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_include(path, path + ZLINES_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(imageMultiNucleusWithSpotsAndZlines, self).__init__(repository, image_path)
        if not self._repository.is_include(image_path, image_path + ZLINES_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_z_lines_masks(self):
        assert False, "This function is not tested"
        descriptor = self._path + ZLINES_PATH_SUFFIX
        if not self._repository.is_include(self._path, descriptor):
            raise LookupError("No zlines for image %s" % self._path)
        raw_value = self._repository.get_multiple(self._path, descriptor)
        z_lines = []
        for zmask in raw_value:
            z_linemask = np.array(zmask)
            z_lines.append(z_linemask)
        return z_lines

    @helpers.checkpoint_decorator(Z_LINE_DISTANCE_PATH_SUFFIX, float)
    def get_minimal_z_line_distance(self, z_line_spacing):
        return self.compute_minimal_z_line_distance(z_line_spacing)

    def compute_minimal_z_line_distance(self, z_line_spacing):
        assert False, "This function is not tested"
        spots = self.get_spots()
        z_lines = self.get_z_lines_masks()
        image_spots_min_distance = np.zeros(len(spots))
        z_lines_idx = helpers.reduce_z_line_mask(z_lines, spots)
        spots_reduced = spots[z_lines_idx[0] <= spots[:, 2]]
        spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
        spots_count = 0
        z_line_distance_profile = np.zeros(z_line_spacing)
        for spot in tqdm.tqdm(spots_reduced, desc="Spots"):
            z_line_mask = z_lines[int(spot[2])]
            if z_line_mask[spot[1], spot[0]] == 0:
                total_segment = np.zeros((360, z_line_spacing))
                for degree in range(360):
                    z_line_segment = np.zeros(z_line_spacing)
                    line = np.zeros((z_line_spacing, 2))
                    angle = degree * 2 * math.pi / 360
                    x_slope, y_slope = math.sin(angle), math.cos(angle)
                    for point in range(z_line_spacing):
                        x = int(round(spot[0] + point * x_slope))
                        y = int(round(spot[1] + point * y_slope))
                        line[point, 0] = x
                        line[point, 1] = y
                        if (x >= 0 and x < z_line_mask.shape[1] and y >= 0 and y < z_line_mask.shape[0]):
                            z_line_segment[point] = z_line_mask[y, x]
                            total_segment[degree, point] = z_line_mask[y, x]
                        else:
                            z_line_segment[point] = 0
                            total_segment[degree, point] = 0

                distance = helpers.compute_minimal_distance(np.sum(total_segment, axis=0))
                image_spots_min_distance[spots_count] = distance
            else:
                image_spots_min_distance[spots_count] = 0

            spots_count += 1

        for i in range(z_line_spacing):
            z_line_distance_profile[i] = float(len(np.where(image_spots_min_distance == i)[0])) / float(
                len(image_spots_min_distance))

        return z_line_distance_profile
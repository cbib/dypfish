#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math

import numpy as np
from loguru import logger

import constants
import helpers
from constants import CELL_MASK_SLICES_PATH_SUFFIX
from constants import HEIGHT_MAP_PATH_SUFFIX
from constants import NUCLEUS_CENTROID_PATH_SUFFIX
from constants import VOLUMES_FROM_PERIPHERY_PATH_SUFFIX
from constants import ZERO_LEVEL_PATH_SUFFIX
from image import Image, ImageWithMTOC
from repository import Repository


class Image3d(Image):
    """ Represents an 3D image, has to have a height map descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        # TODO : check if we need zero level to define a 3D cell.
        # comment ZERO_LEVEL_PATH_SUFFIX is present() for Cultured data to run volume corrected analysis.
        # they do not have ZERO LEVEL descriptors and so was rejected as a 3D image
        return repo.is_present(path + HEIGHT_MAP_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(Image3d, self).__init__(repository, image_path)
        if not self._repository.is_present(image_path + HEIGHT_MAP_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_height_map(self) -> np.ndarray:
        '''This function should not be called by other classes / analyses
         The result is not restricted to the cell_mask and not adjusted to zero_level
         Aways call adjust_height_map instead'''
        descriptor = self._path + HEIGHT_MAP_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No height map for image %s" % self._path)
        return np.array(self._repository.get(descriptor))

    def get_cytoplasm_height_map(self):
        height_map = self.adjust_height_map()
        height_map[self.get_cytoplasm_mask() == 0] = 0
        return height_map

    def adjust_height_map(self, cytoplasm=False):
        '''
        adjust the height_map to respect the zero_level;
        for the coherency sake the periphery of the cell is set
        to be at 0.5 height if the cell_map is wider than the lowest height_map level.
        cytoplasm = True is the same as get_cytoplasm_height_map
        '''
        if cytoplasm == True:
            height_map = self.get_cytoplasm_height_map().astype(float)
            mask = self.get_cytoplasm_mask()
        else:
            height_map = self.get_height_map().astype(float)
            mask = self.get_cell_mask()

        height_diff = int(np.max(height_map)) - self.get_zero_level()
        height_map[np.where(height_map > self.get_zero_level())] -= height_diff # TODO possibly zero_level+1
        height_map[(mask == 1) & (height_map == 0)] = 0.5
        return height_map

    def get_zero_level(self):
        descriptor = self._path + ZERO_LEVEL_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No zero level for image %s" % self._path)
        return np.array(self._repository.get(descriptor))

    @helpers.checkpoint_decorator(CELL_MASK_SLICES_PATH_SUFFIX, dtype=np.float)
    def get_cell_mask_slices(self) :
        return self.compute_cell_mask_slices()

    def compute_cell_mask_slices(self, height_map=None, zero_level=None,
                                 image_width=None, image_height=None) -> np.ndarray:
        """
        Reconstructs the z-slices given a height_map;
        out of focus slices (defined by zero_level) are reconstructed
        index 0 is for the bottom slice, index zero_level is for the top;
        this makes it coherent with the height_map, however z spots coordinates
        are in the inversed relationship with the corresponding slice index
        """
        if height_map is None: height_map = self.adjust_height_map()
        #zero_level = zero_level or self.get_zero_level()
        max_height = np.max(height_map).astype(int)
        image_width = image_width or constants.dataset_config["IMAGE_WIDTH"]
        image_height = image_height or constants.dataset_config["IMAGE_HEIGHT"]

        # Create binary cell masks per slice
        slice_masks = np.zeros((image_width, image_height, max_height+1)) #zero_level+1))
        for slice_num in range(0, slice_masks.shape[2]): # we take the bottom slice but it is adjusted
            slice_mask = np.array(height_map, copy=True)
            slice_mask[height_map < slice_num+1] = 0
            slice_mask[height_map >= slice_num+0.5] = 1 # the 0.5 is because of the adjustment
            slice_masks[:, :, slice_num] = slice_mask # 0 slice the largest
        return slice_masks

    def get_peripheral_cell_volume(self, peripheral_threshold=None):
        peripheral_threshold = peripheral_threshold or constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        volumes = self.get_volumes_from_periphery()
        return volumes[peripheral_threshold]

    def compute_peripheral_cell_volume(self, peripheral_threshold):
        """
        compute volume in the periphery in real pixel size using the cell mask
        """
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = ((cell_mask_dist_map > 0) &
                                  (cell_mask_dist_map <= peripheral_threshold)).astype(int)
        height_map = self.adjust_height_map()
        height_map_periph = np.multiply(height_map, peripheral_binary_mask)
        peripheral_cell_volume = height_map_periph.sum() * helpers.volume_coeff()
        return peripheral_cell_volume

    @helpers.checkpoint_decorator(VOLUMES_FROM_PERIPHERY_PATH_SUFFIX, dtype=np.float)
    def get_volumes_from_periphery(self):
        return self.compute_volumes_from_periphery()

    def compute_volumes_from_periphery(self):
        logger.info("Computing {} volumes from the periphery for {}",
                    constants.analysis_config['NUM_CONTOURS'], self._path)
        volumes = np.zeros(constants.analysis_config['NUM_CONTOURS'])
        for i in range(0, constants.analysis_config['NUM_CONTOURS']):
            volumes[i] = self.compute_peripheral_cell_volume(peripheral_threshold=i+1)
        return volumes

    def compute_cell_volume(self):
        """
        compute cell volume in real pixel size using the cell mask
        """
        cell_mask = self.get_cell_mask()
        height_map = self.adjust_height_map()
        cell_volume = height_map[cell_mask == 1].sum()
        return cell_volume * helpers.volume_coeff()

    def compute_nucleus_volume(self):
        """compute nucleus volume in in real pixel size using nucleus mask"""
        nucleus_mask = self.get_nucleus_mask()
        height_map = self.adjust_height_map()
        nucleus_volume = height_map[nucleus_mask == 1].sum()
        return nucleus_volume * helpers.volume_coeff()

    def compute_cytoplasmic_volume(self):
        """compute volume of the cytoplasm in in real pixel size using cytoplasm mask and adjusted height map"""
        cytoplasm_mask = self.get_cytoplasm_mask()
        height_map = self.adjust_height_map()
        cytoplasm_volume = height_map[cytoplasm_mask == 1].sum()
        return cytoplasm_volume * helpers.volume_coeff()

    def compute_peripheral_volume(self):
         peripheral_cell_volume = self.compute_peripheral_cell_volume()
         if (peripheral_cell_volume <= 0):
             raise (RuntimeError, "peripheral area has inconsistent volumes for image %s" % self._path)
         return peripheral_cell_volume

    def compute_median_cytoplasmic_distance_from_nucleus3d(self, dsAll) -> float:
        '''
        Computes median and maximal distance of a cytoplasmic voxel from the nucleus centroid
        '''
        height_map = self.get_cytoplasm_height_map()
        cytoplasm_mask = self.get_cytoplasm_mask()
        distances = np.array([])
        slices = self.get_cell_mask_slices()
        for slice_num in range(0, slices.shape[2]):
            if slice_num >= dsAll.shape[2]: continue # should not happen, just in case
            slice_mask = slices[:,:,slice_num]
            slice_mask = slice_mask * cytoplasm_mask # not very precise since it is not propagarted through slices
            slice_distances = np.multiply(dsAll[:,:,slice_num], slice_mask)
            distances = np.append(distances, slice_distances[slice_distances > 0])

        return np.median(distances), np.max(distances)

    def compute_cytoplasmic_coordinates_peripheral_distance(self, coordinates) -> np.ndarray:
        """
         Perform the computation of distance to the periphery for cytoplasmic coordinates
         Computation is done in pseudo 3D
         Returns an array of integer distances
         """
        logger.info("Computing 3D peripheral distance for {} coordinates in image {}", len(coordinates), self._path)
        assert coordinates.shape[1] == 3, "3D coordinates needed for distance to the periphery"
        max_height = np.max(self.adjust_height_map())
        peripheral_distance_map = self.get_cell_mask_distance_map()
        cell_area = self.compute_cell_area()

        distances = np.array([])
        slices = self.get_cell_mask_slices()
        for slice_num in range(0, slices.shape[2]):
            slice_mask = slices[:, :, slice_num]
            slice_area = slice_mask.sum() * helpers.surface_coeff()
            assert slice_area >= 0, "Empty slice area"
            coordinates_this_height = coordinates[np.around(coordinates[:, 2]) == max_height - slice_num]
            for c in coordinates_this_height:
                if not self.is_in_cytoplasm(c[0:2][::-1]):
                    distances = np.append(distances, np.nan)
                else:
                    area_ratio = slice_area / cell_area
                    old_periph_distance = peripheral_distance_map[c[1], c[0]]
                    new_periph_distance = math.sqrt(area_ratio) * old_periph_distance
                    if (new_periph_distance > old_periph_distance): # coordinate falls out of the slice
                        distances = np.append(distances, np.nan)
                    else:
                        distances = np.append(distances, int(np.around(new_periph_distance)))

        if (len(distances) < len(coordinates)) or (len(distances[np.isnan(distances)]) > 0):
            num = len(coordinates) - len(distances) + len(distances[np.isnan(distances)])
            logger.info("  found {} coordinates outside of cytoplasm", num)
        return distances


class Image3dWithMTOC(Image3d, ImageWithMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithMTOC.is_a(repo, path) and Image3d.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4):
        raise NotImplementedError

    def compute_peripheral_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4):
        """
        compute volumic density per quadrant;
        return values of density paired with the MTOC presence flag (0/1)
        Note that this is not relative density
        """
        peripheral_fraction_threshold = constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"]
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        quadrant_mask = quadrant_mask * peripheral_binary_mask
        return self.compute_density_per_quadrant(mtoc_quad, quadrant_mask, quadrants_num)

    def compute_density_per_quadrant_and_slices(self, mtoc_quad, quadrant_mask, stripes, quadrants_num=4,
                                                peripheral_flag=False):
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        slices_per_stripe = np.floor(100.0 / stripes)  # number of isolines per stripe
        if peripheral_flag:
            slices_per_stripe = np.floor(constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"] / stripes)
        arr = np.empty((0, 2), float)
        for stripe_num in range(1, stripes + 1):
            stripe_mask = (cell_mask_dist_map > (stripe_num - 1) * slices_per_stripe) & \
                          (cell_mask_dist_map <= stripe_num * slices_per_stripe + 1).astype(int)
            stripe_quadrant_mask = quadrant_mask * stripe_mask
            res = self.compute_density_per_quadrant(mtoc_quad, stripe_quadrant_mask, quadrants_num)
            arr = np.append(arr, res[res[:, 1].argsort()[::-1]], axis=0)  # MTOC quadrant slice always first
        assert arr.shape[0] == quadrants_num * stripes, "Incorrect shape"
        return arr

    def compute_peripheral_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes,
                                                           quadrants_num=4):
        peripheral_fraction_threshold = constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"]
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        quadrant_mask = quadrant_mask * peripheral_binary_mask
        return self.compute_density_per_quadrant_and_slices(mtoc_quad_num, quadrant_mask, stripes,
                                                            quadrants_num, peripheral_flag=True)


class Image3dMultiNucleus(Image3d):
    """
    Represents a generic image with one cell mask and multiple nucleus.
    It has at least a cell_mask, a nucleus_mask and one or more nucleus_centroid
    """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_multiple(path, path + NUCLEUS_CENTROID_PATH_SUFFIX) and Image3d.is_a(repo, path)

    def __init__(self, repository: Repository, image_path: str):
        super(Image3dMultiNucleus, self).__init__(repository, image_path)
        if not self._repository.is_multiple(image_path, image_path + NUCLEUS_CENTROID_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def compute_cell_volume(self):
        """
        computes cell volume in in real pixel size using the cell mask
        """
        cell_mask = self.get_cell_mask()
        height_map = self.adjust_height_map()
        cell_volume = height_map[np.where(cell_mask[:] == 1)].sum()
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return cell_volume / len(self.get_multiple_nucleus_centroid()) * helpers.volume_coeff()
        else:
            return cell_volume * helpers.volume_coeff()

    def compute_peripheral_cell_volume(self):
        """
        computes cell volume in in real pixel size using the cell mask
        """
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        cell_mask = self.get_cell_mask()
        height_map = self.adjust_height_map()
        height_map_periph = np.multiply(height_map, peripheral_binary_mask)
        peripheral_cell_volume = height_map_periph[np.where(cell_mask[:] == 1)].sum()

        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return peripheral_cell_volume / len(self.get_multiple_nucleus_centroid()) * helpers.volume_coeff()
        else:
            return peripheral_cell_volume * helpers.volume_coeff()


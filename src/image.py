#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math
from typing import List
import numpy as np
import numexpr
import tqdm
from loguru import logger
from scipy import signal
import constants
import helpers
import image_processing as ip
from repository import Repository
from sklearn.metrics.pairwise import pairwise_distances

from constants import SPOTS_PATH_SUFFIX
from constants import ZLINES_PATH_SUFFIX
from constants import Z_LINE_DISTANCE_PATH_SUFFIX
from constants import IF_PATH_SUFFIX
from constants import CELL_MASK_PATH_SUFFIX
from constants import PERIPHERAL_MASK_PATH_SUFFIX
from constants import NUCLEUS_MASK_PATH_SUFFIX
from constants import CYTOPLASM_MASK_PATH_SUFFIX
from constants import MTOC_POSITION_PATH_SUFFIX
from constants import MTOC_LEADING_EDGE_SUFFIX

# secondary basic image descriptors
from constants import NUCLEUS_AREA_PATH_SUFFIX
from constants import CELL_AREA_PATH_SUFFIX
from constants import CELL_MASK_DISTANCE_PATH_SUFFIX
from constants import SPOTS_PERIPHERAL_DISTANCE_2D_PATH_SUFFIX
from constants import CYTOPLASMIC_SPOTS_PATH_SUFFIX
from constants import PERIPHERAL_SPOTS_PATH_SUFFIX
from constants import CYTOPLASMIC_SPOT_COUNT_PATH_SUFFIX
from constants import PERIPHERAL_SPOT_COUNT_PATH_SUFFIX
from constants import TOTAL_INTENSITY_PATH_SUFFIX
from constants import CELL_TOTAL_INTENSITY_PATH_SUFFIX
from constants import CYTOPLASMIC_TOTAL_INTENSITY_PATH_SUFFIX
from constants import PERIPHERAL_TOTAL_INTENSITY_PATH_SUFFIX
from constants import CYTOPLASMIC_INTENSITIES_PATH_SUFFIX
from constants import PERIPHERAL_INTENSITIES_PATH_SUFFIX
from constants import CLUSTERING_INDICES_PATH_SUFFIX
from constants import QUADRANT_DENSITIES_PATH_SUFFIX
from constants import PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX
from constants import QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX
from constants import NUCLEUS_CENTROID_PATH_SUFFIX
from constants import PERIPHERAL_QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX

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

    def get_peripheral_cell_mask(self) -> np.ndarray:
        cell_mask = self.get_cell_mask()
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        return np.multiply(cell_mask, peripheral_binary_mask)

    def get_nucleus_mask(self) -> np.ndarray:
        return np.array(self._repository.get(self._path + NUCLEUS_MASK_PATH_SUFFIX)).astype(int)

    def compute_cytoplasm_mask(self) -> np.ndarray:
        cell_mask = self.get_cell_mask()
        nucleus_mask = 1 - self.get_nucleus_mask()
        return np.multiply(cell_mask, nucleus_mask)

    def compute_peripheral_mask(self) -> np.ndarray:
        peripheral_cell_mask = self.get_peripheral_cell_mask()
        nucleus_mask = 1 - self.get_nucleus_mask()
        return np.multiply(peripheral_cell_mask, nucleus_mask)

    def get_nucleus_centroid(self) -> np.ndarray:
        return np.around(
            self._repository.get_by_regex(self._path, self._path + NUCLEUS_CENTROID_PATH_SUFFIX)).flatten().astype(
            np.int)

    def compute_nucleus_area(self):
        """compute nucleus surface in pixel using nucleus mask"""
        nucleus_mask = self.get_nucleus_mask()
        area = nucleus_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        return area

    @helpers.checkpoint_decorator(NUCLEUS_AREA_PATH_SUFFIX, float)
    def get_nucleus_area(self) -> float:
        return self.compute_nucleus_area()

    def compute_cell_area(self):
        """compute cell surface in pixel using cell mask"""
        cell_mask = self.get_cell_mask()
        area = cell_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        return area

    @helpers.checkpoint_decorator(CELL_AREA_PATH_SUFFIX, float)
    def get_cell_area(self) -> float:
        return self.compute_cell_area()

    @helpers.checkpoint_decorator(CYTOPLASM_MASK_PATH_SUFFIX, float)
    def get_cytoplasm_mask(self) -> np.ndarray:
        return self.compute_cytoplasm_mask()

    @helpers.checkpoint_decorator(PERIPHERAL_MASK_PATH_SUFFIX, float)
    def get_peripheral_mask(self) -> np.ndarray:
        return self.compute_peripheral_mask()

    def compute_cell_mask_distance_map(self):
        cell_mask = self.get_cell_mask()
        nucleus_mask = self.get_nucleus_mask()
        nucleus_centroid = self.get_nucleus_centroid().transpose()  # TODO why the transpose ??
        cytoplasm_mask = (cell_mask == 1) & (nucleus_mask == 0)

        # for each degree, analyse the line segment between the nucleus and the periphery
        contour_points = ip.compute_contour_points(nucleus_mask, nucleus_centroid, cytoplasm_mask,
                                                   num_contours=constants.analysis_config['NUM_CONTOURS'])
        cell_mask_distance_map = ip.compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points)
        return cell_mask_distance_map

    @helpers.checkpoint_decorator(CELL_MASK_DISTANCE_PATH_SUFFIX, dtype=np.int)
    def get_cell_mask_distance_map(self):
        return self.compute_cell_mask_distance_map()

    def is_in_cytoplasm(self, position: List) -> bool:
        cytoplasm_mask = self.get_cytoplasm_mask()
        return cytoplasm_mask[position[0], position[1]] == 1

    def compute_peripheral_areas(self) -> List[float]:
        """
        TODO cache this ?
        TODO : why the returned length is 101 and not 100 ?
        :return:
        """
        peripheral_areas = []
        peripheral_distance_map = self.get_cell_mask_distance_map()
        cytoplasm_mask = self.get_cytoplasm_mask()
        peripheral_distance_map[(peripheral_distance_map == 0) & (cytoplasm_mask == 1)] = 1
        nucleus_area = self.get_nucleus_area()
        for i in range(101):
            tmp_mask = np.array(peripheral_distance_map, copy=True)
            tmp_mask[tmp_mask <= i] = 0
            tmp_mask[(tmp_mask > i) & (tmp_mask <= 100)] = 1
            peripheral_areas.append(tmp_mask.sum() * helpers.surface_coeff() + nucleus_area)
        return peripheral_areas


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

    def compute_quadrant_mask(self, degree, slices_num=4, nucleus_centroid=None, mtoc_position=None,
                              cell_mask=None, image_width=None, image_height=None):
        """
        Computes the quadrant mask (slices are quadrants if slices_num == 4) of a cell anchored in
        the MTOC position, rotated by a degree. The original quadrant of the MTOC (before rotation)
        is defined by two lines 45 degrees to the right
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
        # This trunc of pi is used for compatibility between macOS and Linux on how they deal with infinite decimal (0.9999999999999...)
        pi = format(math.pi, '.10f')
        pi = float(pi)

        sliceno = ((pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / ((8 / slices_num) * pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == slices_num + 1] = slices_num  # int conversion sometimes rounds the value
        quadrant_mask[cell_mask == 0] = 0
        return quadrant_mask

    @helpers.checkpoint_decorator(PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_quadrants_densities(self, quadrants_num=4):
        return self.peripheral_split_in_quadrants(quadrants_num=quadrants_num)

    def peripheral_split_in_quadrants(self, quadrants_num=4) -> np.ndarray:
        """
        was : compute_max_density_MTOC_quadrant
        For all possible subdivisions of the cell in quadrants (90 possible)
        computes the normalized density (vs whole cytoplasm) per quadrant
        and keeps the subdivision such that the MTOC containing quadrant is the densiest.
        The anchor for the computation is the MTOC containing quadrant.
        Returns an array with the density values per quadrant and associated MTOC flags
        """
        if not quadrants_num in [2, 3, 4, 5, 6, 8, 9]:  # just in case
            raise (RuntimeError, "Unexpected number of slices (quadrants) %i" % quadrants_num)

        max_density = 0.0
        quadrants_max_MTOC_density = np.zeros((quadrants_num, 2), dtype=float)
        mtoc_position = self.get_mtoc_position()

        degree_span = 360 // quadrants_num
        for degree in range(degree_span):
            quadrant_mask = self.compute_quadrant_mask(degree, quadrants_num)
            mtoc_quad_num = quadrant_mask[mtoc_position[1], mtoc_position[0]]
            # assign each spot to the corresponding quadrant excluding those in the nucleus
            density_per_quadrant = self.compute_density_per_quadrant(mtoc_quad_num, quadrant_mask, quadrants_num)
            if density_per_quadrant[mtoc_quad_num - 1, 0] > max_density:
                max_density = density_per_quadrant[mtoc_quad_num - 1, 0]
                quadrants_max_MTOC_density = density_per_quadrant

        return quadrants_max_MTOC_density

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
            density_per_quadrant = self.get_quadrants_and_slices_densities(quadrants_num, stripes,stripes_flag)
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
            quadrant_mask = self.compute_quadrant_mask(degree, quadrants_num)
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

    def compute_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes, quadrants_num):
        raise NotImplementedError

    def compute_peripheral_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes, quadrants_num):
        raise NotImplementedError

    def compute_peripheral_density_per_quadrant(self, mtoc_quad_num, quadrant_mask, quadrants_num):
        raise NotImplementedError


class ImageWithSpots(Image):
    """ Represents an image with identified spots (e.g. from FISH), has to have spots descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_present(path + SPOTS_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(ImageWithSpots, self).__init__(repository, image_path)
        if not self._repository.is_present(image_path + SPOTS_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_spots(self) -> np.ndarray:
        """
        :return: np.ndarray of int
        """
        descriptor = self._path + SPOTS_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No spots for image %s" % self._path)
        raw_value = self._repository.get(descriptor)
        spots = np.around(np.array(raw_value)).astype(np.int)
        return spots

    # TODO : this will load the mask for each spot
    def compute_cytoplasmic_spots(self) -> np.ndarray:
        """
        :return: np.ndarray of int
        """
        spots = self.get_spots()
        mask = [self.is_in_cytoplasm(s[0:2][::-1]) for s in spots]  # TODO check coordinate coherency for spots
        cytoplasmic_spots = spots[mask]
        return np.asarray(cytoplasmic_spots, dtype=int)

    # TODO : this will load the mask for each spot
    def compute_peripheral_spots(self) -> np.ndarray:
        """
        :return: np.ndarray of int
        """
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        spots = self.get_spots()
        mask = [False if peripheral_binary_mask[s[1], s[0]] == 0 else True for s in
                spots]  # TODO check coordinate coherency for spots
        peripheral_spots = spots[mask]
        return np.asarray(peripheral_spots, dtype=int)

    @helpers.checkpoint_decorator(CYTOPLASMIC_SPOTS_PATH_SUFFIX, np.int)
    def get_cytoplasmic_spots(self) -> np.ndarray:
        return self.compute_cytoplasmic_spots()

    @helpers.checkpoint_decorator(PERIPHERAL_SPOTS_PATH_SUFFIX, np.int)
    def get_peripheral_spots(self) -> np.ndarray:
        return self.compute_peripheral_spots()

    def compute_spots_peripheral_distance_2D(self) -> np.ndarray:
        """
        Return an array of distances to the periphery for cytoplasmic spots (np.ndarray of uint8)
        """
        spots = self.get_cytoplasmic_spots()
        logger.info("Computing {} spots 2D peripheral distance spots for {}",
                    len(spots), self._path)

        peripheral_distance_map = self.get_cell_mask_distance_map()
        spots_distances = [peripheral_distance_map[s[1], s[0]] for s in spots]
        spots_distances = [1 if d == 0 else d for d in spots_distances]  # TODO : why this hack ?
        return np.asarray(spots_distances, dtype=np.uint8)

    @helpers.checkpoint_decorator(SPOTS_PERIPHERAL_DISTANCE_2D_PATH_SUFFIX, np.uint8)
    def get_spots_peripheral_distance(self):
        return self.compute_spots_peripheral_distance_2D()

    def compute_cytoplasmic_total_spots(self):
        cytoplasmic_spots = self.get_cytoplasmic_spots()
        return len(cytoplasmic_spots)

    @helpers.checkpoint_decorator(CYTOPLASMIC_SPOT_COUNT_PATH_SUFFIX, dtype=np.uint)
    def get_cytoplasmic_total_spots(self):
        return self.compute_cytoplasmic_total_spots()

    def compute_peripheral_total_spots(self):
        peripheral_spots = self.get_peripheral_spots()
        return len(peripheral_spots)

    @helpers.checkpoint_decorator(PERIPHERAL_SPOT_COUNT_PATH_SUFFIX, dtype=np.uint)
    def get_peripheral_total_spots(self):
        return self.compute_peripheral_total_spots()

    def compute_average_cytoplasmic_distance_from_nucleus(self, dsAll) -> float:
        cytoplasm_mask = self.get_cytoplasm_mask()
        map_dist = np.multiply(cytoplasm_mask, dsAll)  # TODO : this should be useless
        return map_dist.sum() / cytoplasm_mask.sum()

    def compute_spots_normalized_distance_to_centroid(self) -> float:
        nucleus_centroid = self.get_nucleus_centroid()
        cytoplasmic_spots = self.get_cytoplasmic_spots()[:, 0:2]  # 2d coordinates of 3d spots
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid)  # 2d distances from nucleus_centroid

        dsCytoplasmic = dsAll[cytoplasmic_spots[:, 1], cytoplasmic_spots[:, 0]]
        S = self.compute_average_cytoplasmic_distance_from_nucleus(dsAll)

        # val is average 2D distance from the nucleus centroid of cytoplasmic mRNAs
        # normalized by the cytoplasmic cell spread (taking a value 1 when mRNAs are evenly
        # distributed across the cytoplasm).
        # TODO : this assumption is probably wrong. Should check it
        normalized_average_2d_distance = np.mean(dsCytoplasmic) / S
        return normalized_average_2d_distance

    def compute_spots_normalized_cytoplasmic_spread(self) -> float:
        '''
        Computes the Standard Distance measure of spread as the average distance for all spots
        from the Mean Center. This measures the compactness of a distribution of spots.
        In a Normal Distribution you would expect around 68% of all points to fall within
        the Standard Distance.
        '''
        cytoplasmic_spots = self.get_cytoplasmic_spots()[:, 0:2]
        mu_x = cytoplasmic_spots[:, 0].sum() / len(cytoplasmic_spots)
        mu_y = cytoplasmic_spots[:, 1].sum() / len(cytoplasmic_spots)
        sd = math.sqrt( np.sum((cytoplasmic_spots[:, 0] - mu_x) ** 2) / len(cytoplasmic_spots)+
                           np.sum((cytoplasmic_spots[:, 1] - mu_y) ** 2) / len(cytoplasmic_spots) )
        d = pairwise_distances(cytoplasmic_spots, metric='euclidean')
        return sd / np.mean(d[d != 0])

    def compute_cytoplasmic_density(self):
        # compute mRNA density in the cytoplasm
        cytoplasmic_mrna_count = self.get_cytoplasmic_total_spots()
        if cytoplasmic_mrna_count == 0:
            raise RuntimeError("Image contains no spots %s" % self._path)
        cytoplasmic_area = self.compute_cell_area()
        return cytoplasmic_mrna_count / cytoplasmic_area

    def ripley_k_point_process(self, nuw: float, my_lambda: float, spots=None, r_max: int = None) -> np.ndarray:
        if spots is None: spots = self.get_spots()
        n_spots = len(spots)
        r_max = r_max or constants.analysis_config["MAX_CELL_RADIUS"]
        K = np.zeros(r_max)
        for i in range(n_spots):
            mask = np.zeros((n_spots, 2));
            mask[i, :] = 1
            other_spots = np.ma.masked_where(mask == 1, np.ma.array(spots, mask=False)).compressed().reshape(n_spots - 1, 2)
            x_squared = np.square(spots[i, 0] - other_spots[:, 0])
            y_squared = np.square(spots[i, 1] - other_spots[:, 1])
            ds = np.sqrt(x_squared + y_squared)
            if n_spots - 1 < r_max:
                for m in range(n_spots - 1):
                    K[math.ceil(ds[m]):r_max] = K[math.ceil(ds[m]):r_max] + 1
            else:
                for m in range(r_max):
                    K[m] = K[m] + ds[ds <= m].sum()
        K = K * (1 / (my_lambda ** 2 * nuw))
        return K

    def compute_random_spots(self): # TODO : not tested
        # simulate n list of random spots
        cell_mask = self.get_cell_mask()
        n_spots = len(self.get_spots())
        x, y = np.where(cell_mask == 1)
        idx = np.random.randint(0, len(x), n_spots)  # we chose random indices
        return np.vstack((x[idx], y[idx])).T

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripkey-K computation for disks of radius r < MAX_CELL_RADIUS
        :return: clustering indices for all r
        Was : clustering_index_point_process
        """
        logger.info("Running {} simulations of Ripley-K for {}",
                    constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], self._path)
        spots = self.get_spots()
        n_spots = len(spots)
        cell_mask = self.get_cell_mask()
        nuw = (np.sum(cell_mask[:, :] == 1)) * helpers.surface_coeff()  # whole surface of the cell
        my_lambda = float(n_spots) / float(nuw)  # spot's volumic density

        k = self.ripley_k_point_process(nuw=nuw, my_lambda=my_lambda)  # TODO : first call for _all_ spots while the subsequent only for those in the height_map
        k_sim = np.zeros((constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], constants.analysis_config["MAX_CELL_RADIUS"]))

        # simulate RIPLEY_K_SIMULATION_NUMBER lists of random spots and run ripley_k
        for t in tqdm.tqdm(range(constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]), desc="Simulations"):
            random_spots = self.compute_random_spots()
            tmp_k = self.ripley_k_point_process(spots=random_spots, nuw=nuw, my_lambda=my_lambda).flatten()
            k_sim[t] = tmp_k
        h = np.subtract(np.sqrt(k / math.pi), np.arange(1, constants.analysis_config["MAX_CELL_RADIUS"] + 1).reshape((constants.analysis_config["MAX_CELL_RADIUS"], 1)))

        synth5, synth50, synth95 = helpers.compute_statistics_random_h_star_2d(h_sim=k_sim)
        return helpers.compute_h_star_2d(h, synth5, synth50, synth95)

    @helpers.checkpoint_decorator(CLUSTERING_INDICES_PATH_SUFFIX, dtype=np.float)
    def get_clustering_indices(self):
        return self.compute_clustering_indices()

    def compute_degree_of_clustering(self) -> int:
        h_star = self.get_clustering_indices()
        return np.array(h_star[h_star > 1] - 1).sum()

class ImageWithIntensities(Image):
    """ Represents an image with intensity data (e.g. from IF), has to have IF descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_present(path + IF_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(ImageWithIntensities, self).__init__(repository, image_path)
        if not self._repository.is_present(image_path + IF_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_intensities(self) -> np.ndarray:
        """
          :returns: np.ndarray of floats (intensities of the cell)
          """
        descriptor = self._path + IF_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No intensities for image %s" % self._path)
        raw_value = self._repository.get(descriptor)
        intensities = np.array(raw_value)
        return intensities

    def compute_cell_total_intensity(self) -> float:  # TODO should be redundant with compute_total_intensity, but it is not in V0
        intensities = self.get_intensities()
        cell_mask = self.get_cell_mask()
        cell_intensities = np.multiply(intensities, cell_mask)
        return cell_intensities.sum()

    def compute_cytoplasmic_intensities(self) -> np.ndarray:
        intensities = self.get_intensities()
        cytoplasm_mask = self.get_cytoplasm_mask()
        cytoplasmic_intensities = np.multiply(intensities, cytoplasm_mask)
        return cytoplasmic_intensities

    def compute_total_intensity(self) -> float:
        intensities = self.get_intensities()
        return intensities.sum()

    def compute_cytoplasmic_total_intensity(self) -> float:
        cytoplasmic_intensities = self.get_cytoplasmic_intensities()
        return cytoplasmic_intensities.sum()

    @helpers.checkpoint_decorator(CYTOPLASMIC_TOTAL_INTENSITY_PATH_SUFFIX, dtype=np.float)
    def get_cytoplasmic_total_intensity(self):
        return self.compute_cytoplasmic_total_intensity()

    def compute_peripheral_total_intensity(self) -> float:
        peripheral_intensity_mask = self.get_peripheral_intensity_mask()
        return peripheral_intensity_mask.sum()

    @helpers.checkpoint_decorator(PERIPHERAL_TOTAL_INTENSITY_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_total_intensity(self):
        return self.compute_peripheral_total_intensity()

    @helpers.checkpoint_decorator(CELL_TOTAL_INTENSITY_PATH_SUFFIX, dtype=np.float)
    def get_cell_total_intensity(self):
        return self.compute_cell_total_intensity()

    @helpers.checkpoint_decorator(TOTAL_INTENSITY_PATH_SUFFIX, dtype=np.float)
    def get_total_intensity(self):
        return self.compute_total_intensity()

    @helpers.checkpoint_decorator(CYTOPLASMIC_INTENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_cytoplasmic_intensities(self) -> np.ndarray:
        return self.compute_cytoplasmic_intensities()

    @helpers.checkpoint_decorator(PERIPHERAL_INTENSITIES_PATH_SUFFIX, dtype=np.float)
    def get_peripheral_intensity_mask(self) -> np.ndarray:
        return self.compute_peripheral_intensity_mask()


    def compute_peripheral_intensity_mask(self) -> np.ndarray:
        """
         :returns: np.ndarray of floats (total intensity for each distance percentage from the periphery)
         normalization is the responsibility of the caller
        """
        intensities = self.get_intensities()
        peripheral_mask = self.get_peripheral_mask()
        peripheral_intensity_mask = np.multiply(intensities, peripheral_mask)
        return peripheral_intensity_mask


    def compute_peripheral_intensities(self) -> np.ndarray:
        """
         :returns: np.ndarray of floats (total intensity for each distance percentage from the periphery)
         normalization is the responsibility of the caller
         """
        # cytoplasmic_intensities = self.get_cytoplasmic_intensities()
        intensities = self.get_intensities()
        cell_mask_distance_map = self.get_cell_mask_distance_map()
        intensities_sums = np.zeros(100)  # TODO check evrywhere : maybe NUM_CONTOURS is better ?
        for i in range(0, 100):  # TODO : normalization should be done by the caller
            # TODO : should be intensities_sums[i] = cytoplasmic_intensities[cell_mask_distance_map <= i + 1].sum()
            intensities_sums[i] = intensities[(cell_mask_distance_map <= i + 1) & (cell_mask_distance_map > 0)].sum()
        return intensities_sums

    # TODO : why would we consider that distance can be proportional to intensity, does not make sense
    def compute_average_cytoplasmic_distance_proportional_intensity(self, dsAll):
        cytoplasm_mask = self.get_cytoplasm_mask()
        dsCellular = np.multiply(cytoplasm_mask, dsAll)
        return dsCellular.sum() / cytoplasm_mask.sum()

    def compute_intensities_normalized_spread_to_centroid(self) -> float:
        nucleus_centroid = self.get_nucleus_centroid()
        intensities = self.get_cytoplasmic_intensities()
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid)  # 2d distances from nucleus_centroid
        dsAll = dsAll * self.get_cytoplasm_mask()
        S = self.compute_average_cytoplasmic_distance_proportional_intensity(dsAll)

        proportional_intensities = np.multiply(intensities, dsAll)
        normalized_value = proportional_intensities.sum() / (S * intensities.sum())

        return normalized_value

    def compute_intensities_normalized_cytoplasmic_spread(self):
        IF = self.get_cytoplasmic_intensities()
        cytoplasmic_mask = self.get_cytoplasm_mask()

        # Calculate the spread of signal peaks
        mean_signal = np.mean(IF[cytoplasmic_mask == 1])
        peaks = np.argwhere(IF > mean_signal * 1.5)  # arbitrary choice to reduce the number of peaks
        d = pairwise_distances(peaks, metric='euclidean')

        mu_x = peaks[:, 0].sum() / len(peaks)
        mu_y = peaks[:, 1].sum() / len(peaks)
        sd = math.sqrt(np.sum((peaks[:, 0] - mu_x) ** 2) / len(peaks) +
                       np.sum((peaks[:, 1] - mu_y) ** 2) / len(peaks))

        return sd / np.mean(d[d != 0])

    def signal_to_noise(self) -> float:
        intensities = self.get_intensities()  # the whole image
        snr = np.power(np.mean(intensities), 2) / np.power(np.std(intensities), 2)
        return snr

    def compute_cytoplasmic_density(self):
        # compute signal density of the cytoplasm
        cytoplasmic_intensity_count = self.get_cytoplasmic_total_intensity()
        cytoplasmic_area = self.compute_cell_area()
        return cytoplasmic_intensity_count / cytoplasmic_area

    @helpers.checkpoint_decorator(CLUSTERING_INDICES_PATH_SUFFIX, dtype=np.float)
    def get_clustering_indices(self):
        return self.compute_clustering_indices()

    def compute_degree_of_clustering(self) -> int:
        h_star = self.get_clustering_indices()
        d_of_c = np.array(h_star[h_star > 1] - 1).sum()
        if int(d_of_c) == 0:
            return 0.0001 # TODO this is a hack so that a downstream log does not fail

        return d_of_c

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripkey-K computation for disks of radius r < MAX_CELL_RADIUS
        :return: clustering indices for all r
        Was : clustering_index_point_process
        """
        logger.info("Running {} simulations of Ripley-K for {}",
                    constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], self._path)

        pixels_in_slice = numexpr.evaluate(constants.dataset_config["PIXELS_IN_SLICE"]).item()
        IF = self.get_intensities()
        cell_mask = self.get_cell_mask()
        IF = IF.astype(float) * cell_mask
        # TODO in VO we do not multiply by pixels_in_slice ???
        nuw = (np.sum(cell_mask[:, :] == 1))  # * pixels_in_slice  # whole surface of the cell
        my_lambda = float(np.sum(IF)) / float(nuw)  # volumic density
        k = self.ripley_k_random_measure_2D(IF, my_lambda, nuw)
        k_sim = np.zeros(
            (constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], constants.analysis_config["MAX_CELL_RADIUS"]))
        # simulate RIPLEY_K_SIMULATION_NUMBER list of random intensities and run ripley_k
        indsAll = np.where(cell_mask[:, :] == 1)
        for t in tqdm.tqdm(range(constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]), desc="Simulations"):
            inds_permuted = np.random.permutation(range(len(indsAll[0])))
            I_samp = np.zeros(IF.shape)
            for u in range(len(inds_permuted)):
                I_samp[indsAll[0][inds_permuted[u]], indsAll[1][inds_permuted[u]]] = IF[indsAll[0][u], indsAll[1][u]]
            k_sim[t, :] = self.ripley_k_random_measure_2D(I_samp, my_lambda, nuw).flatten()

        h = np.subtract(np.sqrt(k / math.pi), np.arange(1, constants.analysis_config["MAX_CELL_RADIUS"] + 1).reshape(
            (constants.analysis_config["MAX_CELL_RADIUS"], 1))).flatten()
        synth5, synth50, synth95 = helpers.compute_statistics_random_h_star_2d(k_sim)
        return helpers.compute_h_star_2d(h, synth5, synth50, synth95)

    def ripley_k_random_measure_2D(self, IF, my_lambda, nuw):
        IF_rev = IF[::-1, ::-1]
        P = signal.convolve(IF, IF_rev)
        dMap = np.zeros((P.shape[0], P.shape[1]))
        p, q = np.meshgrid(range(P.shape[0]), range(P.shape[0]))
        dMap = np.sqrt((p - IF.shape[0]) ** 2 + (q - IF.shape[1]) ** 2)
        # sum convolution using dMap
        K = np.zeros((constants.analysis_config["MAX_CELL_RADIUS"], 1))
        for dist in range(constants.analysis_config["MAX_CELL_RADIUS"]):
            K[dist] = P[dMap[:, :] <= dist].sum()
        K = K * (1 / (my_lambda * nuw)) - (1 / my_lambda)

        return K


class ImageWithIntensitiesAndMTOC(ImageWithMTOC, ImageWithIntensities):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithMTOC.is_a(repo, path) and ImageWithIntensities.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, cell_mask, quadrants_num=4) -> np.ndarray:
        """
        Given a quadrant mask and number of MTOC containing quadrant, compute 2D density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        IF = self.get_cytoplasmic_intensities()
        density_per_quadrant = np.zeros((quadrants_num, 2))
        # mark the MTOC quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            density_per_quadrant[quad_num, 0] = np.sum(IF[quadrant_mask == (quad_num + 1)]) / \
                                                (np.sum(
                                                    cell_mask[quadrant_mask == quad_num + 1]) * helpers.surface_coeff())

        if density_per_quadrant[:, 1].sum() != 1.0:
            raise (RuntimeError, "error in the MTOC quadrant detection for image %s" % self._path)

        return density_per_quadrant


class ImageWithSpotsAndMTOC(ImageWithMTOC, ImageWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithMTOC.is_a(repo, path) and ImageWithSpots.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, cell_mask, quadrants_num=4) -> np.ndarray:
        """
        Given an quadrant mask and the number of the MTOC containing quadrant, compute surfacic density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        surface_coeff = helpers.surface_coeff()
        spots = self.get_cytoplasmic_spots()
        density_per_quadrant = np.zeros((quadrants_num, 2))
        for spot in spots:
            spot_quad = quadrant_mask[spot[1], spot[0]]
            density_per_quadrant[spot_quad - 1, 0] += 1

        # mark the mtoc quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            density_per_quadrant[quad_num, 0] = density_per_quadrant[quad_num, 0] / (
                    np.sum(cell_mask[quadrant_mask == quad_num + 1]) * surface_coeff)

        return density_per_quadrant


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
        """compute nucleus surface in pixel using nucleus mask"""
        nucleus_mask = self.get_nucleus_mask()
        area = nucleus_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return area / len(self.get_multiple_nucleus_centroid())
        else:
            return area

    @helpers.checkpoint_decorator(NUCLEUS_AREA_PATH_SUFFIX, float)
    def get_nucleus_area(self) -> float:
        return self.compute_nucleus_area()

    def compute_cell_area(self):
        """compute cell surface in pixel using cell mask"""
        cell_mask = self.get_cell_mask()
        area = cell_mask.sum() * helpers.surface_coeff()  # * by pixel dimensions
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return area / len(self.get_multiple_nucleus_centroid())
        else:
            return area

    @helpers.checkpoint_decorator(CELL_AREA_PATH_SUFFIX, float)
    def get_cell_area(self) -> float:
        return self.compute_cell_area()


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
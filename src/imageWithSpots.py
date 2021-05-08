#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math
import numpy as np
import tqdm
from loguru import logger
import helpers
import constants
from repository import Repository
from sklearn.metrics.pairwise import pairwise_distances

from image import Image, ImageWithMTOC

from constants import SPOTS_PATH_SUFFIX
from constants import CYTOPLASMIC_SPOTS_PATH_SUFFIX
from constants import CLUSTERING_INDICES_PATH_SUFFIX
from constants import CYTOPLASMIC_SPOTS_PERIPHERAL_DISTANCE_PATH_SUFFIX

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

    # TODO : this will load the mask for each spot
    def compute_cytoplasmic_spots(self) -> np.ndarray:
        spots = self.get_spots()
        mask = [self.is_in_cytoplasm(s[0:2][::-1]) for s in spots]  # TODO check coordinate coherency for spots
        cytoplasmic_spots = spots[mask]
        return np.asarray(cytoplasmic_spots, dtype=int)

    @helpers.checkpoint_decorator(CYTOPLASMIC_SPOTS_PERIPHERAL_DISTANCE_PATH_SUFFIX, dtype=np.int)
    def get_cytoplasmic_spots_peripheral_distance(self):
        return self.compute_cytoplasmic_spots_peripheral_distance()

    def compute_cytoplasmic_spots_peripheral_distance(self) -> np.ndarray:
        """
        Return an array of distances to the periphery for cytoplasmic spots
        returns an array of int
        """
        spots = self.get_cytoplasmic_spots()
        distances = self.compute_cytoplasmic_coordinates_peripheral_distance(spots[:,0:2])
        return distances

    def compute_cytoplasmic_total_spots(self):
        cytoplasmic_spots = self.get_cytoplasmic_spots()
        return len(cytoplasmic_spots)

    def compute_median_cytoplasmic_distance_from_nucleus(self, dsAll) -> float:
        cytoplasm_mask = self.get_cytoplasm_mask()
        map_dist = np.multiply(cytoplasm_mask, dsAll)
        return np.median(map_dist[map_dist != 0])

    def compute_spots_normalized_distance_to_nucleus(self, quantile=.68) -> float:
        dists = constants.analysis_config['NUM_CONTOURS'] - self.get_cytoplasmic_spots_peripheral_distance()
        assert np.all(dists >= 0), "Negative distance to nucleus"
        normalized_dist_to_nucleus = np.quantile(dists, quantile) / constants.analysis_config['NUM_CONTOURS']

        return normalized_dist_to_nucleus

    def compute_spots_cytoplasmic_spread_entropy(self) -> float:
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
        return sd / np.median(d[d != 0])

    def compute_cytoplasmic_density(self):
        # compute mRNA density in the cytoplasm
        cytoplasmic_mrna_count = self.compute_cytoplasmic_total_spots()
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
        return: clustering indices for all r
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

    def compute_signal_from_periphery(self) -> np.ndarray:
        """
         np.ndarray of floats (total spots for each distance percentage from the periphery)
         normalization is the responsibility of the caller
         """
        spots_distances = self.get_cytoplasmic_spots_peripheral_distance()
        spots_counts = np.zeros(constants.analysis_config['NUM_CONTOURS'])
        for i in range(0, constants.analysis_config['NUM_CONTOURS']):
            spots_counts[i] = len(spots_distances[spots_distances <= i + 1])
        return spots_counts

    def compute_peripheral_total_spots(self):
        all_counts = self.compute_signal_from_periphery()
        return all_counts[constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']]


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


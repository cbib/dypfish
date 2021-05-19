#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math

import numexpr
import numpy as np
import tqdm
from loguru import logger

import constants
import helpers

from image3d import Image3d, Image3dWithMTOC
from imageWithSpots import ImageWithSpots
from repository import Repository

from constants import CLUSTERING_INDICES_PATH_SUFFIX


class Image3dWithSpots(Image3d, ImageWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3d.is_a(repo, path) and ImageWithSpots.is_a(repo, path)

    def compute_cytoplasmic_spots(self) -> np.ndarray:
        spots = self.get_spots()
        max_height = np.max(self.adjust_height_map()) # TODO this is a hack
        mask = [self.is_in_cytoplasm(s[0:2][::-1]) for s in spots]  # TODO check coordinate coherency for spots
        cytoplasmic_spots = spots[mask]
        cytoplasmic_spots = cytoplasmic_spots[cytoplasmic_spots[:,2] <= max_height]
        logger.info("Keeping {} cytoplasmic spots out of {} for {}",
                    len(cytoplasmic_spots), len(spots), self._path)
        return np.asarray(cytoplasmic_spots, dtype=int)

    def compute_cytoplasmic_spots_peripheral_distance(self)  -> np.ndarray:
        """
        Perform the computation in pseudo 3D using the height map
        Return an array of distances, values are np.nan if the spots is out of cytoplasm
        """
        spots = self.get_cytoplasmic_spots()
        logger.info("Computing 3D peripheral distance for {} spots in image {}", len(spots), self._path)
        distances = self.compute_cytoplasmic_coordinates_peripheral_distance(spots)
        return distances

    def ripley_k_point_process(self, nuw: float, my_lambda: float, spots=None, r_max: int = None) -> np.ndarray:
        if spots is None: spots = self.get_cytoplasmic_spots()
        n_spots = len(spots)
        r_max = r_max or constants.analysis_config["MAX_CELL_RADIUS"]
        pixels_in_slice = numexpr.evaluate(constants.dataset_config["PIXELS_IN_SLICE"]).item()

        K = np.zeros(r_max)
        for i in range(n_spots):
            # TODO : why only z_squared is computed in real size ?
            mask = np.zeros((n_spots, 3));
            mask[i, :] = 1
            other_spots = np.ma.masked_where(mask == 1, np.ma.array(spots, mask=False)).compressed().reshape(
                n_spots - 1, 3)
            x_squared = np.square(spots[i, 0] - other_spots[:, 0])
            y_squared = np.square(spots[i, 1] - other_spots[:, 1])
            z_squared = np.square(pixels_in_slice * (spots[i, 2] - other_spots[:, 2]))
            ds = np.sqrt(x_squared + y_squared + z_squared)
            if n_spots - 1 < r_max:
                for m in range(n_spots - 1):
                    K[math.ceil(ds[m]):r_max] = K[math.ceil(ds[m]):r_max] + 1
            else:
                for m in range(r_max):
                    K[m] = K[m] + ds[ds <= m].sum()
        K = K * (1 / (my_lambda ** 2 * nuw))
        return K

    def compute_random_cytoplasmic_spots_in_slices(self, num_spots, factor=100):
        slices = self.get_cell_mask_slices()
        height_map = self.get_cytoplasm_height_map()
        max_height = np.max(height_map)

        # generate random spot coordinates in a sphere and convert to int for compatibility
        # with masks; we generate 50 times more than the number of spots due to consecutive
        # selection criteria (within the cytoplasm and within slices)
        nucleus_centroid = self.get_nucleus_centroid()
        nucleus_centroid_z = height_map[nucleus_centroid[0], nucleus_centroid[1]] // 2
        center = np.array([nucleus_centroid[0], nucleus_centroid[1], nucleus_centroid_z])
        radius = self.compute_cell_diameter() // 2
        random_spots = helpers.random_points_in_sphere(center, radius, num_spots*factor).astype(int)
        random_spots_constrained_x = random_spots[(random_spots[:, 0] >= 1) &
                                                  (random_spots[:, 0] < height_map.shape[0])]
        random_spots_constrained_y = random_spots_constrained_x[(random_spots_constrained_x[:, 1] >= 1) &
                                                                (random_spots_constrained_x[:, 1] < height_map.shape[0])]
        random_spots_constrained_z = random_spots_constrained_y[(random_spots_constrained_y[:, 2] >= 0) &
                                                                (random_spots_constrained_y[:, 2] < max_height)]
        cytoplasm_mask = [self.is_in_cytoplasm(s[::-1]) for s in random_spots_constrained_z[:,0:2]]
        random_spots_in_cytoplasm = random_spots_constrained_z[cytoplasm_mask]
        slices_mask = slices[random_spots_in_cytoplasm[:, 0],
                             random_spots_in_cytoplasm[:, 1],
                             random_spots_in_cytoplasm[:, 2]]
        random_spots_in_slices = random_spots_in_cytoplasm[slices_mask == 1]

        if (random_spots_in_slices.shape[0] < num_spots): # 100 times has not been enough
            logger.warning("Was not able to generate {} random spots in {}", num_spots, self._path)
            return random_spots_in_slices

        return random_spots_in_slices[0:num_spots]

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripley-K computation in 3D for spheres of radius r < MAX_CELL_RADIUS
        return: clustering indices for all r
        """
        logger.info("Running {} simulations of Ripley-K for {} in 3D",
                    constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], self._path)
        spots = self.get_cytoplasmic_spots()
        n_spots = len(spots)
        nuw = self.compute_cell_volume()
        my_lambda = float(n_spots) / float(nuw)  # spot's volumic density

        k = self.ripley_k_point_process(nuw=nuw, my_lambda=my_lambda)  # TODO : first call for _all_ spots while the subsequent only for those in the height_map
        k_sim = np.zeros((constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], constants.analysis_config["MAX_CELL_RADIUS"]))

        # simulate RIPLEY_K_SIMULATION_NUMBER lists of random spots and run ripley_k
        for t in tqdm.tqdm(range(constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]), desc="Simulations"):
            random_spots = self.compute_random_cytoplasmic_spots_in_slices(len(spots), factor=200 )
            tmp_k = self.ripley_k_point_process(spots=random_spots, nuw=nuw, my_lambda=my_lambda).flatten()
            k_sim[t] = tmp_k

        h = np.subtract(np.power(((k * 3) / (4 * math.pi)), 1. / 3),
                        np.arange(1, constants.analysis_config["MAX_CELL_RADIUS"] + 1))
        synth5, synth50, synth95 = helpers.compute_statistics_random_h_star(h_sim=k_sim)
        return helpers.compute_h_star(h, synth5, synth50, synth95)

    @helpers.checkpoint_decorator(CLUSTERING_INDICES_PATH_SUFFIX, dtype=np.float)
    def get_clustering_indices(self):
        return self.compute_clustering_indices()

    def compute_degree_of_clustering(self) -> int:
        h_star = self.get_clustering_indices()
        d_of_c = np.array(h_star[h_star > 1] - 1).sum()
        if int(d_of_c) == 0:
            return 0.0001 # TODO this is a hack so that a downstream log does not fail

        return d_of_c

    def compute_peripheral_density(self):
        ''' compute mRNA density in the peripheral area '''
        peripheral_mrna_count = self.compute_peripheral_total_spots()
        if peripheral_mrna_count == 0:
            raise RuntimeError("Image contains no spots in periphery %s" % self._path)
        peripheral_volume = self.compute_peripheral_volume()
        return peripheral_mrna_count / peripheral_volume

    def compute_cytoplasmic_density(self):
        ''' compute mRNA density in the cytoplasm '''
        cytoplasmic_mrna_count = self.compute_cytoplasmic_total_spots()
        if cytoplasmic_mrna_count == 0:
            raise RuntimeError("Image contains no spots %s" % self._path)
        cytoplasmic_volume = self.compute_cytoplasmic_volume()
        return cytoplasmic_mrna_count / cytoplasmic_volume

    def compute_spots_cytoplasmic_spread_entropy(self) -> float:
        '''
        Computes entropy of the spatial distribution of cytoplasmic spots using
        Kozachenko-Leonenko entropy estimate
        '''
        cytoplasmic_spots = self.get_cytoplasmic_spots()
        # the number of neighbors k is chosen arbitrarily
        entropy = helpers.compute_entropy(cytoplasmic_spots, k=15, norm='euclidean')
        return entropy


class Image3dWithSpotsAndMTOC(Image3dWithMTOC, Image3dWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithMTOC.is_a(repo, path) and Image3dWithSpots.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4) -> np.ndarray:
        """
        compute volumic density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        volume_coeff = helpers.volume_coeff()
        height_map = self.adjust_height_map(cytoplasm=True)
        spots = self.get_cytoplasmic_spots()
        density_per_quadrant = np.zeros((quadrants_num, 2))
        for spot in spots:
            spot_quad = quadrant_mask[spot[1], spot[0]]
            if spot_quad == 0: continue
            density_per_quadrant[spot_quad - 1, 0] += 1

        # mark the mtoc quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            quadrant_volume = np.sum(height_map[quadrant_mask == quad_num + 1]) * volume_coeff
            density_per_quadrant[quad_num, 0] = density_per_quadrant[quad_num, 0] / quadrant_volume

        if density_per_quadrant[:, 1].sum() != 1.0:
            raise (RuntimeError, "error in the MTOC quadrant detection for image %s" % self._path)

        return density_per_quadrant


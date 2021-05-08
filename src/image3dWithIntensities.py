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
from imageWithIntensities import ImageWithIntensities
from repository import Repository


class Image3dWithIntensities(Image3d, ImageWithIntensities):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3d.is_a(repo, path) and ImageWithIntensities.is_a(repo, path)

    def compute_cytoplasmic_density(self):
        # compute signal density of the cytoplasm
        cytoplasmic_intensity_count = self.compute_cytoplasmic_total_intensity()
        cytoplasmic_volume = self.compute_cytoplasmic_volume()
        return cytoplasmic_intensity_count / cytoplasmic_volume

    def compute_intensities_normalized_distance_to_nucleus(self, quartile=.68) -> float:
        height_map = self.get_cytoplasm_height_map()
        IF = self.compute_cytoplasmic_intensities()
        mean_signal = np.mean(IF[IF > 0])
        peaks = np.argwhere(IF > mean_signal * 1.8) # arbitrary value
        peaks_z = np.array([height_map[p[0], p[1]] // 2 for p in peaks])
        peaks = np.hstack((peaks, peaks_z[:, None])).astype(int)
        dists = self.compute_cytoplasmic_coordinates_peripheral_distance(peaks)
        valid_dists = dists[~np.isnan(dists)]
        dists = constants.analysis_config['NUM_CONTOURS'] - valid_dists
        assert np.all(dists >= 0), "Negative distance to nucleus"
        normalized_dist_to_nucleus = np.quantile(dists, quartile) / constants.analysis_config['NUM_CONTOURS']

        return normalized_dist_to_nucleus

    def compute_intensities_normalized_cytoplasmic_spread(self):
        IF = self.compute_cytoplasmic_intensities()
        height_map = self.adjust_height_map(cytoplasm=True)
        IF = np.multiply(IF, height_map) # factoring in the 3D
        cytoplasm_mask = self.get_cytoplasm_mask()

        # Calculate the spread of signal peaks
        mean_signal = np.mean(IF[cytoplasm_mask == 1])
        peaks = np.argwhere(IF > mean_signal * 2)  # arbitrary choice to reduce the number of peaks

        mu_x = peaks[:, 0].sum() / len(peaks)
        mu_y = peaks[:, 1].sum() / len(peaks)
        sd = math.sqrt(np.sum((peaks[:, 0] - mu_x) ** 2) / len(peaks) +
                       np.sum((peaks[:, 1] - mu_y) ** 2) / len(peaks))

        diameter = self.compute_cell_diameter()
        return sd / (0.68 * diameter / 2)

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripkey-K computation for disks of radius r < MAX_CELL_RADIUS
        return: clustering indices for all r
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



class Image3dWithIntensitiesAndMTOC(Image3dWithMTOC, Image3dWithIntensities):

    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithMTOC.is_a(repo, path) and Image3dWithIntensities.is_a(repo, path)

    def compute_cytoplasmic_density(self):
        # compute signal density of the cytoplasm
        cytoplasmic_intensity_count = self.compute_cytoplasmic_total_intensity()
        cytoplasmic_volume = self.compute_cytoplasmic_volume()
        return cytoplasmic_intensity_count / cytoplasmic_volume

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4) -> np.ndarray:
        """
        compute volumic density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        IF = self.compute_cytoplasmic_intensities()
        height_map = self.adjust_height_map(cytoplasm=True)
        density_per_quadrant = np.zeros((quadrants_num, 2))
        # mark the MTOC quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            height = np.sum(height_map[quadrant_mask == quad_num + 1])
            quadrant_volume = height * helpers.volume_coeff()
            density_per_quadrant[quad_num, 0] = np.sum(IF[quadrant_mask == (quad_num + 1)]) / quadrant_volume

        if density_per_quadrant[:, 1].sum() != 1.0:
            raise (RuntimeError, "error in the MTOC quadrant detection for image %s" % self._path)

        return density_per_quadrant


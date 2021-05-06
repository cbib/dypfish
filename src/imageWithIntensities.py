#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math

import numexpr
import numpy as np
import tqdm
from loguru import logger
from scipy import signal
from sklearn.metrics.pairwise import pairwise_distances

import constants
import helpers
import image_processing as ip
from constants import CLUSTERING_INDICES_PATH_SUFFIX
from constants import IF_PATH_SUFFIX
from image import Image, ImageWithMTOC
from repository import Repository


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
        returns: np.ndarray of floats (intensities of the cell)
        """
        descriptor = self._path + IF_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No intensities for image %s" % self._path)
        raw_value = self._repository.get(descriptor)
        intensities = np.array(raw_value)
        cell_mask = self.get_cell_mask()
        return np.multiply(intensities, cell_mask)

    def compute_cell_total_intensity(self) -> float:
        intensities = self.get_intensities()
        cell_mask = self.get_cell_mask()
        cell_intensities = np.multiply(intensities, cell_mask)
        return cell_intensities.sum()

    def compute_cytoplasmic_intensities(self) -> np.ndarray:
        intensities = self.get_intensities()
        cytoplasm_mask = self.get_cytoplasm_mask()
        cytoplasmic_intensities = np.multiply(intensities, cytoplasm_mask)
        return cytoplasmic_intensities

    def compute_cytoplasmic_total_intensity(self) -> float:
        cytoplasmic_intensities = self.compute_cytoplasmic_intensities()
        return cytoplasmic_intensities.sum()

    def compute_peripheral_total_intensity(self) -> float:
        intensities = self.compute_cytoplasmic_intensities()
        peripheral_mask = self.compute_peripheral_mask()
        peripheral_intensity_mask = np.multiply(intensities, peripheral_mask)
        return peripheral_intensity_mask.sum()

    def compute_signal_from_periphery(self) -> np.ndarray:
        """
         np.ndarray of floats (total intensity for each distance percentage from the periphery)
         normalization is the responsibility of the caller
         """
        intensities = np.multiply(self.get_intensities(), self.get_cytoplasm_mask())
        cell_mask_distance_map = self.get_cell_mask_distance_map()
        intensities_sums = np.zeros(constants.analysis_config['NUM_CONTOURS'])
        for i in range(0, constants.analysis_config['NUM_CONTOURS']):
            intensities_sums[i] = intensities[cell_mask_distance_map <= i + 1].sum()
        return intensities_sums

    def compute_median_cytoplasmic_distance_from_nucleus2d(self, dsAll):
        cytoplasm_mask = self.get_cytoplasm_mask()
        distances = np.multiply(dsAll, cytoplasm_mask)
        return np.median(distances[distances != 0]), np.max(distances[distances != 0])

    def compute_intensities_normalized_spread_to_centroid(self) -> float:
        nucleus_centroid = self.get_nucleus_centroid()
        IF = self.compute_cytoplasmic_intensities()
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid)  # 2d distances from nucleus_centroid
        dsAll = dsAll * self.get_cytoplasm_mask()
        median_dist, max_dist = self.compute_median_cytoplasmic_distance_from_nucleus2d(dsAll)

        # Calculate the distances of signal peaks to nucleus_centroid
        mean_signal = np.mean(IF[IF > 0])
        peaks = np.argwhere(IF > mean_signal * 1.5)  # arbitrary choice to reduce the number of peaks
        dsPeaks = np.sqrt(np.sum((peaks - [nucleus_centroid[0], nucleus_centroid[1]]) ** 2, axis=1))

        spread_to_centroid = np.median(dsPeaks) / median_dist
        return spread_to_centroid

    def compute_intensities_normalized_cytoplasmic_spread(self):
        IF = self.compute_cytoplasmic_intensities()
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
        cytoplasmic_intensity_count = self.compute_cytoplasmic_total_intensity()
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
        Point process Ripkey-K computation in 2D for disks of radius r < MAX_CELL_RADIUS
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

    def ripley_k_random_measure_2D(self, IF, my_lambda, nuw):
        IF_rev = IF[::-1, ::-1]
        P = signal.convolve(IF, IF_rev)
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
        IF = self.compute_cytoplasmic_intensities()
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


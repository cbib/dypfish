#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

from loguru import logger

import constants
import path
from image3dWithIntensities import Image3dWithIntensities
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithIntensities(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithIntensities(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 2052091.52470807, places=5)

    def test_compute_intensities_normalized_distance_to_nucleus(self):
        normalized_value = self.img.compute_intensities_normalized_distance_to_nucleus()
        self.assertAlmostEqual(normalized_value, 0.67, places=3)

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 1771979.4115436368, places=5)

    def test_compute_intensities_normalized_cytoplasmic_spread(self):
        spread = self.img.compute_intensities_normalized_cytoplasmic_spread()
        self.assertAlmostEqual(spread, 0.611584381061, places=5)

    def test_compute_cytoplasmic_intensities(self):
        result = self.img.compute_cytoplasmic_intensities()
        self.assertEqual(result.sum(), 1178283127.0)

    def test_compute_signal_from_periphery(self):
        peripheral_intensities = self.img.compute_signal_from_periphery()
        self.assertEqual(peripheral_intensities.shape[0], 100)
        # test an arbitrary value
        self.assertEqual(peripheral_intensities[10], 22178801.0)
        self.assertTrue(all(
            peripheral_intensities[i] <= peripheral_intensities[i + 1] for i in range(len(peripheral_intensities) - 1)))
        self.assertAlmostEqual(peripheral_intensities.sum(), 50951395310.0)

    def test_compute_clustering_indices(self):
        logger.error("This function is not tested")
        self.fail()

    def test_ripley_k_random_measure_2d(self):
        IF = self.img.get_intensities()
        test_k = np.array([[2.53917575e+04, 1.26888405e+05, 3.29580740e+05, 7.33753607e+05, 1.23718579e+06]])
        nuw = (np.sum(self.img.get_cell_mask()[:, :] == 1))  # * pixels_in_slice  # whole surface of the cell
        my_lambda = float(np.sum(IF)) / float(nuw)  # volumic density
        K = self.img.ripley_k_random_measure_2D(IF, my_lambda, nuw)
        self.assertAlmostEqual(np.sum(K[0:5]), np.sum(test_k), places=2)

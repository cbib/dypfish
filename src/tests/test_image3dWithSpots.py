#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

import numpy as np
from loguru import logger

import constants
import path
from image3dWithSpots import Image3dWithSpots
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)

class TestImage3dWithSpots(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def test_compute_cytoplasmic_spots(self):
        # image 10 has very high z coordinates for spots
        self.img2 = Image3dWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/10")
        result = self.img2.compute_cytoplasmic_spots()
        self.assertEqual(len(result), 33)

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cytoplasmic_spots_peripheral_distance(self):
        results = self.img.compute_cytoplasmic_spots_peripheral_distance()
        self.assertEqual(len(results), 153)
        # arbitrarily check two values
        self.assertEqual(results[3], 29)
        self.assertEqual(results[65], 28)
        self.assertEqual(results.sum(), 7538.0)

    def test_compute_spots_normalized_distance_to_nucleus(self):
        relative_to_centroid = self.img.compute_spots_normalized_distance_to_nucleus()
        self.assertAlmostEqual(relative_to_centroid, 0.58, places=5)

    def test_compute_spots_cytoplasmic_spread_entropy(self):
        spread = self.img.compute_spots_cytoplasmic_spread_entropy()
        self.assertAlmostEqual(spread, 16.334849331, places=5)

    def test_clustering_index_point_process(self):
        np.random.seed(0)
        h_star = self.img.compute_clustering_indices()
        self.assertEqual(len(h_star), 300)
        self.assertAlmostEqual(h_star.sum(), 35.014805372158634)  # might not work since random

    def test_ripley_k_point_process(self):
        K = self.img.ripley_k_point_process(nuw=1158349.7249999999, my_lambda=0.00018819877563315348)
        self.assertAlmostEqual(K.sum(), 110417705.92985436, places=3)
        self.assertEqual(K.shape, (300,))

    def test_compute_degree_of_clustering(self):
        np.random.seed(0)
        self.assertAlmostEqual(self.img.compute_degree_of_clustering(), 23.95041457661648, places=5)

    def test_compute_mrna_density(self):
        mrna_density = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(mrna_density, 0.15130537974, places=5)

    def test_compute_random_cytoplasmic_spots_in_slices(self):
        np.random.seed(0)
        random_spots = self.img.compute_random_cytoplasmic_spots_in_slices(100)
        self.assertTrue(np.all([self.img.is_in_cytoplasm(s[::-1]) for s in random_spots[:,0:2]]))
        self.assertEqual(random_spots.shape[0], 100)

    def test_compute_signal_from_periphery(self):
        # calls the code in the super class (ImageWithSpots)
        # image 10 has spots with very hight z coordinates
        self.img2 = Image3dWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/10")
        spot_counts = self.img2.compute_signal_from_periphery()
        self.assertEqual(spot_counts[99], self.img2.compute_cytoplasmic_total_spots())



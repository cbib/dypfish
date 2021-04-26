#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

import numpy as np

import constants
import path
from image3d import Image3dWithSpots
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithSpots(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_spots_peripheral_distance_3d(self):
        results = self.img.compute_spots_peripheral_distance_3d()
        self.assertGreater(len(results), 100)  # more than 100 spots
        # arbitrarily check two values
        self.assertEqual(results[0], 75)
        self.assertEqual(results[65], 16)
        self.assertEqual(results.sum(), 7313) # TODO: double check this

    def test_compute_spots_cytoplasmic_spread(self):
        normalized_average_2d_distance = self.img.compute_spots_cytoplasmic_spread()
        self.assertAlmostEqual(normalized_average_2d_distance, 1.1262934840267351)

    def test_clustering_index_point_process(self):
        np.random.seed(0)
        h_star = self.img.compute_clustering_indices()
        self.assertEqual(len(h_star), 300)
        self.assertAlmostEqual(h_star.sum(), 289.98034888178876)  # I think this is ok, but might not work since random

    def test_ripley_k_point_process(self):
        K = self.img.ripley_k_point_process(nuw=1158349.7249999999, my_lambda=0.00018819877563315348)
        self.assertAlmostEqual(K.sum(), 234044508.31346154)
        self.assertEqual(K.shape, (300,))

    def test_compute_degree_of_clustering(self):
        np.random.seed(0)
        self.assertAlmostEqual(self.img.compute_degree_of_clustering(), 71.93654959822206)

    def test_spots_in_slices(self):
        spots = self.img.get_spots()
        spots_in_slices = self.img.compute_spots_in_slices()
        self.assertLess(len(spots_in_slices), len(spots))

    def test_compute_mrna_density(self):
        mrna_density = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(mrna_density, 0.128494541373)

    def test_compute_spots_cytoplasmic_spread(self):
        result = self.img.compute_spots_cytoplasmic_spread()
        self.assertEqual(result, 0.9253082533496478)
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import warnings
from unittest import TestCase

import constants
import image_processing as ip
import path
from imageWithSpots import ImageWithSpots
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImageWithSpots(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = ImageWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_get_spots(self):
        self.assertEqual(self.img.get_spots().shape, (218, 3))

    def test_get_cytoplasmic_spots(self):
        cytoplasmic_spots = self.img.get_cytoplasmic_spots()
        self.assertEqual(cytoplasmic_spots.shape, (155, 3))
        self.assertTrue((cytoplasmic_spots[0] == [302, 123, 12]).all())

    def test_compute_peripheral_total_spots(self):
        spots_num = self.img.compute_peripheral_total_spots()
        self.assertEqual(spots_num, 18.0)

    def test_compute_spots_peripheral_distance_2d(self):
        peripheral_distance_2D = self.img.compute_cytoplasmic_spots_peripheral_distance()
        self.assertEqual(peripheral_distance_2D[0], 46)
        self.assertEqual(peripheral_distance_2D.sum(), 9247)
        self.assertEqual(peripheral_distance_2D.size, len(self.img.get_cytoplasmic_spots()))

    def test_compute_cytoplasmic_spots(self):
        self.assertEqual(len(self.img.compute_cytoplasmic_spots()), 155)

    def test_compute_cytoplasmic_total_spots(self):
        self.assertEqual(self.img.compute_cytoplasmic_total_spots(), 155)

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 0.3012427677, places=5)

    def test_compute_median_cytoplasmic_distance_from_nucleus(self):
        warnings.warn("This function is not sufficiently tested", RuntimeWarning)
        nucleus_centroid = self.img.get_nucleus_centroid()
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid)
        result = self.img.compute_median_cytoplasmic_distance_from_nucleus(dsAll)
        self.assertAlmostEqual(result, 109.17875251164, places=3)

    def test_compute_spots_normalizaed_distance_to_nucleus(self):
        normalized_average_2d_distance = self.img.compute_spots_normalized_distance_to_nucleus()
        self.assertAlmostEqual(normalized_average_2d_distance, 0.49720000000000014, places=5)

    def test_compute_spots_normalized_cytoplasmic_spread(self):
        result = self.img.compute_spots_cytoplasmic_spread_entropy()
        self.assertAlmostEqual(result, 0.783341863054, places=5)

    def test_compute_random_spots(self):
        random_spots = self.img.compute_random_spots()
        self.assertTrue(random_spots.shape == (218, 2))

    def test_compute_signal_from_periphery(self):
        peripheral_spots = self.img.compute_signal_from_periphery()
        self.assertEqual(peripheral_spots.shape[0], 100)
        # test an arbitrary value
        self.assertEqual(peripheral_spots[30], 18.0)
        self.assertTrue(all(
            peripheral_spots[i] <= peripheral_spots[i + 1] for i in range(len(peripheral_spots) - 1)))
        self.assertEqual(peripheral_spots.sum(), 6408.0)
        self.assertEqual(peripheral_spots[99], self.img.compute_cytoplasmic_total_spots())

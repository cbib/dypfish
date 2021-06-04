#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import unittest
from unittest import TestCase

import constants
import image_processing as ip
import path
from imageWithIntensities import ImageWithIntensities
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)  # TODO this is annoying


class TestImageWithIntensities(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = ImageWithIntensities(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_get_intensities(self):
        self.assertEqual(self.img.get_intensities().shape, (512, 512))
        total_intensities = self.img.compute_cell_total_intensity()
        self.assertEqual(total_intensities, 1746049018.0)

    def test_compute_cytoplasmic_intensities(self):
        self.assertEqual(self.img.compute_cytoplasmic_intensities().shape, (512, 512))
        self.assertEqual(self.img.compute_cytoplasmic_total_intensity(), 1373183555.0)

    def test_compute_signal_from_periphery(self):
        peripheral_intensities = self.img.compute_signal_from_periphery()
        self.assertEqual(peripheral_intensities.shape[0], 100)
        # test an arbitrary value
        self.assertEqual(peripheral_intensities[10], 132752104.0)
        self.assertTrue(all(
            peripheral_intensities[i] <= peripheral_intensities[i + 1] for i in range(len(peripheral_intensities) - 1)))
        self.assertAlmostEqual(peripheral_intensities.sum(), 68177643289.0)

    def test_compute_cell_total_intensity(self):
        cell_intensity = self.img.compute_cell_total_intensity()
        self.assertEqual(cell_intensity, 1746049018.0)

    def test_compute_median_cytoplasmic_distance_from_nucleus2d(self):
        nucleus_centroid = self.img.get_nucleus_centroid()
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid)
        median, max = self.img.compute_median_cytoplasmic_distance_from_nucleus2d(dsAll)
        self.assertAlmostEqual(median, 97.0824391947, places=3)
        self.assertAlmostEqual(max, 149.9333185119, places=3)

    def test_compute_intensities_normalized_distance_to_nucleus(self):
        normalized_value = self.img.compute_intensities_normalized_distance_to_nucleus()
        self.assertAlmostEqual(normalized_value, 0.58, places=5)

    def test_compute_intensities_cytoplasmic_spread(self):
        spread = self.img.compute_intensities_normalized_cytoplasmic_spread()
        self.assertAlmostEqual(spread, 0.78921627304, places=5)

    def test_compute_peripheral_total_intensity(self):
        intensity = self.img.compute_peripheral_total_intensity()
        self.assertAlmostEqual(intensity, 382727324.0)
        self.assertLess(intensity, self.img.compute_cytoplasmic_total_intensity())

    def test_compute_cytoplasmic_total_intensity(self):
        intensity = self.img.compute_cytoplasmic_total_intensity()
        self.assertAlmostEqual(intensity, 1373183555.0)

if __name__ == '__main__':
    unittest.main()

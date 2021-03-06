#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import numpy as np
import path
import constants

from repository import H5RepositoryWithCheckpoint
from image import ImageWithSpots

constants.init_config(analysis_config_js_path=path.test_config_path)  # TODO this is annoying


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

    def test_compute_cell_mask_distance_map(self):
        distance_mask = self.img.compute_cell_mask_distance_map()
        self.assertEqual(distance_mask.shape, (512, 512))
        self.assertEqual(np.min(distance_mask), 0)
        self.assertEqual(np.max(distance_mask), 100)
        self.assertEqual(np.sum(distance_mask), 2024717)

    def test_compute_spots_peripheral_distance_2d(self):
        peripheral_distance_2D = self.img.compute_spots_peripheral_distance_2D()
        self.assertEqual(peripheral_distance_2D[0], 46)
        self.assertEqual(peripheral_distance_2D.size, 155)

    def test_compute_cytoplasmic_spots(self):
        self.assertEqual(len(self.img.compute_cytoplasmic_spots()), 155)

    def test_compute_cytoplasmic_total_spots(self):
        self.assertEqual(self.img.compute_cytoplasmic_total_spots(), 155)

    def test_compute_spots_cytoplasmic_spread(self):
        normalized_average_2d_distance = self.img.compute_spots_cytoplasmic_spread()
        self.assertAlmostEqual(normalized_average_2d_distance, 0.9366763988321943)

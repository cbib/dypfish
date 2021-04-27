#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import warnings
from unittest import TestCase

import numpy as np

import constants
import path
from image3d import Image3dWithSpotsAndMTOC
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithSpotsAndMTOC(TestCase):

    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithSpotsAndMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertEqual(0.1284945413733293, result)

    def test_compute_density_per_quadrant(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_density_per_quadrant(mtoc_quad, quadrant_mask)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 0.13705567289739, places=3)
        self.assertAlmostEqual(result[:, 0].sum(), 0.5481894332973, places=3)
        self.assertEqual(result[:, 1].sum(), 1.0)

    def test_compute_quadrant_densities(self):
        result = self.img.compute_quadrant_densities()
        self.assertAlmostEqual(result[:, 0].sum(), 0.529914777)
        self.assertEqual(result[3, 1], 1)

    # def test_split_in_quadrants(self):
    #     warnings.warn(
    #         "test_compute_max_density_MTOC_quadrant test not well tested",
    #         RuntimeWarning
    #     )
    #     test_array = np.array([(0.12947032, 0.), (0.15919367, 0.), (0.06432269, 0.), (0.17692811, 1.)])
    #     results = self.img.split_in_quadrants()
    #     self.assertEqual(np.shape(results), np.shape(test_array))
    #     self.assertAlmostEqual(np.sum(results[:,0].sum()), np.sum(test_array[:,0].sum()), places=7)

    def test_compute_density_per_quadrant_and_slices(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        result = self.img.compute_density_per_quadrant_and_slices(quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result.sum(), 54.2215439855)
        self.assertAlmostEqual(result[2], 3.07907363333)

    def test_compute_peripheral_density_per_quadrant_and_slices(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        result = self.img.compute_peripheral_density_per_quadrant_and_slices(quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result.sum(), 6.08385417)
        self.assertAlmostEqual(result[7], 1.1286306931)

    def test_split_in_quadrants_and_slices(self, quadrants_num=4, stripes=3):
        result = self.img.split_in_quadrants_and_slices()
        self.assertAlmostEqual(result.sum(), 61.0562912512)
        self.assertAlmostEqual(result[2], 2.952970008)

    def test_compute_peripheral_density_per_quadrant(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        height_map = self.img.adjust_height_map()
        result = self.img.compute_peripheral_density_per_quadrant(mtoc_quad, quadrant_mask, height_map)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 0.039297451478, places=3)
        self.assertAlmostEqual(result[:, 0].sum(), 0.3030269392, places=3)
        self.assertEqual(result[:, 1].sum(), 1.0)
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

    def test_compute_peripheral_density_per_quadrant(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_peripheral_density_per_quadrant(mtoc_quad, quadrant_mask)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 0.0392974514788, places=3)
        self.assertAlmostEqual(result[:, 0].sum(), 0.3030269392, places=3)
        self.assertEqual(result[:, 1].sum(), 1.0)

    def test_compute_quadrant_densities(self):
        result = self.img.compute_quadrant_densities()
        self.assertAlmostEqual(result[:, 0].sum(), 0.529914777)
        self.assertEqual(result[3, 1], 1)

    def test_compute_density_per_quadrant_and_slices(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_density_per_quadrant_and_slices(mtoc_quad, quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result[:,0].sum(), 10.7932895271)
        self.assertAlmostEqual(result[:,1].sum(), 3)
        self.assertAlmostEqual(result[3,0], 0.4972141871094)

    def test_compute_peripheral_density_per_quadrant_and_slices(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_peripheral_density_per_quadrant_and_slices(mtoc_quad, quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result[:,0].sum(), 7.074859430)
        self.assertAlmostEqual(result[7,0], 0.305829734546)

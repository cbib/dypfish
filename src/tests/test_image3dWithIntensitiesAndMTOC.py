#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski


import pathlib
from unittest import TestCase

import constants
import path
from image3dWithIntensities import Image3dWithIntensitiesAndMTOC
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithIntensitiesAndMTOC(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithIntensitiesAndMTOC(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_density_per_quadrant(self):
        quadrant_mask = self.img.rotate_quadrant_mask(45, 4).astype(int)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_density_per_quadrant(mtoc_quad, quadrant_mask)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 2162328.90115262, places=3)
        self.assertEqual(result[:, 1].sum(), 1.0)

    def test_compute_peripheral_density_per_quadrant(self):
        # implemented in the parent class, calls compute_density_per_quadrant
        quadrant_mask = self.img.rotate_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_peripheral_density_per_quadrant(mtoc_quad, quadrant_mask)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 3626622.27561184, places=3)
        self.assertAlmostEqual(result[:, 0].sum(), 22299061.50705328, places=3)
        self.assertEqual(result[:, 1].sum(), 1.0)

    def test_compute_quadrant_densities(self):
        # implemented in the parent class, calls compute_density_per_quadrant
        result = self.img.compute_quadrant_densities()
        self.assertAlmostEqual(result[:, 0].sum(), 9092726.97793292, places=3)
        self.assertEqual(result[2, 1], 1)

    def test_compute_density_per_quadrant_and_slices(self):
        # implemented in the parent class, calls compute_density_per_quadrant
        quadrant_mask = self.img.rotate_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_density_per_quadrant_and_slices(mtoc_quad, quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result[:,0].sum(), 36121655.010594, places=3)
        self.assertAlmostEqual(result[3,0], 9911765.260131167)
        self.assertEqual(result[:, 1].sum(), 3)

    def test_compute_peripheral_density_per_quadrant_and_slices(self):
        # implemented in the parent class, calls compute_density_per_quadrant
        quadrant_mask = self.img.rotate_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        result = self.img.compute_peripheral_density_per_quadrant_and_slices(mtoc_quad, quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result.sum(), 79737408.88673353, places=3)
        self.assertAlmostEqual(result[7,0], 10437447.76977793, places=3)
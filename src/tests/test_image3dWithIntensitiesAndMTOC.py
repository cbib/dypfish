#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski


import pathlib
from unittest import TestCase
import numpy as np
import warnings
import constants
import path
from repository import H5RepositoryWithCheckpoint
from image3d import Image3dWithIntensitiesAndMTOC

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithIntensitiesAndMTOC(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithIntensitiesAndMTOC(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 1759171.4466687627)

    def test_split_in_quadrants(self):
        test_array = np.array([(2928540.93937, 0), (2614853.08760, 0), (1675368.87807, 1), (1619814.17746, 0)])
        result = self.img.split_in_quadrants()
        self.assertEqual(np.shape(result), np.shape(test_array))
        self.assertAlmostEqual(np.sum(result[:,0]), np.sum(test_array[:,0]), places=2)

    def test_split_in_quadrants_and_slices(self, quadrants_num=4, stripes=3):
        result = self.img.split_in_quadrants_and_slices()
        self.assertAlmostEqual(result.sum(), 25.3735566411)
        self.assertAlmostEqual(result[2], 1.90160983)

    def test_compute_density_per_quadrant(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        mtoc_position = self.img.get_mtoc_position()
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        height_map = self.img.adjust_height_map(cytoplasm=True)
        result = self.img.compute_density_per_quadrant(mtoc_quad, quadrant_mask, height_map)

        self.assertEqual(result.shape, (4, 2))
        self.assertAlmostEqual(result[1, 0], 2114080.4358608127)
        self.assertEqual(result[:, 1].sum(), 1.0)

    def test_compute_peripheral_density_per_quadrant_and_slices(self):
        quadrant_mask = self.img.compute_quadrant_mask(45, 4)
        result = self.img.compute_peripheral_density_per_quadrant_and_slices(quadrant_mask, stripes=3, quadrants_num=4)
        self.assertAlmostEqual(result.sum(), 14.56655481)
        self.assertAlmostEqual(result[7], 1.9488238603)
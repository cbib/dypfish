#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import constants
import path
from repository import H5RepositoryWithCheckpoint
from image3d import Image3dWithIntensities

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithIntensities(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithIntensities(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_intensities_normalized_spread_to_centroid(self):
        normalized_value = self.img.compute_intensities_normalized_spread_to_centroid()
        self.assertAlmostEqual(normalized_value, 1.073387706817907, places=3)

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertEqual(1759171.4466687627, result)

    def test_compute_intensities_normalized_cytoplasmic_spread(self):
        spread = self.img.compute_intensities_normalized_cytoplasmic_spread()
        self.assertAlmostEqual(spread, 0.787203983556, places=5)

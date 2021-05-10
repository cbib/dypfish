#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

from loguru import logger

import constants
import path
from image3dWithIntensities import Image3dWithIntensities
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImage3dWithIntensities(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithIntensities(repository=self.repo, image_path="protein/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 2052091.52470807, places=5)

    def test_compute_intensities_normalized_distance_to_nucleus(self):
        normalized_value = self.img.compute_intensities_normalized_distance_to_nucleus()
        self.assertAlmostEqual(normalized_value, 0.7, places=3)

    def test_compute_cytoplasmic_density(self):
        result = self.img.compute_cytoplasmic_density()
        self.assertAlmostEqual(result, 2065083.4523325, places=5)

    def test_compute_intensities_normalized_cytoplasmic_spread(self):
        spread = self.img.compute_intensities_normalized_cytoplasmic_spread()
        self.assertAlmostEqual(spread, 0.609855834172, places=5)

    def test_compute_clustering_indices(self):
        logger.error("This function is not tested")
        self.fail()

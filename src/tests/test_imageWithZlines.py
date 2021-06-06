#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

from loguru import logger
import numpy as np

import constants
import path
from imageWithZlines import imageMultiNucleusWithSpotsAndZlines, ImageMultiNucleus, imageWithSpotsAndZlines, ImageMultiNucleusWithSpots

from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImageMultiNucleus(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = ImageMultiNucleus(repository=self.repo, image_path="mrna/actn2/immature/02")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_get_multiple_nucleus_centroid(self):
        centroids = self.img.get_multiple_nucleus_centroid()
        valid_centroids = [[873, 97], [108, 109], [1608, 121]]
        self.assertTrue(np.array_equal(centroids, valid_centroids))
        for i, centroid in enumerate(centroids):
            self.assertTrue(np.array_equal(centroid, valid_centroids[i]))

    def test_compute_nucleus_area(self):
        nucleus_area = self.img.compute_nucleus_area()
        self.assertEqual(122.20381328073636, nucleus_area)

    def test_compute_cell_area(self):
        cell_area = self.img.compute_cell_area()
        self.assertEqual(530.3914091606399, cell_area)


class TestimageWithSpotsAndZlines(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = imageWithSpotsAndZlines(repository=self.repo, image_path="mrna/actn2/immature/02")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_get_z_lines_masks(self):
        z_lines_masks = self.img.get_z_lines_masks()
        whole_area = 0
        for z_lines_mask in z_lines_masks:
            whole_area += np.sum(z_lines_mask)
        self.assertEqual(1137023, whole_area)

    def test_compute_minimal_z_line_distance(self):
        minimal_z_line_distance = self.img.compute_minimal_z_line_distance(z_line_spacing=15)
        print(np.sum(minimal_z_line_distance))
        minimal_z_line_distance_test = [0.10612387, 0.01322199, 0.01043841, 0.01461378, 0.01078636, 0.00904662,
                                        0.01078636, 0.0118302, 0.00974252, 0.0059151, 0.00382742, 0.0045233,
                                        0.00243563, 0.00487126, 0.00417537
                                        ]
        self.assertAlmostEqual(np.sum(minimal_z_line_distance), np.sum(minimal_z_line_distance_test), places=6)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import numpy as np

from image3d import Image3dWithSpots
from repository import H5RepositoryWithCheckpoint
import path
import constants
import helpers
import image_processing as ip

constants.init_config(analysis_config_js_path=path.test_config_path)  # TODO this is annoying


class Test(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3dWithSpots(repository=self.repo, image_path="mrna/arhgdia/2h/1")
        self.nucleus_mask = self.img.get_nucleus_mask()
        cell_mask = self.img.get_cell_mask()
        self.cytoplasm_mask = (cell_mask == 1) & (self.nucleus_mask == 0)

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_line_segments(self):  # TODO optimize for performance
        nucleus_centroid = self.img.get_nucleus_centroid()
        x_slope = -7.0
        y_slope = 10.0
        segments = ip.compute_nucleus_and_cytoplasm_line_segments(self.nucleus_mask, self.cytoplasm_mask,
                                                                  nucleus_centroid, x_slope, y_slope)
        self.assertEqual(segments[0].sum(), 4)  # nucleus
        self.assertEqual(segments[1].sum(), 8)  # cytoplasm

    def test_compute_cell_mask_distance_map(self):
        nucleus_mask = self.img.get_nucleus_mask()
        nucleus_centroid = self.img.get_nucleus_centroid()
        cytoplasm_mask = self.img.get_cytoplasm_mask()
        contour_points = ip.compute_contour_points(nucleus_mask, nucleus_centroid, cytoplasm_mask)
        dm = ip.compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points)
        self.assertEqual(dm.sum(), 2024717)  # self.assertEqual(dm.sum(), 2137740) # new code version

        cell_mask = helpers.unit_circle(45, 12.5)
        nucleus_mask = helpers.unit_circle(45, 5.5)
        nucleus_centroid = [22, 22]
        cytoplasm_mask = cell_mask - nucleus_mask
        contour_points = ip.compute_contour_points(nucleus_mask, nucleus_centroid, cytoplasm_mask, num_contours=4,
                                                   max_cell_radius=13, image_width=45, image_height=45)
        dm = ip.compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points, num_contours=4)
        self.assertEqual(dm.sum(), 572)  # self.assertEqual(dm.sum(), 1160) #  new code version

    def test_compute_all_distances(self):
        nucleus_centroid = [6, 6]
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid, 11, 11)
        self.assertEqual(np.sum(dsAll).astype(int),
                         526)  # TODO : should be 507, but we keep the 1 offset for the V0 compatibility
        with self.assertRaises(IndexError):
            nucleus_centroid = [6, 16]
            ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid, 11, 31)

        nucleus_centroid = [2, 6]
        dsAll = ip.compute_all_distances_to_nucleus_centroid(nucleus_centroid, 11, 11)
        self.assertEqual(np.sum(dsAll).astype(int), 602)

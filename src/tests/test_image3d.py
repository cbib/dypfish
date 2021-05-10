#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

import numpy as np
from loguru import logger

import constants
import helpers
import image_processing as ip
import path
from image3d import Image3d
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)  # TODO this is annoying


class TestImage3d(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image3d(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_compute_cell_mask_slices(self):
        cell_mask = helpers.unit_circle(15, 7)
        slice1 = helpers.unit_circle(15, 5)
        slice2 = helpers.unit_circle(15, 3)
        slice3 = helpers.unit_circle(15, 1)
        height_map = cell_mask.copy()
        height_map[slice1 == 1] = 2
        height_map[slice2 == 1] = 3
        height_map[slice3 == 1] = 4
        zero_level = 3
        slices = self.img.compute_cell_mask_slices(height_map=height_map,
                                                   zero_level=zero_level, image_width=15, image_height=15)
        self.assertEqual(slices[:, :, 0][7, 7], 1)  # the center is always 1
        self.assertEqual(slices[:, :, 1][7, 7], 1)
        self.assertEqual(slices[:, :, 0].sum(), cell_mask.sum()) # bottom slice is the same as cell mask
        self.assertEqual(slices[:, :, 3].sum(), 1)  # the slice at the top the smallest
        self.assertGreater(slices[:, :, 0].sum(), slices[:, :, 1].sum())  # they grow in size from the bottom
        self.assertEqual(slices.sum(), height_map.sum())
        slices = self.img.compute_cell_mask_slices()
        self.assertEqual(slices[:, :, 0].sum(), self.img.get_cell_mask().sum()) # bottom slice same as cell_mask
        self.assertEqual(slices.shape, (512, 512, 12))
        self.assertEqual(slices.sum(), 426034.0)

    def test_get_zero_level(self):
        logger.warning("get_zero_level should return an int, to be cheched thoughout the code")
        self.assertEqual(self.img.get_zero_level(), 11)

    def test_compute_cytoplasmic_volume(self):
        cyt_vol = self.img.compute_cytoplasmic_volume()
        self.assertAlmostEqual(cyt_vol, 1206.276923076, places=5)
        self.assertAlmostEqual(cyt_vol,
                               self.img.compute_cell_volume() - self.img.compute_nucleus_volume()(), places=5)

    def test_compute_cell_volume(self):
        cell_vol = self.img.get_cell_volume()
        self.assertAlmostEqual(cell_vol, 1325.570019723, places=5)

    def test_compute_nucleus_volume(self):
        cell_vol = self.img.compute_nucleus_volume()
        self.assertAlmostEqual(cell_vol, 371.528205128, places=5)

    def test_adjust_height_map(self):
        height_map_original = self.img.get_height_map()
        height_map = self.img.adjust_height_map()
        self.assertGreater(height_map_original.sum(), height_map.sum())
        adj_height_map = self.img.adjust_height_map(cytoplasm=True)
        self.assertGreater(height_map.sum(), adj_height_map.sum())

    def test_compute_median_cytoplasmic_distance_from_nucleus(self):
        nucleus_centroid = self.img.get_nucleus_centroid()
        height_map = self.img.get_cytoplasm_height_map()
        dsAll = ip.compute_all_distances_to_nucleus_centroid3d(height_map, nucleus_centroid)
        median, max = self.img.compute_median_cytoplasmic_distance_from_nucleus3d(dsAll)
        self.assertAlmostEqual(median, 123.84264209, places=5)
        self.assertAlmostEqual(max, 212.492352803, places=5)

    def test_compute_cell_volume(self):
        volume = self.img.compute_cell_volume()
        self.assertAlmostEqual(volume, 1325.570019723, places=5)

    def test_compute_nucleus_volume(self):
        volume = self.img.compute_nucleus_volume()
        self.assertAlmostEqual(volume, 314.370019723, places=5)

    def test_compute_cytoplasmic_volume(self):
        volume = self.img.compute_cytoplasmic_volume()
        self.assertAlmostEqual(volume, 1011.19999999, places=5)
        self.assertAlmostEqual(volume, self.img.compute_cell_volume() - self.img.compute_nucleus_volume(), places=5)

    def test_compute_peripheral_cell_volume(self):
        threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        volume = self.img.compute_peripheral_cell_volume(threshold)
        self.assertAlmostEqual(volume, 133.3096646942, places=5)

    def test_compute_volumes_from_periphery(self):
        volumes = self.img.compute_volumes_from_periphery()
        self.assertEqual(len(volumes), 100)
        self.assertTrue(all(volumes[i] <= volumes[i + 1] for i in range(len(volumes) - 1)))
        self.assertAlmostEqual(volumes.sum(), 43899.15581854, places=5)
        self.assertAlmostEqual(volumes[99], self.img.compute_cytoplasmic_volume(), places=5)

    def test_compute_cytoplasmic_coordinates_peripheral_distance(self):
        coordinates = np.array([(20, 20, 10), (178, 178, 6)])
        distances = self.img.compute_cytoplasmic_coordinates_peripheral_distance(coordinates)
        self.assertTrue(np.isnan(distances[0]))
        self.assertEqual(distances[1], 29)

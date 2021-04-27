#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import numpy as np
import constants
import helpers
from repository import H5Repository, H5RepositoryWithCheckpoint
from image import Image, ImageWithMTOC
import path

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImageWithMTOC(TestCase):

    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/1")
        self.img1 = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/1")
        self.img10 = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/10")
        self.img11 = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/11")
        self.img12 = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/12")
        self.img13 = ImageWithMTOC(repository=self.repo, image_path="mrna/arhgdia/2h/13")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_is_a(self):
        self.assertTrue(self.img.is_a(self.repo, self.img._path))

    def test_mtoc_is_in_leading_edge(self):
        result1 = self.img.mtoc_is_in_leading_edge()
        result10 = self.img10.mtoc_is_in_leading_edge()
        result11 = self.img11.mtoc_is_in_leading_edge()
        result12 = self.img12.mtoc_is_in_leading_edge()
        result13 = self.img13.mtoc_is_in_leading_edge()
        self.assertEqual(result1, 1)
        self.assertEqual(result10, 0)
        self.assertEqual(result11, 1)
        self.assertEqual(result12, 0)
        self.assertEqual(result13, 1)

    def test_get_mtoc_quadrant(self):
        self.assertTrue(self.img1.mtoc_is_in_leading_edge())
        self.assertFalse(self.img10.mtoc_is_in_leading_edge())

    def test_get_mtoc_position(self):
        cell_mask = self.img.get_cell_mask()
        mtoc_position = self.img.get_mtoc_position()
        self.assertEqual(mtoc_position.tolist(), [265, 200])

    def test_get_quadrant_mask(self):
        cell_mask = helpers.unit_circle(11, 5)
        nucleus_centroid = [5, 5]
        mtoc_position = [5, 2]

        quadrant_mask = self.img.compute_quadrant_mask(degree=45, slices_num=4, nucleus_centroid=nucleus_centroid,
                                                       image_width=11, image_height=11, cell_mask=cell_mask,
                                                       mtoc_position=mtoc_position)
        self.assertEqual(quadrant_mask.sum(), 179)

        quadrant_mask = self.img.compute_quadrant_mask(degree=45, slices_num=4)
        self.assertEqual(quadrant_mask.sum(), 133618)


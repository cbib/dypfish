#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski


from unittest import TestCase

from loguru import logger

import constants
import path

constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImageMultiNucleus(TestCase):
    def test_get_multiple_nucleus_centroid(self):
        logger.error("this function is not tested")
        self.fail()

    def test_compute_nucleus_area(self):
        logger.error("this function is not tested")
        self.fail()

    def test_compute_cell_area(self):
        logger.error("this function is not tested")
        self.fail()


class TestimageWithSpotsAndZlines(TestCase):
    def test_get_z_lines_masks(self):
        logger.error("this function is not tested")
        self.fail()

    def test_compute_minimal_z_line_distance(self):
        logger.error("this function is not tested")
        self.fail()

    def test_get_z_lines_masks(self):
        logger.error("this function is not tested")
        self.fail()

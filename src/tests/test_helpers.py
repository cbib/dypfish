#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

from unittest import TestCase

import helpers
import mpi_calculator
import constants
import path
import numpy as np

constants.init_config(analysis_config_js_path=path.test_config_path)


class Test(TestCase):
    def test_compute_statistics_random_h_star(self):
        a = np.array([[0.51, 0.13, 0.7, 0.01, 1.9],
                      [0.4, 1.54, 0.3, 0.2, 0.8],
                      [0.67, 0.93, 1.76, 0.06, 1.2]])
        s5, s50, s95 = helpers.compute_statistics_random_h_star(h_sim=a, max_cell_radius=5, simulation_number=3)
        self.assertAlmostEqual(s5.sum(), -13.103848980244447)
        self.assertAlmostEqual(s50.sum(), -12.445951075315389)
        self.assertAlmostEqual(s95.sum(), -11.860677638869003)

    def test_rotate_point(self):
        nucleus_centroid = [0, 0]
        mtoc_position = [1, 1]
        degree = 45
        result = helpers.rotate_point(self, nucleus_centroid, mtoc_position, degree)
        self.assertEqual(result[0], 0)
        self.assertEqual(result[1], 1)

    def test_slope_from_points(self):
        point1 = np.array([0, 0])
        point2 = np.array([1, 1])
        result = helpers.slope_from_points(self, point1, point2)
        self.assertEqual(result, 1.0)
        self.assertEqual(np.degrees(np.arctan(1)), 45.0)

        rotate_point = np.array([268.0, 202.0])
        nucleus_centroid = np.array([252.98498233215548, 219.89697438162545])
        result = helpers.slope_from_points(self, nucleus_centroid, rotate_point)
        self.assertEqual(result, -0.8389696128335637)

    def test_rotate_meshgrid(self):
        radians = np.arctan(1)
        xx, yy = np.meshgrid(np.array(range(0, 3)) - 0, np.array(range(0, 3)) - 0)
        result_xx, result_yy = helpers.rotate_meshgrid(xx, yy, -radians)
        self.assertAlmostEqual(result_xx.sum(), 1.1102230246251565e-15)
        self.assertAlmostEqual(result_yy.sum(), 12.727922061357855)

        radians = np.arctan(-0.8389696128335637)
        xx, yy = np.meshgrid(np.array(range(0, 512)) - 252.98498233215548, np.array(range(0, 512)) - 219.89697438162545)
        result_xx, result_yy = helpers.rotate_meshgrid(xx, yy, -radians)
        self.assertEqual(np.shape(result_xx), (512, 512))
        self.assertEqual(np.shape(result_yy), (512, 512))

    def test_calculate_mpi(self):
        mpi = mpi_calculator.calculate_mpi([1, 1, 1, 1, 1, 1, 1], [7, 7, 7, 7, 7, 7, 7])
        self.assertEqual(mpi.index, -1.0)
        self.assertEqual(mpi.pvalue, 0.00041247833638552163)

    def test_median_confidence_interval(self):
        a = np.array([24, 38, 61, 22, 16, 57, 31, 29, 35])
        l, h = helpers.median_confidence_interval(a, cutoff=0.8)
        self.assertEqual(l, 29)
        self.assertEqual(h, 57)

    def test_sem(self):
        a = np.array([0, 0.1, 0, 3, 0.001, 10])
        sem = helpers.sem(a)
        self.assertAlmostEqual(sem, 1.49447, places=5)

    def test_compute_entropy(self):
        center = np.array([5, 5, 5])
        radius = 3
        points1 = helpers.random_points_in_sphere(center, radius, 400)
        entropy1 = helpers.compute_entropy(points1, k=15, norm='euclidean')
        radius = 5
        points2 = helpers.random_points_in_sphere(center, radius, 400)
        entropy2 = helpers.compute_entropy(points2, k=15, norm='euclidean')
        self.assertLess(entropy1, entropy2) # increases with volume

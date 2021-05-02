#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image3d import Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC
from image_set import ImageSet
import path
import constants

constants.init_config(analysis_config_js_path=path.test_config_path)  # TODO this is annoying


class TestImageSet(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)

    def tearDown(self) -> None:
        self.repo.clear()

    def test_can_build_mono_type_fish(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        for i in image_set.images:
            self.assertTrue(type(i) == Image3dWithSpotsAndMTOC)
            self.assertEqual(i.get_spots().shape[1], 3)

    def test_can_build_mono_type_if(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/3h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        for i in image_set.images:
            self.assertTrue(type(i) == Image3dWithIntensitiesAndMTOC)
            with self.assertRaises(AttributeError):  # 'ImageWithIntensities' object has no attribute 'get_spots'
                self.assertEqual(i.get_spots().shape[1], 3)
            self.assertEqual(i.get_intensities().shape, (512, 512))

    # TODO : major bottleneck identified in draw polygon
    def test_compute_spots_fractions_per_periphery(self):
        #self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/'])
        self.assertEqual(len(image_set.images), 20, "Expected 20 images")
        peripheral_fractions = image_set.compute_spots_fractions_per_periphery()
        self.assertEqual(peripheral_fractions.shape, (20, 100))
        self.assertAlmostEqual(peripheral_fractions.sum(), 1985.585078961, places=5)

    def test_compute_intensities_fractions_from_periphery(self):
        #self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/'])
        self.assertEqual(len(image_set.images), 20, "Expected 20 images")
        peripheral_fractions = image_set.compute_intensities_fractions_from_periphery()
        self.assertEqual(peripheral_fractions.shape, (20, 100))
        self.assertAlmostEqual(peripheral_fractions.sum(), 947.7339075841, places=5)

    def test_compute_cytoplasmic_spots_counts(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        spots_counts = image_set.compute_cytoplasmic_spots_counts()
        self.assertEqual(sorted(spots_counts), [76, 81, 81, 101, 155])

    def test_compute_degree_of_clustering(self):
        # TODO : raises RuntimeWarning: divide by zero encountered in true_divide
        np.random.seed(0)
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        clustering_indices = image_set.compute_degree_of_clustering()
        self.assertGreater(clustering_indices,
                           [0, 0, 0, 0, 0])  # TODO : how to better test this: np.random.seed(0) does not seem to work
        self.assertEqual(len(clustering_indices), 5)

    def test_compute_normalised_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        expected = np.array([[1.74626861, 1.], [0.64688921, 0.], [0.65049158, 0.], [1.54168467, 0.],
                             [1.14561588, 1.], [2.77314831, 0.], [1.26125903, 0.], [0.50137217, 0.],
                             [1.38473332, 1.], [0.51311625, 0.], [1.26063698, 0.], [1.02385656, 0.],
                             [1.17236082, 1.], [0.4422494, 0.], [2.20267986, 0.], [1.69777197, 0.],
                             [1.08763936, 1.], [3.26577748, 0.], [2.42107527, 0.], [0.46270397, 0.]])

        self.assertEqual(res.shape[0], expected.shape[0])
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        expected_mtoc_density = expected[expected[:, 1] == 1].sum() / len(expected[expected[:, 1] == 1])
        self.assertAlmostEqual(mtoc_density, expected_mtoc_density)
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        expected_non_mtoc_density = expected[expected[:, 1] == 0].sum() / len(expected[expected[:, 1] == 0])
        self.assertAlmostEqual(non_mtoc_density, expected_non_mtoc_density)

    def test_compute_normalised_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(mtoc_density, 2.1244031950903275, places=3)
        self.assertAlmostEqual(non_mtoc_density, 1.244101903437526, places=3)

    def test_mtoc_is_in_leading_edge(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.mtoc_is_in_leading_edge()
        self.assertListEqual(list(np.sort(result)), [False, False, True, True, True])

    def test_compute_cytoplasmic_spots_spread(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_spread())
        self.assertAlmostEqual(np.sum(result), 3.7360605912, places=5)
        self.assertAlmostEqual(result[2], 0.75042495387, places=5)

    def test_compute_spots_cytoplasmic_centrality(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_centrality())
        self.assertAlmostEqual(np.sum(result), 4.445184544252, places=5)
        self.assertAlmostEqual(result[2], 0.86929413285, places=5)

    def test_compute_intensities_cytoplasmic_centrality(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = np.sort(image_set.compute_intensities_cytoplasmic_centrality())
        self.assertAlmostEqual(np.sum(result), 3.35035785762, places=5)
        self.assertAlmostEqual(result[1], 0.60617975966, places=5)

    def test_compute_intensities_cytoplasmic_spread(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = np.sort(image_set.compute_intensities_cytoplasmic_spread())
        self.assertAlmostEqual(np.sum(result), 3.38954627718, places=5)
        self.assertAlmostEqual(result[1], 0.6615370616, places=5)

    def test_compute_surface_corrected_nm(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_surface_corrected_nm()
        self.assertAlmostEqual(result, 0.10109035262, places=5)

    def test_compute_normalized_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images, 1.1983781428,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7), 1.38104578,
                               places=5)  # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.218815736,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.371284986, places=5)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3,
                                                                  stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.0011051805,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 0.8443642562, places=5)

    def test_compute_normalized_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images, 1.0387081758,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7), 1.2690913,
                               places=5)  # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.1715709,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.7474789, places=5)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3,
                                                                  stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 3.56607132,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 4.35721309, places=5)

    def test_compute_zline_distance(self):
        self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/actn2/immature/'])
        result = image_set.compute_zline_distance(20)
        test = [[0.24374599, 0.03463759, 0.0365619, 0.0436177, 0.03207184,
                 0.02758178, 0.0365619, 0.03207184, 0.03014753, 0.02694035,
                 0.02309173, 0.01860167, 0.02180885, 0.02758178, 0.01988454,
                 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.38974359, 0.01230769, 0.01025641, 0.02871795, 0.02153846,
                 0.02461538, 0.01948718, 0.02666667, 0.03076923, 0.01333333,
                 0.0225641, 0.01025641, 0.01538462, 0.01333333, 0.00923077,
                 0.0, 0.0, 0.0, 0.0, 0.0]]
        self.assertEqual(len(result), len(test))
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=5)


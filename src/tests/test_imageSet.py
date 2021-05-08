#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase

import numpy as np
from loguru import logger

import constants
import path
from image3dWithIntensities import Image3dWithIntensitiesAndMTOC
from image3dWithSpots import Image3dWithSpotsAndMTOC
from image_set import ImageSet
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=path.test_config_path)


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

    # major bottleneck identified in draw polygon in 2D
    def test_compute_spots_signal_from_periphery(self):
        # self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/'])
        self.assertEqual(len(image_set.images), 20, "Expected 20 images")
        peripheral_signals = image_set.compute_signal_from_periphery()
        self.assertEqual(peripheral_signals.shape, (20, 100))
        self.assertAlmostEqual(peripheral_signals.sum(), 119196.0)

    def test_compute_areas_from_periphery(self):
        # self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/'])
        self.assertEqual(len(image_set.images), 20, "Expected 20 images")
        peripheral_areas = image_set.compute_areas_from_periphery()
        self.assertEqual(peripheral_areas.shape, (20, 100))
        self.assertAlmostEqual(peripheral_areas.sum(), 429444.207758053, places=5)

    def test_compute_intensities_signal_from_periphery(self):
        # self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/'])
        self.assertEqual(len(image_set.images), 20, "Expected 20 images")
        peripheral_signals = image_set.compute_signal_from_periphery()
        self.assertEqual(peripheral_signals.shape, (20, 100))
        self.assertAlmostEqual(peripheral_signals.sum(), 426775759669.0)

    def test_compute_cytoplsamic_spots_fractions_per_periphery(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/'])
        peripheral_fractions = image_set.compute_cytoplsamic_spots_fractions_per_periphery()
        self.assertEqual(peripheral_fractions.shape, (20, 100))
        self.assertAlmostEqual(peripheral_fractions.sum(), 1909.56435359227, places=5)
        self.assertTrue(np.all(peripheral_fractions[:, 99] == 1))

    def test_compute_cytoplsamic_intensities_fractions_per_periphery(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/'])
        peripheral_fractions = image_set.compute_cytoplsamic_intensities_fractions_per_periphery()
        self.assertEqual(peripheral_fractions.shape, (20, 100))
        self.assertAlmostEqual(peripheral_fractions.sum(), 1076.1379240779, places=5)
        self.assertTrue(np.all(peripheral_fractions[:, 99] == 1))

    def test_compute_cytoplasmic_spots_counts(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        spots_counts = image_set.compute_cytoplasmic_spots_counts()
        self.assertEqual(sorted(spots_counts), [33, 62, 64, 68, 151])

    def test_compute_degree_of_clustering(self):
        logger.error("This function has not been tested with cytoplasmic spots and new random spots")
        self.fail()
        self.skipTest() # skipping because not tested, see above
        np.random.seed(0)
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        clustering_indices = image_set.compute_degree_of_clustering()
        self.assertGreater(clustering_indices,
                           [0, 0, 0, 0, 0])  # TODO : how to better test this: np.random.seed(0) does not seem to work
        self.assertEqual(len(clustering_indices), 5)

    def test_compute_normalised_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        self.assertEqual(res.shape[0], 20)
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        self.assertAlmostEqual(mtoc_density, 2.21729906014287, places=5)
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(non_mtoc_density, 1.1013546485392969, places=5)

    def test_compute_normalised_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(mtoc_density, 1.9772754505021, places=3)
        self.assertAlmostEqual(non_mtoc_density, 1.0734285610952, places=3)

    def test_mtoc_is_in_leading_edge(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.mtoc_is_in_leading_edge()
        self.assertListEqual(list(np.sort(result)), [False, False, True, True, True])

    def test_compute_cytoplasmic_spots_spread(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_spread())
        self.assertAlmostEqual(np.sum(result), 84.32560104311, places=5)
        self.assertAlmostEqual(result[2], 16.983566893, places=5)

    def test_compute_spots_cytoplasmic_centrality(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_centrality())
        self.assertAlmostEqual(np.sum(result), 2.9176, places=3)
        self.assertAlmostEqual(result[2], 0.90269627884, places=5)

    def test_compute_intensities_normalized_distance_to_nucleus(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = np.sort(image_set.compute_intensities_normalized_distance_to_nucleus())
        self.assertAlmostEqual(np.sum(result), 4.97198077749, places=5)
        self.assertAlmostEqual(result[1], 0.92969214548, places=5)

    def test_compute_intensities_cytoplasmic_spread(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = np.sort(image_set.compute_intensities_cytoplasmic_spread())
        self.assertAlmostEqual(np.sum(result), 2.739389070456, places=5)
        self.assertAlmostEqual(result[1], 0.6615370616, places=5)

    def test_compute_surface_corrected_nm(self):
        logger.error("child class functions not tested?")
        self.fail()
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_surface_corrected_nm()
        self.assertAlmostEqual(result, 0.089117434154, places=5)

    def test_compute_normalized_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images, 1.074792686952,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7), 1.1392495492,
                               places=5)  # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.0146171384,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.093835714152, places=5)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3,
                                                                  stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 0.8149847897169,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 0.60899299769, places=5)

    def test_compute_normalized_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images, 0.919149036299,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7), 1.095215502,
                               places=5)  # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.028008679,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.51199182,
                               places=5)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3,
                                                                  stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 3.036024433241,
                               places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 3.7152829163, places=5)

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

    def test_compute_cell_mask_between_nucleus_centroids(self):
        logger.error("This function has not been tested")
        self.fail()

    def test_compute_volume_corrected_nm(self):
        logger.error("This function has not been tested")
        self.fail()

    def test_compute_spots_peripheral_distance(self):
        logger.error("This function has not been tested")
        self.fail()

    def test_compute_cytoplasmic_density(self):
        logger.error("This function has not been tested")
        self.fail()

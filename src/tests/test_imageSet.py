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
        for image in image_set.images:
            self.assertTrue(type(image) == Image3dWithIntensitiesAndMTOC)
            with self.assertRaises(AttributeError):  # 'ImageWithIntensities' object has no attribute 'get_spots'
                self.assertEqual(image.get_spots().shape[1], 3)
            self.assertEqual(image.get_intensities().shape, (512, 512))

    def test_compute_areas_from_periphery(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/3h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        peripheral_areas = image_set.compute_areas_from_periphery()
        self.assertEqual(peripheral_areas.shape, (5, 100))
        cytoplasmic_areas = [img.compute_cell_area() - img.compute_nucleus_area() for img in image_set.images]
        self.assertTrue(np.allclose(peripheral_areas[:, 99], np.array(cytoplasmic_areas)))
        self.assertAlmostEqual(peripheral_areas.sum(), 142574.96909927676, places=5)

    def test_compute_volumes_from_periphery(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/3h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        peripheral_volumes = image_set.compute_volumes_from_periphery()
        self.assertEqual(peripheral_volumes.shape, (5, 100))
        cytoplasmic_volumes = [img.compute_cytoplasmic_volume() for img in image_set.images]
        self.assertTrue(np.allclose(peripheral_volumes[:, 99], np.array(cytoplasmic_volumes)))
        self.assertAlmostEqual(peripheral_volumes.sum(), 121874.96804733, places=5)

    # major bottleneck identified in draw polygon in 2D
    def test_compute_spots_signal_from_periphery(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        peripheral_signals = image_set.compute_signal_from_periphery()
        self.assertEqual(peripheral_signals.shape, (5, 100))
        spots_counts = [img.compute_cytoplasmic_total_spots() for img in image_set.images]
        self.assertTrue(np.all(peripheral_signals[:, 99] == np.array(spots_counts)))
        self.assertAlmostEqual(peripheral_signals.sum(), 20166.0)

    def test_compute_intensities_signal_from_periphery(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/4h/'])
        self.assertEqual(len(image_set.images), 5, "Expected 5 images")
        peripheral_signals = image_set.compute_signal_from_periphery()
        cytoplasmic_intensities = [img.compute_cytoplasmic_total_intensity() for img in image_set.images]
        self.assertEqual(peripheral_signals.shape, (5, 100))
        self.assertTrue(np.allclose(peripheral_signals[:, 99], np.array(cytoplasmic_intensities)))
        self.assertAlmostEqual(peripheral_signals.sum(), 29384734352.0)

    def test_compute_cytoplasmic_spots_fractions_per_periphery(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/3h/'])
        peripheral_fractions = image_set.compute_cytoplasmic_spots_fractions_per_periphery()
        self.assertEqual(peripheral_fractions.shape, (5, 100))
        self.assertTrue(np.all(peripheral_fractions[:, 99] == 1))
        self.assertAlmostEqual(peripheral_fractions.sum(), 719.995179047, places=5)

        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/4h/'], force2D=True)
        peripheral_fractions = image_set.compute_cytoplasmic_spots_fractions_per_periphery(force2D=True)
        self.assertEqual(peripheral_fractions.shape, (5, 100))
        self.assertTrue(np.all(peripheral_fractions[:, 99] == 1))

    def test_compute_cytoplasmic_intensities_fractions_per_periphery(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/4h/'])
        peripheral_fractions = image_set.compute_cytoplasmic_intensities_fractions_per_periphery()
        self.assertEqual(peripheral_fractions.shape, (5, 100))
        self.assertTrue(np.all(peripheral_fractions[:, 99] == 1))
        self.assertAlmostEqual(peripheral_fractions.sum(), 126.65345873037, places=5)

    def test_compute_cytoplasmic_spots_counts(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        spots_counts = image_set.compute_cytoplasmic_spots_counts()
        self.assertEqual(sorted(spots_counts), [33, 68, 71, 73, 153])

    def test_compute_degree_of_clustering(self):
        # logger.error("This function has not been tested with cytoplasmic spots and new random spots")
        # self.fail()
        self.skipTest()  # skipping because not tested, see above
        np.random.seed(0)
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        clustering_indices = image_set.compute_degree_of_clustering()
        self.assertEqual(len(clustering_indices), 5)
        total_sum = 0
        for ci in clustering_indices:
            total_sum += np.sum(ci)
        print(total_sum)
        self.assertAlmostEqual(total_sum, 5012.542842389925, places=5)
        #self.assertGreater(clustering_indices,[0, 0, 0, 0, 0])  # TODO : how to better test this: np.random.seed(0) does not seem to work


    def test_compute_normalised_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        self.assertEqual(res.shape[0], 20)
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        self.assertAlmostEqual(mtoc_density, 2.22758047513, places=5)
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(non_mtoc_density, 1.089256691215, places=5)

    def test_compute_normalised_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(mtoc_density, 1.98375276421, places=3)
        self.assertAlmostEqual(non_mtoc_density, 1.0445809579196, places=3)

    def test_mtoc_is_in_leading_edge(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.mtoc_is_in_leading_edge()
        self.assertListEqual(list(np.sort(result)), [False, False, True, True, True])

    def test_compute_cytoplasmic_spots_spread(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_spread())
        self.assertAlmostEqual(np.sum(result), 84.307443952, places=5)
        self.assertAlmostEqual(result[2], 16.9820520014, places=5)

    def test_compute_spots_cytoplasmic_centrality(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = np.sort(image_set.compute_cytoplasmic_spots_centrality())
        self.assertAlmostEqual(np.sum(result), 2.325, places=3)
        self.assertAlmostEqual(result[3], 0.51, places=5)

    def test_compute_cytoplasmic_intensities_centrality(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/3h/'])
        result = np.sort(image_set.compute_cytoplasmic_intensities_centrality())
        self.assertAlmostEqual(np.sum(result), 3.89, places=2)
        self.assertAlmostEqual(result[1], 0.73, places=5)

    def test_compute_intensities_cytoplasmic_spread(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = np.sort(image_set.compute_intensities_cytoplasmic_spread())
        self.assertAlmostEqual(result.sum(), 3.0812597654, places=5)
        self.assertAlmostEqual(result[1], 0.611584381, places=5)

    def test_compute_surface_corrected_nm(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_surface_corrected_nm()
        self.assertAlmostEqual(result, 0.1010903526279469, places=5)

    def test_compute_normalized_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = image_set.__sizeof__()
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images,
                               1.0780173273, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7),
                               1.1263246459, places=5)  # non MTOC quadrant density
        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images),
                               1.0225376792, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3),
                               1.08160649427, places=5)  # non MTOC quadrant density
        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True,
                                                                  stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images),
                               0.79154809935, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3),
                               0.59408432134, places=5)  # non MTOC quadrant density

    def test_compute_normalized_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/3h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = image_set.__sizeof__()
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:, 1] == 1][:, 0].sum() / num_images,
                               1.098770606, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:, 1] == 0][:, 0].sum() / (num_images * 7),
                               0.865466580, places=5)  # non MTOC quadtant density
        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3 * num_images),
                               0.9803183396, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3),
                               0.6801687989, places=5)  # non MTOC quadtant density
        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True,
                                                                  stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images),
                               0.126123087, places=5)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3),
                               0.011356787, places=5)  # non MTOC quadtant density

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
        image_set = ImageSet(self.repo, path_list=['mrna/actn2/immature/'])
        nuc_dist, nucs_dist, cell_masks, nucs_pos = image_set.compute_cell_mask_between_nucleus_centroids()
        self.assertEqual(np.sorted([754, 483, 526]), np.sorted(nuc_dist))
        self.assertEqual([[754], [483, 526]], nucs_dist)
        self.assertEqual([[[110, 864]], [[154, 637], [637, 1163]]], nucs_pos)

        #logger.error("This function has not been tested")
        #self.fail()

    def test_compute_volume_corrected_nm(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        vcnm = image_set.compute_volume_corrected_nm()
        self.assertAlmostEqual(0.029713582957013474, vcnm, places=5)


    def test_compute_spots_peripheral_distance(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_spots_peripheral_distance()
        total_sum=0
        for res in result:
            total_sum += np.sum(res)
        self.assertEqual(total_sum, 20030.0)

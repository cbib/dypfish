#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from unittest import TestCase
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image3d import Image3dWithIntensities, Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC
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
    def test_compute_histogram_spots_peripheral_counts_perf(self):
        # self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        peripheral_fractions = image_set.compute_histogram_spots_peripheral_counts()
        test = [0.2679738562091503, 0.2191780821917808, 0.12121212121212122, 0.3380281690140845, 0.16176470588235295]
        self.assertEqual(sorted(peripheral_fractions), sorted(test))

    # this can be long : ~5min
    def test_compute_histogram_spots_peripheral_counts(self):
        self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['mrna/'])
        self.assertEqual(len(image_set.images), 40, "Expected 40 images")
        peripheral_fractions = image_set.compute_histogram_spots_peripheral_counts()
        self.assertEqual(len(peripheral_fractions), 40, "Expected 40 fractions")
        # TODO : test the sum?

    def test_compute_histogram_intensities_peripheral_fractions(self):
        self.skipTest("Skipping for inefficiency reasons")
        image_set = ImageSet(self.repo, path_list=['protein/'])
        self.assertEqual(len(image_set.images), 40, "Expected 40 images")
        peripheral_fractions = image_set.compute_histogram_intensities_peripheral_fractions()
        self.assertEqual(len(peripheral_fractions), 40, "Expected 40 fraction arrays")
        # TODO : test the sum?

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

    def test_compute_mtoc_dependent_degree_of_clustering(self):
        # TODO : raises RuntimeWarning: divide by zero encountered in true_divide
        np.random.seed(0)
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        clustering_indices = image_set.compute_mtoc_dependent_degree_of_clustering()
        self.assertEqual(len(clustering_indices), 3)

    def test_compute_normalised_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        expected = np.array([[1.37693091, 1.],
                    [0.50058693, 0.],
                    [1.23891387, 0.],
                    [1.00759389, 0.],
                    [1.12683717, 1.],
                    [2.72051787, 0.],
                    [1.23250137, 0.],
                    [0.49768729, 0.],
                    [1.15910767, 1.],
                    [0.43691438, 0.],
                    [2.15659572, 0.],
                    [1.64046996, 0.],
                    [1.72992133, 1.],
                    [0.63856315, 0.],
                    [0.63143629, 0.],
                    [1.48932139, 0.],
                    [1.07688498, 1.],
                    [3.12587233, 0.],
                    [2.36768991, 0.],
                    [0.45863156, 0.]])

        self.assertEqual(res.shape[0], expected.shape[0])
        mtoc_density = res[res[:,1]==1].sum() / len(res[res[:,1]==1])
        expected_mtoc_density = expected[expected[:,1]==1].sum() / len(expected[expected[:,1]==1])
        self.assertAlmostEqual(mtoc_density, expected_mtoc_density)
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        expected_non_mtoc_density = expected[expected[:, 1] == 0].sum() / len(expected[expected[:, 1] == 0])
        self.assertAlmostEqual(non_mtoc_density, expected_non_mtoc_density)

    def test_compute_normalised_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        res = image_set.compute_normalised_quadrant_densities()
        mtoc_density = res[res[:, 1] == 1].sum() / len(res[res[:, 1] == 1])
        non_mtoc_density = res[res[:, 1] == 0].sum() / len(res[res[:, 1] == 0])
        self.assertAlmostEqual(mtoc_density, 2.10689868198936, places=3)
        self.assertAlmostEqual(non_mtoc_density, 1.2117248932548452, places=3)

    def test_mtoc_is_in_leading_edge(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.mtoc_is_in_leading_edge()
        self.assertListEqual(list(np.sort(result)), [False, False, True, True, True])

    def test_compute_spots_cytoplasmic_spread(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_spots_cytoplasmic_spread()
        test = [0.9674088178799569, 0.8851708145804912, 0.9470329525805975, 0.9253082533496478, 1.079142761044456]
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=5)

    def test_compute_intensities_cytoplasmic_spread(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = image_set.compute_intensities_cytoplasmic_spread()
        test = [2.5397563772060727, 2.5923623840262615, 2.586344214162044, 2.291650989581199, 2.6410827964730803]
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=5)

    def test_compute_surface_corrected_nm(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_surface_corrected_nm()
        self.assertAlmostEqual(result, 0.101090352627947, places=5)

    def test_compute_normalized_quadrant_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:,1]==1][:,0].sum() / num_images, 1.190650759626, places=3) # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:,1]==0][:,0].sum() / (num_images*7), 1.345688920915, places=3) # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3*num_images), 1.218815736, places=3)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.371284986, places=3)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 1.0011051805, places=3)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 0.8443642562, places=3)

    def test_compute_normalized_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result1 = image_set.compute_normalised_quadrant_densities(quadrants_num=8)
        num_images = int(result1.shape[0] / 8)
        self.assertEqual(num_images, image_set.__sizeof__())
        self.assertAlmostEqual(result1[result1[:,1]==1][:,0].sum() / num_images, 1.02719092608, places=3) # MTOC quadrant density
        self.assertAlmostEqual(result1[result1[:,1]==0][:,0].sum() / (num_images*7), 1.233073655, places=3) # non MTOC quadtant density

        result2 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result2[result2[:, 1] == 1][:, 0].sum() / (3*num_images), 1.171570955, places=3)  # MTOC quadrant density
        self.assertAlmostEqual(result2[result2[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 1.7474789442, places=3)

        result3 = image_set.compute_normalised_quadrant_densities(quadrants_num=8, peripheral_flag=True, stripes=3, stripes_flag=True)
        self.assertAlmostEqual(result3[result3[:, 1] == 1][:, 0].sum() / (3 * num_images), 3.5660713211, places=3)  # MTOC quadrant density
        self.assertAlmostEqual(result3[result3[:, 1] == 0][:, 0].sum() / (num_images * 7 * 3), 4.3572130927, places=3)


    def test_compute_zline_distance(self):
        image_set = ImageSet(self.repo, path_list=['mrna/actn2/immature/'])
        result = image_set.compute_zline_distance(20)
        test= [[0.24374599, 0.03463759, 0.0365619 , 0.0436177 , 0.03207184,
                0.02758178, 0.0365619 , 0.03207184, 0.03014753, 0.02694035,
                0.02309173, 0.01860167, 0.02180885, 0.02758178, 0.01988454,
                0.0, 0.0, 0.0, 0.0, 0.0],
               [0.38974359, 0.01230769, 0.01025641, 0.02871795, 0.02153846,
                0.02461538, 0.01948718, 0.02666667, 0.03076923, 0.01333333,
                0.0225641 , 0.01025641, 0.01538462, 0.01333333, 0.00923077,
                0.0, 0.0, 0.0, 0.0, 0.0 ]]
        self.assertEqual(len(result), len(test))
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=5)


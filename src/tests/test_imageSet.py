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
        dict = image_set.compute_normalised_quadrant_densities(mtoc_quadrant_label='MTOC',
                                                               quadrant_labels=['Non MTOC1', 'Non MTOC2', 'Non MTOC3'])

        dict_to_test = {'Non MTOC1': [0.49768729218137575, 1.489321390402567, 1.0075938930481634, 1.6404699643896201, 0.45863155775318776],
                        'Non MTOC2': [1.2325013740513415, 0.6314362884578822, 1.2389138731555531, 2.1565957192717073, 2.367689909949143],
                        'Non MTOC3': [2.7205178674369273, 0.6385631502060634, 0.5005869269604428, 0.43691438495349, 3.1258723291183923],
                        'MTOC': [1.1268371735918787, 1.7299213344562772, 1.3769309057937469, 1.159107670776693, 1.0768849783625516]}

        self.assertEqual(len(dict), len(dict_to_test))
        self.assertAlmostEqual(np.sum(dict['MTOC']), np.sum(dict_to_test['MTOC']), places=2)
        self.assertAlmostEqual(np.sum(dict['Non MTOC1']), np.sum(dict_to_test['Non MTOC1']), places=3)
        self.assertAlmostEqual(np.sum(dict['Non MTOC2']), np.sum(dict_to_test['Non MTOC2']), places=3)
        self.assertAlmostEqual(np.sum(dict['Non MTOC3']), np.sum(dict_to_test['Non MTOC3']), places=3)

    def test_compute_normalised_quadrant_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        dict = image_set.compute_normalised_quadrant_densities(mtoc_quadrant_label='MTOC',
                                                               quadrant_labels=['Non MTOC1', 'Non MTOC2', 'Non MTOC3'])
        dict_to_test = {
            'Non MTOC1': [1.3771057445405632, 1.5614476627026483, 1.664727417509405, 1.0827726652063892,
                          1.229341850783425],
            'Non MTOC2': [1.085037753765736, 1.297011182965118, 1.4864117380673487, 1.0279187806700973,
                          1.312671364648296],
            'Non MTOC3': [1.011541659193154, 0.9375504850744625, 0.9207203335157421, 1.2291889509826397,
                          0.9526461385119073],
            'MTOC': [1.1147102515553686, 1.0748140925013452, 0.9524250785256029, 1.1768732783823252,
                     1.2157700467455073]}

        self.assertEqual(len(dict), len(dict_to_test))
        self.assertAlmostEqual(np.sum(dict['MTOC']), np.sum(dict_to_test['MTOC']), places=3)
        self.assertAlmostEqual(np.sum(dict['Non MTOC1']), np.sum(dict_to_test['Non MTOC1']), places=3)
        self.assertAlmostEqual(np.sum(dict['Non MTOC2']), np.sum(dict_to_test['Non MTOC2']), places=3)
        self.assertAlmostEqual(np.sum(dict['Non MTOC3']), np.sum(dict_to_test['Non MTOC3']), places=3)

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

    def test_compute_normalized_quadrant_and_slice_densities_mrna(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=8, stripes = 3)

        test=[[0.0, 0.0, 0.0, 0.0002607, 0.00053841, 0.0,
            0.0, 0.00215078, 0.00229417, 0.00113145, 0.00208612, 0.00149977,
            0.00401398, 0.00365251, 0.00426, 0.00269258, 0.00227383, 0.0,
            0.00295789, 0.00310227, 0.00802069, 0.01203337, 0.00661882, 0.0],
            [0.00179082, 0.0001572, 0.00068915, 0.0, 0.0, 0.0009817,
             0.0, 0.00045781, 0.00486055, 0.00398069, 0.00132457, 0.002896,
             0.002572, 0.00379238, 0.00125048, 0.00129015, 0.00286275, 0.00430695,
             0.00232102, 0.00550598, 0.00278462, 0.00844987, 0.00043668, 0.00613927],
            [0.0, 0.0, 0.0, 0.00070836, 0.00036358, 0.00041137,
             0, 0.00030849, 0.00258737, 0.00412314, 0.00210611, 0.00359072,
             0.00130157, 0.00229797, 0.00089768, 0.00127565, 0.00234351, 0.01404688,
            0.01005627, 0.00230535, 0.0035237, 0.00084205, 0.00090538, 0.00437455],
            [0.00056598, 0.00042134, 0.00026082, 0.0003623, 0.00047486, 0.0,
            0.00073601, 0.00261043, 0.00202661, 0.0023504, 0.00109015, 0.00254217,
            0.0020163, 0.00147113, 0.00673483, 0.0, 0.00873729, 0.00091848,
            0.00226075, 0.00628799, 0.00508859, 0.00518854, 0.00639555, 0.0042712],
            [0.00323102, 0.0, 0.0, 0.0, 0.00145618, 0.0011335,
            0.0, 0.0, 0.00439512, 0.00372517, 0.00078775, 0.00293404,
            0.00512146, 0.00609113, 0.00649155, 0.0, 0.00791599, 0.00812198,
            0.0, 0.00136125, 0.0, 0.00214355, 0.0054675, 0.0043093]]

        self.assertEqual(len(result), len(test))
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=3)

    def test_compute_normalized_quadrant_and_slice_densities_protein(self):
        image_set = ImageSet(self.repo, path_list=['protein/arhgdia/2h/'])
        result = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=8, stripes=3)
        self.assertAlmostEqual(result[2, 3], 0.11328982823)
        self.assertAlmostEqual(result.sum(), 16.109926486)

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


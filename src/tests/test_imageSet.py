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
        dict_to_test = {'Non MTOC1': [0.45863155775318776, 1.4090463349586326, 0.9506388740096401, 0.49768729218137575,
                                      1.6404699643896201],
                        'Non MTOC2': [2.36822678485314, 0.6286779012909666, 1.246944937763505, 1.233215694432939,
                                      2.15916541093823],
                        'Non MTOC3': [3.124796815655304, 0.6657557729895914, 0.6025062459115876, 2.71873528147916,
                                      0.4366406563223134],
                        'MTOC': [1.0768849783625516, 1.7501883153478572, 1.3718400011554692, 1.1268371735918787,
                                 1.159107670776693]}

        self.assertEqual(len(dict), len(dict_to_test))
        self.assertAlmostEqual(np.sum(dict['MTOC']), np.sum(dict_to_test['MTOC']), places=3)
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
        self.assertAlmostEquals(np.sum(result), np.sum(test), places=5)

    def test_compute_surface_corrected_nm(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_surface_corrected_nm()
        self.assertAlmostEquals(result, 0.101090352627947, places=5)


    def test_compute_normalized_quadrant_and_slice_densities(self):
        image_set = ImageSet(self.repo, path_list=['mrna/arhgdia/2h/'])
        result = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=8, stripes = 3)
        print(result)
        test=[[2.30111172e-04, 2.01989822e-05, 8.85521631e-05, 0.00000000e+00,
       0.00000000e+00, 1.26142488e-04, 0.00000000e+00, 5.88262094e-05,
       6.24553515e-04, 5.11497212e-04, 1.70200571e-04, 3.72120340e-04,
       3.30487499e-04, 4.87300445e-04, 1.60680029e-04, 1.65777539e-04,
       3.67847366e-04, 5.53418951e-04, 2.98238524e-04, 7.07488054e-04,
       3.57808019e-04, 1.08576227e-03, 5.61113604e-05, 7.88863248e-04],[0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 7.67194000e-05,
       3.93783744e-05, 4.45535598e-05, 0.00000000e+00, 3.34113383e-05,
       2.80227237e-04, 4.46559818e-04, 2.28104181e-04, 3.88895938e-04,
       1.40967722e-04, 2.48883435e-04, 9.72239844e-05, 1.38160399e-04,
       2.53815576e-04, 1.52135898e-03, 1.08915233e-03, 2.49683189e-04,
       3.81637565e-04, 9.11988361e-05, 9.80583144e-05, 4.73789075e-04],[5.26597973e-05, 3.92026993e-05, 2.42671356e-05, 3.37091786e-05,
       4.41815494e-05, 0.00000000e+00, 6.84800635e-05, 2.42879882e-04,
       1.88560015e-04, 2.18685851e-04, 1.01429951e-04, 2.36528420e-04,
       1.87600528e-04, 1.36877221e-04, 6.26622402e-04, 0.00000000e+00,
       8.12935622e-04, 8.54577361e-05, 2.10345274e-04, 5.85048305e-04,
       4.73453036e-04, 4.82752284e-04, 5.95055710e-04, 3.97401001e-04],[0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.25106322e-05,
       4.64901895e-05, 0.00000000e+00, 0.00000000e+00, 1.85712716e-04,
       1.98093564e-04, 9.76965049e-05, 1.80129436e-04, 1.29500443e-04,
       3.46593290e-04, 3.15381949e-04, 3.67836780e-04, 2.32495308e-04,
       1.96337698e-04, 0.00000000e+00, 2.55403651e-04, 2.67870584e-04,
       6.92559668e-04, 1.03904087e-03, 5.71513096e-04, 0.00000000e+00], [2.74582564e-04, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       1.23750882e-04, 9.63282644e-05, 0.00000000e+00, 0.00000000e+00,
       3.73511828e-04, 3.16577544e-04, 6.69452784e-05, 2.49344673e-04,
       4.35238466e-04, 5.17644829e-04, 5.51673848e-04, 0.00000000e+00,
       6.72727281e-04, 6.90233196e-04, 0.00000000e+00, 1.15683944e-04,
       0.00000000e+00, 1.82166033e-04, 4.64646470e-04, 3.66218402e-04]]

        self.assertEqual(len(result), len(test))
        self.assertAlmostEqual(np.sum(result), np.sum(test), places=3)

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


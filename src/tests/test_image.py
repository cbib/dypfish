import pathlib
from unittest import TestCase

import numpy as np
import constants
from repository import H5Repository, H5RepositoryWithCheckpoint
from image import Image
import path


class TestImage(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5Repository(repo_path=self.h5_sample_path)
        self.img = Image(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def test_get_cell_mask(self):
        cell_mask = self.img.get_cell_mask()
        self.assertEqual(cell_mask.shape, (512, 512))
        expected_values = [0, 1]
        self.assertTrue(np.any(np.isin(cell_mask, expected_values)))
        self.assertEqual(cell_mask.sum(), 57969)

    def test_get_nucleus_mask(self):
        nucleus_mask = self.img.get_nucleus_mask()
        self.assertEqual(nucleus_mask.shape, (512, 512))
        expected_values = [0, 1]
        self.assertTrue(np.any(np.isin(nucleus_mask, expected_values)))
        self.assertEqual(nucleus_mask.sum(), 9056)

    def test_compute_cytoplasm_mask(self):
        cytoplasm_mask = self.img.compute_cytoplasm_mask()
        self.assertEqual(cytoplasm_mask.shape, (512, 512))
        expected_values = [0, 1]
        self.assertTrue(np.any(np.isin(cytoplasm_mask, expected_values)))
        self.assertEqual(cytoplasm_mask.sum(), 48913)

    def test_get_nucleus_centroid(self):
        nucleus_centroid = self.img.get_nucleus_centroid()
        self.assertEqual(nucleus_centroid.tolist(), [253, 220])

    def test_get_cell_area(self):
        with self.assertRaises(PermissionError):
            # primary h5 repo should not allow write
            area = self.img.get_cell_area()


constants.init_config(analysis_config_js_path=path.test_config_path)


class TestImageSecondary(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(path.global_example_data, "basic.h5")
        self.repo = H5RepositoryWithCheckpoint(repo_path=self.h5_sample_path)
        self.img = Image(repository=self.repo, image_path="mrna/arhgdia/2h/1")

    def tearDown(self) -> None:
        self.repo.clear()

    def test_get_cell_area_write(self):
        area = self.img.get_cell_area()
        self.assertIsNotNone(area)
        # second call
        area = self.img.get_cell_area()
        self.assertIsNotNone(area)

    def test_get_cytoplasm_mask(self):
        self.assertEqual(self.img.get_cytoplasm_mask().sum(), 48913)

    def test_compute_peripheral_areas_in_cytoplasm(self):
        areas = self.img.compute_peripheral_areas()
        self.assertEqual(len(areas), 101)
        # test arbitrary values
        self.assertEqual(areas[2], 585.8251150558842)  # without nucleus : 490.5614727153188
        self.assertEqual(areas[4], 571.4766600920447)  # without nucleus : 476.21301775147924

    def test_is_in_cytoplasm(self):
        # test a random coord
        self.assertEqual(self.img.is_in_cytoplasm([50, 70]), False)
        # test some coords for coherency i.e. that they are not swapped (at the cell's border)
        self.assertEqual(self.img.is_in_cytoplasm([266, 147]), True)
        self.assertEqual(self.img.is_in_cytoplasm([147, 398]), True)
        # test some coords in the nucleus
        self.assertEqual(self.img.is_in_cytoplasm([266, 265]), False)

    def test_compute_cell_diameter(self):
        d = self.img.compute_cell_diameter()
        self.assertAlmostEqual(d, 311.192866242, places=5)

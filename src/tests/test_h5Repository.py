import pathlib
from unittest import TestCase

import h5py

from repository import H5Repository
from path import global_example_data


class TestH5Repository(TestCase):
    def setUp(self) -> None:
        self.h5_sample_path = pathlib.Path(global_example_data, "basic.h5")
        self.repo = H5Repository(repo_path=self.h5_sample_path)

    def test_image_is_present(self):
        self.assertTrue(self.repo.is_present("mrna/arhgdia/2h/1"), "Expected image to be found")
        self.assertFalse(self.repo.is_present("mrna/arhgdia/2h/11111"), "Expected image to be absent")

    def test_gene_is_present(self):
        self.assertTrue(self.repo.is_present("mrna/arhgdia"), "Expected gene to be found")

    def test_get_spots(self):
        self.assertEqual(type(self.repo.get("mrna/arhgdia/2h/1/spots")), h5py.Dataset)
        self.assertEqual(type(self.repo.get("mrna/arhgdia/2h/1")), h5py.Group)

    def test_get_image_list_inexistant(self):
        with self.assertRaises(LookupError):
            retrieved_paths = self.repo.get_image_path_list(path_list=['mrna/arghdia/2h/'])

    def test_get_image_list_wrong_format(self):
        with self.assertRaises(ValueError):
            retrieved_paths = self.repo.get_image_path_list(path_list=['mrna/arghdia/2h'])

    def test_get_image_list(self):
        # does not care if there is a "/" or not in the beginning
        retrieved_paths = self.repo.get_image_path_list(path_list=['/mrna/arhgdia/2h/', 'protein/arhgdia/3h/'])
        self.assertEqual(len(retrieved_paths), 10)
        self.assertIn("/mrna/arhgdia/2h/1", retrieved_paths) # have to add "/" in the beginning
        self.assertIn("/protein/arhgdia/3h/1", retrieved_paths)
        self.assertNotIn("/mrna/arhgdia/2h/1/spots", retrieved_paths)

    def test_get_image_list_final_level(self):
        with self.assertRaises(AttributeError):
            retrieved_paths = self.repo.get_image_path_list(path_list=['mrna/arhgdia/2h/1/spots/'])

    def test_spots_are_present(self):
        self.assertTrue(self.repo.is_present('mrna/arhgdia/2h/1/spots'), "should have found spots")
        self.assertFalse(self.repo.is_present('mrna/arhgdia/2h/1/spottttts'), "should have not found absent path")

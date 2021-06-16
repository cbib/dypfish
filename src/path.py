#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import os
import pathlib


def get_or_create_dir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname


global_root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
global_data_dir = pathlib.Path(global_root_dir, "data")
global_example_data = pathlib.Path(global_data_dir, 'example_hdf5')
example_config_path = pathlib.Path(global_example_data, "example_datasets_config.json")
test_config_path = pathlib.Path(global_root_dir, "src/tests/test_config.json")
basic_file_basename = 'basic.h5'

# -*- coding: utf-8 -*-

import os

# TODO
def get_or_create_dir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname

root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
#config_dir = os.path.join(root_dir, 'config')
log_dir = get_or_create_dir(os.path.join(root_dir, 'logs'))
src_dir = get_or_create_dir(os.path.join(root_dir, 'src'))
analysis_dir = get_or_create_dir(os.path.join(root_dir, 'analysis/'))
analysis_data_dir = get_or_create_dir(os.path.join(root_dir, 'data/'))
data_dir = get_or_create_dir(os.path.join(root_dir, 'data/'))


test_data_dir = get_or_create_dir(os.path.join(root_dir, 'test_data/'))
raw_data_dir = '/home/ben/dypfish/dypFISH_data/raw_image_data'
scratch_data_dir = '/home/ben/dypfish/Basics_Descriptors/preprocessed_data/new_scratch_data/scratch_if'  # TODO
basic_file_basename = 'basic'
sec_file_basename = 'secondary'
path_data = raw_data_dir
basic_file_path = analysis_data_dir + 'basic.h5'
secondary_file_path = analysis_data_dir + 'secondary.h5'
mtoc_file_path = analysis_data_dir + 'mtoc.h5'
h_star_file_path = analysis_data_dir + 'h_star.h5'

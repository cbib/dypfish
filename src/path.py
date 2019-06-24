# -*- coding: utf-8 -*-

import os

def get_or_create_dir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname

root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
log_dir = get_or_create_dir(os.path.join(root_dir, 'logs'))
src_dir = get_or_create_dir(os.path.join(root_dir, 'src'))
analysis_dir = get_or_create_dir(os.path.join(root_dir, 'analysis/'))
analysis_data_dir = get_or_create_dir(os.path.join(root_dir, 'data/'))
data_dir = get_or_create_dir(os.path.join(root_dir, 'data/'))
raw_data_dir = get_or_create_dir(os.path.join(root_dir, 'raw_image_data/'))
scratch_data_dir = get_or_create_dir(os.path.join(root_dir, 'Basics_Descriptors/preprocessed_data/new_scratch_data/scratch_if'))
test_data_dir = get_or_create_dir(os.path.join(root_dir, 'test_data/'))

basic_file_basename = 'basic'
sec_file_basename = 'secondary'
path_data = raw_data_dir
#basic_file_path = analysis_data_dir + 'genuine_basic.h5'
basic_file_path = analysis_data_dir + 'basic_scipy_0.19.0.h5'
secondary_file_path = analysis_data_dir + 'test_secondary.h5'
mtoc_file_path = analysis_data_dir + 'mtoc.h5'
h_star_file_path = analysis_data_dir + 'h_star.h5'

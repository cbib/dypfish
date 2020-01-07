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
config_data_dir = get_or_create_dir(os.path.join(root_dir, 'config/'))
data_dir = get_or_create_dir(os.path.join(root_dir, 'data/'))
raw_data_dir = get_or_create_dir(os.path.join(root_dir, 'raw_image_data/'))
test_data_dir = get_or_create_dir(os.path.join(root_dir, 'test_data/'))

basic_file_basename = 'basic'
sec_file_basename = 'secondary'
# TODO replace path_data everywhere by raw_data_dir
path_data = raw_data_dir


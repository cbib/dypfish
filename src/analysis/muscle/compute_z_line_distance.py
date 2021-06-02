#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import numpy as np
import pandas as pd
from loguru import logger

import constants
from image_set import ImageSet
from path import global_root_dir
from plot import profile
from repository import H5RepositoryWithCheckpoint


def compute_zline_distance(repo, molecule_list, timepoints, z_line_spacing):
    all_median_profiles = []
    for molecule in molecule_list:
        for timepoint in timepoints:
            image_set = ImageSet(repo, [f"{'mrna'}/{molecule}/{timepoint}/"])
            if image_set.__sizeof__() < 5:
                logger.warning("Image set is small for {}", molecule)
            total_profile = image_set.compute_zline_distance(z_line_spacing)
            all_median_profiles.append(np.median(total_profile, axis=0))
    return all_median_profiles


constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/muscle/config_muscle.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
z_line_spacing = constants.analysis_config['Z_LINE_SPACING']

molecule_list = ['actn2', 'gapdh']
timepoints = ['mature']
all_median_profiles_mature = compute_zline_distance(analysis_repo, molecule_list, timepoints, z_line_spacing)
df_mature = pd.DataFrame(all_median_profiles_mature)

molecule_list = ['actn2']
timepoints = ['immature']
all_median_profiles_immature = compute_zline_distance(analysis_repo, molecule_list, timepoints, z_line_spacing)
df_immature = pd.DataFrame(all_median_profiles_immature)


all_median_profiles = []
figure_title = 'z line spots distance profile'
genes = ["actn2 mature", "gapdh mature", "actn2 immature"]
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_GRAPH']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)

all_median_profiles.append(df_mature.loc[0].values)
all_median_profiles.append(df_mature.loc[1].values)
all_median_profiles.append(df_immature.loc[0].values)
profile(all_median_profiles, tgt_fp, keep=0)

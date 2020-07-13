#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
from plot import profile
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint
import numpy as np
import pandas as pd


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

""" 
Figure 7B - mRNA distance profiles. 
For each mRNA we computed its distance to the closest Z-lines, 
which allowed us to count the number of mRNAs having a certain distance to Z-lines. 
Normalized median counts are represented on the y axis.  
A higher  number of actn2 immature mRNA falls inside or close to Z-lines compared to mature fibers, 
suggesting greater clustering of mRNA between Z-lines for mature actn2. 
The Z-line distance is a Euclidean distance.
"""
logger.info("muscle mRNA distance profile for the mRNA muscle data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/muscle/config_muscle.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

logger.info("compute muscle mRNA distance profile for the mRNA Actn2 and Gapdh muscle mature cells")
mature_cells=compute_zline_distance(analysis_repo,['actn2','gapdh'], ['mature'], constants.analysis_config['Z_LINE_SPACING'])

logger.info("compute muscle mRNA distance profile for the mRNA Actn2 muscle immature cells")
immature_cells = compute_zline_distance(analysis_repo, ['actn2'],  ['immature'], constants.analysis_config['Z_LINE_SPACING'])
all_median_profiles = np.concatenate((mature_cells,immature_cells))

logger.info("generate plot for muscle mRNA distance profile for the mRNA Actn2 and Gapdh muscle mature cells")
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_GRAPH']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
profile(all_median_profiles, constants.analysis_config['MRNA_GENES_LABEL'], constants.analysis_config['Z_LINE_SPACING'], tgt_fp)
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import constants
from random import *
import numpy as np
from plot import plot_figure
from image_set import ImageSet
from path import global_root_dir
from helpers import mean_absolute_deviation, median_absolute_deviation
from repository import H5RepositoryWithCheckpoint


def get_spots_peripheral_distance_2D(spots, nucleus_mask, cell_mask, periph_dist_map):
    spots_peripheral_distance_2D = []
    for spot in spots:
        if nucleus_mask[spot[1], spot[0]] == 0:
            if cell_mask[spot[1], spot[0]] == 1 and periph_dist_map[spot[1], spot[0]] == 0:
                spots_peripheral_distance_2D.append(int(1))
            else:
                spots_peripheral_distance_2D.append(periph_dist_map[spot[1], spot[0]])
    return np.array(spots_peripheral_distance_2D)

def compute_stability(gene, force2D =True):
    total_mads =[]
    imageset = ImageSet(analysis_repo, ["mrna/"+gene+"/3h"], force2D=force2D)
    spots_peripheral_distances = imageset.compute_spots_peripheral_distance()
    # spots_peripheral_distances = np.array([get_spots_peripheral_distance_2D(image.get_spots(), image.get_nucleus_mask(),image.get_cell_mask(),image.get_cell_mask_distance_map()) for image inimageset.get_images()])
    peripheral_profiles = np.zeros((len(spots_peripheral_distances), 10))
    for j in range(len(spots_peripheral_distances)):
        for i in range(0, 10):
            peripheral_profiles[j, i] = float(len(np.where(
                (spots_peripheral_distances[j] >= ((i * 10) + 1)) & (spots_peripheral_distances[j] <= (i + 1) * 10))[
                                                      0]) / float(len(spots_peripheral_distances[j])))
    for j in range(500):
        mads = []
        for i in range(1, len(imageset.images)-1):
            arr = peripheral_profiles[np.random.choice(peripheral_profiles.shape[0], i, replace=True)]
            rand_idx = randint(0, peripheral_profiles.shape[0] - 1)
            mean_arr = np.mean(arr, axis=0)
            arr_diff = mean_arr - peripheral_profiles[rand_idx, :]
            mse = mean_absolute_deviation(arr_diff)
            mads.append(mse)
        total_mads.append(mads)
    return total_mads


"""
Figure 1.E 
Mean Absolute Deviation of Arhgdia mRNA distribution for peripheral fraction descriptors of a randomly selected cell 
from a pooled average of up to ~40 cells for cultured and micropatterned cells.
"""

constants.init_config(
    analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/stability/config_original.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

total_mads = []
for gene in constants.analysis_config['MRNA_GENES']:
    total_mads.append(compute_stability(gene, force2D =True))

tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
plot_figure(total_mads[0], total_mads[1], tgt_fp)


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import constants
import helpers
from loguru import logger
from plot import spline_graph, heatmap
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def plot_spline(grid_mat, mask_count, im):
    # spline graph density by vertical quadrants
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_GRAPH_STRIPE'].format(
        image=str(constants.analysis_config['STRIPE_NUM']) + im._path.replace("/", "_") + "_" + str(mask_count))
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    spline_graph(grid_mat, tgt_fp, constants.analysis_config['STRIPE_NUM'])


def plot_heatmap(grid_mat, mask_count, im):
    # heatmap density by vertical quadrants
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_HEATMAP'].format(
        image=str(constants.analysis_config['STRIPE_NUM']) + im._path.replace("/", "_") + "_" + str(mask_count))
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    heatmap(grid_mat, tgt_fp, constants.analysis_config['STRIPE_NUM'])


def compute_hist_mode(dist):
    # compute histogram mod
    hx, hy, _ = plt.hist(dist)
    bin_max = np.where(hx == hx.max())[0]
    hist_mod = hy[bin_max]
    if len(hist_mod) > 1:
        hist_mod = np.mean(hist_mod)
    return hist_mod


constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/muscle/config_muscle.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

"""
Figure 7.C 
The mRNA local density was computed between two nuclei. 
Each cell was quantized in vertical quadrants (constant STRIPE_NUM) and relative concentration of mRNA in each quadrant was computed 
by normalizing the counts by the relevant surface.  
"""

logger.info("Compute mRNA local density for the mRNA muscle data")
for g in constants.analysis_config['MRNA_GENES_LABEL']:
    [gene, timepoint] = g.split(" ")
    image_set = ImageSet(analysis_repo, [f"{'mrna'}/{gene}/{timepoint}/"])
    nuc_dist, nucs_dist, cell_masks, nucs_pos = image_set.compute_cell_mask_between_nucleus_centroid()
    # compute histogram mod
    hist_mod = compute_hist_mode(nuc_dist)
    # compute density by stripe and build spline wave graph
    image_counter = 0
    for i, im in enumerate(image_set.get_images()):
        nucleus_centroids = im.get_multiple_nucleus_centroid()
        mask_count = 0
        for j in range(len(cell_masks[i])):
            if hist_mod - 250 < nucs_dist[i][j] < hist_mod + 250:
                mask_c = cell_masks[i][j].copy()
                mask_reduced = mask_c[:, nucs_pos[i][j][0] - 50:nucs_pos[i][j][1] + 50]
                spots_reduced = im.keep_cell_mask_spots(cell_masks[i][j])
                spots_reduced = np.array(spots_reduced).reshape((len(spots_reduced), 3))
                spots_reduced[:, 0] -= nucs_pos[i][j][0] - 50
                grid_mat = im.build_density_by_stripe(spots_reduced, mask_reduced,stripe_num=constants.analysis_config['STRIPE_NUM'])
                plot_spline(grid_mat, mask_count, im)
                plot_heatmap(grid_mat, mask_count, im)

                mask_count += 1
        image_counter += 1

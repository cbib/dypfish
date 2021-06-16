#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import matplotlib.pyplot as plt
import numpy as np

import constants
import helpers
import plot
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/muscle/config_muscle.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
z_line_spacing = constants.analysis_config['Z_LINE_SPACING']
band_n = constants.analysis_config['STRIPE_NUM']

molecule_type = ['/mrna']
genes = ['actn2-mature', 'gapdh-mature', 'actn2-immature']

# figure 8C degree of clustering for muscle data

for g in genes:
    [gene, timepoint] = g.split("-")
    image_set = ImageSet(analysis_repo, [f"{'mrna'}/{gene}/{timepoint}/"])
    nuc_dist, nucs_dist, cell_masks, nucs_pos = image_set.compute_cell_mask_between_nucleus_centroids()

    # compute histogram mod
    hx, hy, _ = plt.hist(nuc_dist)
    bin_max = np.where(hx == hx.max())[0]
    hist_mod = hy[bin_max]
    if len(hist_mod) > 1:
        hist_mod = np.mean(hist_mod)

    # compute density by stripe and build spline wave graph
    image_counter = 0
    for im in image_set.get_images():
        nucleus_mask = im.get_nucleus_mask()
        nucleus_centroids = im.get_multiple_nucleus_centroid()
        spots = im.get_spots()
        z_lines = im.get_z_lines_masks()
        cell_masks_im = cell_masks[image_counter]
        nuc_dist_im = nucs_dist[image_counter]
        nuc_pos_im = nucs_pos[image_counter]
        mask_count = 0
        for i in range(len(cell_masks_im)):
            mask = cell_masks_im[i]
            nuc_d = nuc_dist_im[i]
            nuc_pos = nuc_pos_im[i]
            if hist_mod - 250 < nuc_d < hist_mod + 250:
                spots_reduced = helpers.keep_cell_mask_spots(spots, mask)
                spots_reduced = np.array(spots_reduced).reshape((len(spots_reduced), 3))
                mask_c = mask.copy()
                mask_reduced = mask_c[:, nuc_pos[0] - 50:nuc_pos[1] + 50]
                spots_reduced[:, 0] -= nuc_pos[0] - 50
                grid_mat = helpers.build_density_by_stripe(spots_reduced, z_lines, mask_reduced, band_n=band_n)

                # spline graph density by band_n
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_GRAPH_STRIPE'].format(
                    image=str(band_n) + im._path.replace("/", "_") + "_" + str(mask_count))
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                      tgt_image_name)
                plot.spline_graph(grid_mat, tgt_fp, band_n)

                # heatmap density by band_n
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_HEATMAP'].format(
                    image=str(band_n) + im._path.replace("/", "_") + "_" + str(mask_count))
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                      tgt_image_name)
                plot.heatmap(grid_mat, tgt_fp, band_n)
                mask_count += 1
        image_counter += 1

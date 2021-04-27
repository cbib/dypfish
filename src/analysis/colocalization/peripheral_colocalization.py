#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import pprint as pp
import pandas as pd
from loguru import logger
import constants
import plot
from plot import compute_heatmap
import helpers
from helpers import calculate_colocalization_score, open_repo
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir

def compute_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, molecule_type, quadrants_num=4):
    cs_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    if molecule_type == 'mrna':
        timepoints = constants.dataset_config['TIMEPOINTS_MRNA']
    else:
        timepoints = constants.dataset_config['TIMEPOINTS_PROTEIN']
    for gene in constants.analysis_config['PROTEINS']:
        mean_densities = []
        for timepoint in timepoints:
            image_set = ImageSet(analysis_repo, [molecule_type + "/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_peripheral_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num,
                                                                                       stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            mean_densities.append(mrna_tp_df.mean(axis=0).values)
        cs_dict[gene] = mean_densities
    return cs_dict

# Figure 5C Analysis peripheral Colocalization Score (CS) for original data (5 figures)
if __name__ == '__main__':
    logger.info("Colocalization Score")
    conf_full_path = pathlib.Path(global_root_dir, "src/analysis/colocalization/config_original_periph.json")
    constants.init_config(analysis_config_js_path=conf_full_path)
    repo = open_repo()
    mrna_cs_dict = compute_peripheral_relative_density_per_quadrants_and_slices(repo, 'mrna', quadrants_num=8)
    prot_cs_dict = compute_peripheral_relative_density_per_quadrants_and_slices(repo, 'protein', quadrants_num=8)

    css, p_vals = [], []
    for gene in constants.analysis_config['PROTEINS']:
        cs, p, ranking = calculate_colocalization_score(mrna_cs_dict[gene], prot_cs_dict[gene],
                                                            constants.dataset_config['TIMEPOINTS_NUM_MRNA'],
                                                            constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
        css.append(cs)
        p_vals.append(p)
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS'].format(gene=gene)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        compute_heatmap(ranking, gene, tgt_fp)

    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS_HISTOGRAM']
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
    plot.bar_profile(css, tgt_fp, constants.analysis_config['PLOT_COLORS'])




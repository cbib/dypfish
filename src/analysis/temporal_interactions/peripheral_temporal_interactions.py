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
from helpers import calculate_temporal_interaction_score
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir
import itertools
import matplotlib.pyplot as plt


def compute_mrna_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=4):
    mrna_tis_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['MRNA_GENES']:
        mrna_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_peripheral_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num,
                                                                                       stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            mrna_median.append(mrna_tp_df.mean(axis=0).values)
        mrna_tis_dict[gene] = mrna_median

    return mrna_tis_dict


def compute_protein_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=4):
    prot_tis_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['PROTEINS']:
        prot_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(analysis_repo, ["protein/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_peripheral_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num,
                                                                                       stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            prot_median.append(mrna_tp_df.mean(axis=0).values)
        prot_tis_dict[gene] = prot_median

    return prot_tis_dict


# Figure 5C
logger.info("Temporal interaction score for the mRNA original data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,
                                                           "src/analysis/temporal_interactions/config_original_periph.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_tis_dict = compute_mrna_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
prot_tis_dict = compute_protein_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)

tiss = []
p_vals = []
for gene in constants.analysis_config['PROTEINS']:
    mrna_list = mrna_tis_dict[gene]
    prot_list = prot_tis_dict[gene]
    (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list,
                                                             constants.dataset_config['TIMEPOINTS_NUM_MRNA'],
                                                             constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
    tiss.append(tis)
    p_vals.append(p)
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS'].format(gene=gene)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    compute_heatmap(ranking, gene, tgt_fp)

tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS_HISTOGRAM']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
plot.bar_profile(tiss, constants.analysis_config['PROTEINS'], tgt_fp, compute_median_and_error=False)

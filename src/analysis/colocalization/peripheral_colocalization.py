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
import itertools
import matplotlib.pyplot as plt


def compute_mrna_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=4):
    _mrna_cs_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for g in constants.analysis_config['MRNA_GENES']:
        mrna_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(g, timepoint)])
            arr = image_set.compute_peripheral_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num,
                                                                                       stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            mrna_median.append(mrna_tp_df.mean(axis=0).values)
        _mrna_cs_dict[g] = mrna_median

    return _mrna_cs_dict


def compute_protein_peripheral_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=4):
    protein_cs_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for g in constants.analysis_config['PROTEINS']:
        prot_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(analysis_repo, ["protein/{0}/{1}/".format(g, timepoint)])
            arr = image_set.compute_peripheral_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num,
                                                                                       stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            prot_median.append(mrna_tp_df.mean(axis=0).values)
        protein_cs_dict[g] = prot_median

    return protein_cs_dict


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/colocalization/config_original_periph.json", []]
]

# Figure 5C Analysis peripheral Colocalization Score (CS) for original data (5 figures)
if __name__ == '__main__':
    logger.info("Colocalization Score")
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        mrna_cs_dict = compute_mrna_peripheral_relative_density_per_quadrants_and_slices(repo, quadrants_num=8)
        prot_cs_dict = compute_protein_peripheral_relative_density_per_quadrants_and_slices(repo, quadrants_num=8)

        css = []
        p_vals = []
        for gene in constants.analysis_config['PROTEINS']:
            mrna_list = mrna_cs_dict[gene]
            prot_list = prot_cs_dict[gene]
            (cs, p, ranking) = calculate_colocalization_score(mrna_list, prot_list,
                                                                     constants.dataset_config['TIMEPOINTS_NUM_MRNA'],
                                                                     constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
            css.append(cs)
            p_vals.append(p)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS'].format(gene=gene)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            compute_heatmap(ranking, gene, tgt_fp)

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS_HISTOGRAM']
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.bar_profile(css, tgt_fp, constants.analysis_config['PLOT_COLORS'])




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
from helpers import open_repo
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir
import itertools
import matplotlib.pyplot as plt


def compute_protein_relative_density_per_quadrants_and_slices(_analysis_repo, _quadrants_num=4):
    _prot_tis_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['PROTEINS']:
        prot_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(_analysis_repo, ["protein/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=_quadrants_num,
                                                                            stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            prot_median.append(mrna_tp_df.mean(axis=0).values)
        _prot_tis_dict[gene] = prot_median

    return _prot_tis_dict


def compute_mrna_relative_density_per_quadrants_and_slices(_analysis_repo, _quadrants_num=4):
    _mrna_tis_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['MRNA_GENES']:
        mrna_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(_analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=_quadrants_num,
                                                                            stripes=stripes)
            mrna_tp_df = pd.DataFrame(arr)
            mrna_median.append(mrna_tp_df.mean(axis=0).values)
        _mrna_tis_dict[gene] = mrna_median

    return _mrna_tis_dict


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/temporal_interactions/config_original.json", []],
    ["src/analysis/temporal_interactions/config_nocodazole_arhgdia.json", ["arhgdia", "Nocodazole+"]],
    ["src/analysis/temporal_interactions/config_nocodazole_pard3.json", ["pard3", "Nocodazole+"]]
]


# Figure 5D Analysis TIS for original data (5 figures)
# Figure 6E Analysis TIS for nocodazole arhgdia data (3 figures)
# Figure 6E Analysis TIS for nocodazole pard3 data (3 figures)
if __name__ == '__main__':

    for conf in configurations:
        logger.info("Temporal interaction score")
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        # Use annot=True if you want to add stats annotation in plots
        mrna_tis_dict = compute_mrna_relative_density_per_quadrants_and_slices(repo, _quadrants_num=8)
        prot_tis_dict = compute_protein_relative_density_per_quadrants_and_slices(repo, _quadrants_num=8)

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
            if len(conf[1])==0:
                compute_heatmap(ranking, gene, tgt_fp)
            else:
                compute_heatmap(ranking, gene, tgt_fp, size=2, xtickslabel=['3h', '5h'], ytickslabel=['3h', '5h'])
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS_HISTOGRAM']
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.bar_profile_simple(tiss, tgt_fp)


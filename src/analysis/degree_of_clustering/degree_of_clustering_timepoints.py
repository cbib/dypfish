#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import math
from helpers import open_repo

from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


"""
for each timepoint, process_data function compute median, low and up envelope
normalize data by mean of median timepoint result
"""


def process_data(gene, molecule_type, repo, timepoints):
    data = np.zeros((3, len(timepoints)))
    for i in range(len(timepoints)):
        image_set = ImageSet(repo, ["{0}/{1}/{2}/".format(molecule_type, gene, timepoints[i])])
        # degree_of_clustering = np.array(image_set.compute_mtoc_dependent_degree_of_clustering())
        degree_of_clustering = np.array(image_set.compute_degree_of_clustering())
        results_median = np.median(degree_of_clustering)
        err = np.median(
            np.abs(np.tile(np.median(degree_of_clustering), (1, len(degree_of_clustering))) - degree_of_clustering))
        upp_env = results_median + err
        low_env = results_median - err
        data[0, i] = results_median
        data[1, i] = upp_env
        data[2, i] = low_env
    normalized_data = data / np.mean(data[0, :])

    return normalized_data


configurations = [
    ["src/analysis/degree_of_clustering/config_original.json", "", "", ""]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()

        # Figure 2F  : Dynamic profile of degree of clustering for original data
        plot_colors = constants.analysis_config['PLOT_COLORS']
        for i, gene in enumerate(constants.analysis_config['MRNA_GENES']):
            mrna_data = process_data(gene, "mrna", repo, constants.dataset_config['TIMEPOINTS_MRNA'])
            if gene in constants.analysis_config['PROTEINS']:
                protein_data = process_data(gene, "protein", repo, constants.dataset_config['TIMEPOINTS_PROTEIN'])
            # generate image
            tgt_image_name = constants.analysis_config['DYNAMIC_FIGURE_NAME_FORMAT'].format(gene=gene)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.dynamic_profiles(mrna_data, protein_data, gene, 'Time(hrs)', 'Degree of clustering(*)', tgt_fp,
                                  plot_colors[i])

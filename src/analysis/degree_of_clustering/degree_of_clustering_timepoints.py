#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import math

import pandas as pd
from helpers import open_repo, color_variant

from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def plot_dynamic_barplot(analysis_repo):
    plot_colors = constants.analysis_config['PLOT_COLORS']

    for i, gene in enumerate(constants.analysis_config['PROTEINS']):
        df = pd.DataFrame(columns=["Gene", "Molecule_type", "Timepoint", "d_of_c"])
        dict = {'Gene': [], 'Molecule_type': [], 'Timepoint': [], 'DoC': []}
        cpt = 0
        for tp in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(analysis_repo, ["{0}/{1}/{2}/".format("mrna", gene, tp)])
            degree_of_clustering = np.array(image_set.compute_degree_of_clustering())
            for d_of_c in degree_of_clustering:
                dict["Gene"].append(gene)
                dict["Molecule_type"].append("mrna")
                dict["Timepoint"].append(tp)
                dict["d_of_c"].append(d_of_c)
            cpt += 1
        dict["d_of_c"] = dict["d_of_c"] / np.mean(dict["d_of_c"])
        df = pd.concat([df, pd.DataFrame(dict)])

        dict = {'Gene': [], 'Molecule_type': [], 'Timepoint': [], 'DoC': []}
        cpt = 0
        for tp in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(analysis_repo, ["{0}/{1}/{2}/".format("protein", gene, tp)])
            # degree_of_clustering = np.array(image_set.compute_mtoc_dependent_degree_of_clustering())
            degree_of_clustering = np.array(image_set.compute_degree_of_clustering())
            for d_of_c in degree_of_clustering:
                dict["Gene"].append(gene)
                dict["Molecule_type"].append("protein")
                dict["Timepoint"].append(tp)
                dict["d_of_c"].append(d_of_c)
            cpt += 1
        dict["d_of_c"] = dict["d_of_c"] / np.mean(dict["d_of_c"])
        df = pd.concat([df, pd.DataFrame(dict)])
        tgt_image_name = constants.analysis_config['DYNAMIC_FIGURE_NAME_FORMAT'].format(gene=gene)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)

        my_pal = {"mrna": str(plot_colors[i]), "protein": str(color_variant(plot_colors[i], +80))}
        plot.sns_barplot_simple(df, my_pal, tgt_fp, x="Timepoint", y="d_of_c", hue="Molecule_type")

configurations = [
    ["src/analysis/degree_of_clustering/config_original.json", "", "", ""]
]

# # Figure 2F  : Dynamic profile of degree of clustering for original data
if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        plot_dynamic_barplot(repo)
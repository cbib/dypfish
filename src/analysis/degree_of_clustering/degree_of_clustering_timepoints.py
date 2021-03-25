#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import helpers

import pandas as pd
from helpers import open_repo, color_variant

from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def plot_dynamic_barplot(analysis_repo):
    plot_colors = constants.analysis_config['PLOT_COLORS']

    # paired mRNA-protein barplots, so we go through proteins (we have less proteins than mRNA)
    for i, gene in enumerate(constants.analysis_config['PROTEINS']):
        df = pd.DataFrame(columns=["Gene", "Molecule", "Timepoint", "DoC"])
        for molecule, timepoints in zip(["mrna", "protein"], ['TIMEPOINTS_MRNA', 'TIMEPOINTS_PROTEIN']):
            dict = {"Gene": [], "Molecule": [], "Timepoint": [], "DoC": []}
            for tp in constants.dataset_config[timepoints]:
                image_set = ImageSet(analysis_repo, ["{0}/{1}/{2}/".format(molecule, gene, tp)])
                degree_of_clustering = np.array(image_set.compute_degree_of_clustering())
                for d_of_c in degree_of_clustering:
                    dict["Gene"].append(gene)
                    dict["Molecule"].append(molecule)
                    dict["Timepoint"].append(tp)
                    dict["DoC"].append(d_of_c)
            dict["DoC"] = dict["DoC"] / np.mean(dict["DoC"])
            df = pd.concat([df, pd.DataFrame(dict)])

        my_pal = {"mrna": str(plot_colors[i]), "protein": str(color_variant(plot_colors[i], +80))}
        plot.sns_barplot_simple(df, palette=my_pal, gene=gene, x="Timepoint", y="DoC", hue="Molecule")

configurations = [
    ["src/analysis/degree_of_clustering/config_original.json", "", "", ""]
]

# # Figure 2F  : Dynamic profile of degree of clustering for original data, both mRNA and proteins
if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        plot_dynamic_barplot(repo)
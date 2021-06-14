#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
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


def mean_column(group):
    group["d_of_c"] = group["d_of_c"].apply(np.median)
    return group


data = [{'beta_actin': [0.5, 1, 0.8, 1, 1], 'arhgdia': [1, 1, 1.5, 1, 1],
         'gapdh': [1, 0.9, 0.7, 0.7, 1], 'pard3': [0.8, 1, 0.8, 1.2, 1]},
        {'beta_actin': [0.6, 1, 1, 1, 0.6], 'arhgdia': [0.7, 0.7, 1, 0.7, 0.7],
         'gapdh': [1, 0.5, 1, 0.6, 0.8], 'pard3': [0.8, 0.8, 1, 0.6, 0.7]}]
factor = pd.DataFrame(data, index=['mrna', 'protein'])

def plot_dynamic_barplot(analysis_repo):
    '''
    Formats the data and calls the plotting function
    '''
    plot_colors = constants.analysis_config['PLOT_COLORS']

    # paired mRNA-protein barplots, so we go through proteins (we have less proteins than mRNA)
    tp_mrna = constants.dataset_config['TIMEPOINTS_MRNA']
    tp_proteins = constants.dataset_config['TIMEPOINTS_PROTEIN']
    all_timepoints = np.sort(list(set(tp_mrna) | set(tp_proteins)))
    for i, gene in enumerate(constants.analysis_config['PROTEINS']):
        df = pd.DataFrame(columns=["Molecule", "Timepoint", "d_of_c", "error", "CI"])
        for molecule, timepoints in zip(["mrna", "protein"],[tp_mrna, tp_proteins]):
            for j, tp in enumerate(all_timepoints):
                if tp not in timepoints:
                    df = df.append({"Molecule": molecule, "Timepoint": tp, "error": 0, "CI": [0, 0],
                                    "d_of_c": 0}, ignore_index=True)
                    continue
                image_set = ImageSet(analysis_repo, ["{0}/{1}/{2}/".format(molecule, gene, tp)])
                degree_of_clustering = np.log(image_set.compute_degree_of_clustering()) #* factor[gene][molecule][j]
                err = helpers.sem(degree_of_clustering, factor=6)
                lower, higher = helpers.median_confidence_interval(degree_of_clustering)
                df = df.append({"Molecule": molecule, "Timepoint": tp, "error": err, "CI": [lower, higher],
                                "d_of_c": degree_of_clustering}, ignore_index=True)
        df = df.sort_values('Timepoint')
        df = df.groupby('Molecule').apply(mean_column)
        my_pal = {"mrna": str(plot_colors[i]), "protein": str(color_variant(plot_colors[i], +80))}
        tgt_image_name = constants.analysis_config['DYNAMIC_FIGURE_NAME_FORMAT'].format(gene=gene)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        plot.bar_profile_median_timepoints(df, palette=my_pal, figname=tgt_fp, fixed_yscale=15)


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
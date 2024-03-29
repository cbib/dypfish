#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import collections
import constants
import plot
import numpy as np
import helpers
from helpers import open_repo
from image_set import ImageSet
from loguru import logger
# this should be called as soon as possible
from path import global_root_dir


def compute_degree_of_clustering(genes_list, analysis_repo, molecule_type):
    gene2_degree_of_clustering = {}
    gene2median_degree_of_clustering = {}
    gene2error_degree_of_clustering = {}
    gene2confidence_interval = {}
    degrees_of_clustering = []

    for gene in genes_list:
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)])
        d_of_c = np.array(image_set.compute_degree_of_clustering())
        degrees_of_clustering.append(d_of_c)

    for gene, degree_of_clustering in zip(genes_list, degrees_of_clustering):
        degree_of_clustering = np.log(degree_of_clustering)
        gene2_degree_of_clustering[gene] = degree_of_clustering
        gene2median_degree_of_clustering[gene] = np.median(degree_of_clustering)
        # Standard error and CI computation
        gene2error_degree_of_clustering[gene] = helpers.sem(degree_of_clustering, factor=0)
        lower, higher = helpers.median_confidence_interval(degree_of_clustering)
        gene2confidence_interval[gene] = [lower, higher]

    return gene2_degree_of_clustering, gene2median_degree_of_clustering, gene2error_degree_of_clustering, gene2confidence_interval


def plot_bar_profile_median_and_violin(molecule_type, median_d_of_c, d_of_c,
                                       errors, confidence_interval, annot=False):
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    if molecule_type == 'mrna':
        xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    else:
        xlabels = constants.analysis_config['PROTEINS_LABEL']

    # generate the bar profile plot
    plot.bar_profile_median(median_d_of_c, errors.values(), molecule_type, xlabels,
                            tgt_fp, confidence_interval, annot=annot, data_to_annot=d_of_c)
    logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.violin_profile(d_of_c, tgt_fp, xlabels, rotation=0, annot=annot)


''' 
Figure 3D left panel: plots the log mRNA degree of clustering for original data
Figure 3D right panel: plots the log protein degree of clustering for original data
Figure 3F : plots the log protein degree of clustering for CHX
Figure 5B left panel: plots the log mRNA degree of clustering for prrc2c
Figure 5B right panel: plots the log protein degree of clustering for prrc2c
Figure S2E  log mRNA degree of clustering for CHX

'''

# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/degree_of_clustering/config_original.json", ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/degree_of_clustering/config_chx.json",
     ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/degree_of_clustering/config_prrc2c.json",
     ['arhgdia/control', "arhgdia/prrc2c_depleted"]]
]


def build_plots(analysis_repo, conf, annot=False):
    for molecule_type, molecules in zip(["mrna", "protein"], ['MRNA_GENES', 'PROTEINS']):
        if "original" in conf[0]:
            annot = False

        genes_list = constants.dataset_config[molecules]
        d_of_c, median_d_of_c, err, confidence_interval = compute_degree_of_clustering(genes_list, analysis_repo, molecule_type=molecule_type)
        # sort everything in the same way for plotting
        keyorder = conf[1]
        median_d_of_c = collections.OrderedDict(sorted(median_d_of_c.items(), key=lambda i: keyorder.index(i[0])))
        d_of_c = collections.OrderedDict(sorted(d_of_c.items(), key=lambda i: keyorder.index(i[0])))

        err = collections.OrderedDict(sorted(err.items(), key=lambda i: keyorder.index(i[0])))
        confidence_interval = collections.OrderedDict(sorted(confidence_interval.items(), key=lambda i: keyorder.index(i[0])))
        plot_bar_profile_median_and_violin(molecule_type, median_d_of_c, d_of_c, err, confidence_interval, annot=annot)


if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        # Use annot=True if you want to add stats annotation in your figures
        build_plots(repo, conf, annot=False)

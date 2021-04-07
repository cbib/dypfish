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
# this should be called as soon as possible
from path import global_root_dir

def compute_degree_of_clustering(genes_list, repo, molecule_type):
    gene2median_degree_of_clustering = {}
    gene2error_degree_of_clustering = {}
    gene2CI = {}
    degrees_of_clustering = []

    for gene in genes_list:
        image_set = ImageSet(repo, ['{0}/{1}/'.format(molecule_type, gene)])
        d_of_c = np.array(image_set.compute_degree_of_clustering())
        degrees_of_clustering.append(d_of_c)

    for gene, degree_of_clustering in zip(genes_list, degrees_of_clustering):
        degree_of_clustering = np.log(degree_of_clustering)
        gene2median_degree_of_clustering[gene] = np.median(degree_of_clustering)
        # Standard error and CI computation
        gene2error_degree_of_clustering[gene] = helpers.sem(degree_of_clustering, factor=0)
        lower, higher = helpers.median_confidence_interval(degree_of_clustering)
        gene2CI[gene] = [lower, higher]

    return gene2median_degree_of_clustering, gene2error_degree_of_clustering, gene2CI

''' 
Figure 2E left panel: plots the log mRNA degree of clustering normalized by log(0.5) for original
Figure 2E right panel: plots the log protein degree of clustering normalized by log(0.01) for original
Figure 2E left panel: plots the log mRNA degree of clustering normalized by log(0.5) for prrc2c
Figure 4B right panel: plots the log protein degree of clustering normalized by log(0.01) for prrc2c
Figure 2G : plots the log protein degree of clustering normalized by log(0.01) for CHX
Figure S2D right panel: plots the log protein degree of clustering normalized by log(0.01) for CHX
'''

# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/degree_of_clustering/config_original.json",
     ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/degree_of_clustering/config_chx.json",
     ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/degree_of_clustering/config_prrc2c.json",
     ['arhgdia/control', "arhgdia/prrc2c_depleted"]]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()

        for molecule_type, molecules in zip(["mrna", "protein"], ['MRNA_GENES', 'PROTEINS']):
            genes_list = constants.dataset_config[molecules]
            d_of_c, err, CI = compute_degree_of_clustering(genes_list, repo, molecule_type=molecule_type)
            # sort everything in the same way for plotting
            keyorder = conf[1]
            d_of_c = collections.OrderedDict(sorted(d_of_c.items(), key=lambda i: keyorder.index(i[0])))
            err = collections.OrderedDict(sorted(err.items(), key=lambda i: keyorder.index(i[0])))
            CI = collections.OrderedDict(sorted(CI.items(), key=lambda i: keyorder.index(i[0])))
            plot.bar_profile_median(d_of_c.keys(), d_of_c.values(), err.values(), CI.values(), molecule_type=molecule_type, fixed_yscale=15)


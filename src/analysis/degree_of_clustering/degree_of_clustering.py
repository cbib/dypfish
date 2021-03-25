#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import math
import helpers
import collections
import statsmodels.stats.api as sms
from helpers import open_repo
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir

def compute_degree_of_clustering(genes_list, repo, molecule_type):
    gene2median_degree_of_clustering = {}
    gene2error_degree_of_clustering = {}
    gene2CI = {}
    for gene in genes_list:
        image_set = ImageSet(repo, ['{0}/{1}/'.format(molecule_type, gene)])
        degree_of_clustering = np.array(image_set.compute_degree_of_clustering())
        median_degree_of_clustering = np.median(degree_of_clustering)
        gene2median_degree_of_clustering[gene] = math.log(median_degree_of_clustering)
        #Standard error and CI computation
        gene2error_degree_of_clustering[gene] = np.std(np.log(degree_of_clustering))/math.sqrt(len(degree_of_clustering))
        lower, higher = helpers.median_confidence_interval(degree_of_clustering)
        gene2CI[gene]= [gene2median_degree_of_clustering[gene] - math.log(lower),
                   math.log(higher) - gene2median_degree_of_clustering[gene]]

    return gene2median_degree_of_clustering, gene2error_degree_of_clustering, gene2CI


''' 
Figure 2E left panel: plots the log mRNA degree of clustering for original
Figure 2E right panel: plots the log protein degree of clustering for original
Figure 4B left panel: plots the log mRNA degree of clustering for prrc2c
Figure 4B right panel: plots the log protein degree of clustering for prrc2c
Figure 2G bottom panel: plots the log protein degree of clustering for CHX
Figure S2D right panel: plots the log mRNA degree of clustering  for CHX
'''

configurations = [
    ["src/analysis/degree_of_clustering/config_original.json",
     ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13'], ['beta_actin', 'arhgdia', 'gapdh', 'pard3']], #, "", "", ""],
    ["src/analysis/degree_of_clustering/config_chx.json",
     ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX'], ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']], #"arhgdia", "chx+", "Gene"],
    ["src/analysis/degree_of_clustering/config_prrc2c.json",
     ['arhgdia/control', "arhgdia/prrc2c_depleted"], ["arhgdia/control", "arhgdia/prrc2c_depleted"]]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()

        ## mRNA
        genes_list = constants.dataset_config['MRNA_GENES']
        gene2median_degree_of_clustering, gene2error_degree_of_clustering, gene2CI = compute_degree_of_clustering(genes_list, repo,
                                                                                                                  molecule_type="mrna")

        keyorder = conf[1]
        gene2median_degree_of_clustering = collections.OrderedDict(sorted(gene2median_degree_of_clustering.items(),
                                                                          key=lambda i: keyorder.index(i[0])))
        gene2error_degree_of_clustering = collections.OrderedDict(sorted(gene2error_degree_of_clustering.items(),
                                                                         key=lambda i: keyorder.index(i[0])))
        gene2CI = collections.OrderedDict(sorted(gene2CI.items(),
                                                 key=lambda i: keyorder.index(i[0])))
        plot.bar_profile_median(gene2median_degree_of_clustering.values(),
                                gene2median_degree_of_clustering.keys(),
                                gene2error_degree_of_clustering.values(), gene2CI, molecule_type="mrna")

        ## Proteins
        protein_time_points = constants.dataset_config['TIMEPOINTS_PROTEIN']
        protein_list = constants.dataset_config['PROTEINS']

        gene2median_degree_of_clustering, gene2error_degree_of_clustering, gene2CI = compute_degree_of_clustering(protein_list, repo,
                                                                                                                  molecule_type="protein")
        keyorder = conf[2]
        gene2median_degree_of_clustering = collections.OrderedDict(sorted(gene2median_degree_of_clustering.items(),
                                                                          key=lambda i: keyorder.index(i[0])))
        gene2error_degree_of_clustering = collections.OrderedDict(sorted(gene2error_degree_of_clustering.items(),
                                                                         key=lambda i: keyorder.index(i[0])))
        gene2CI = collections.OrderedDict(sorted(gene2CI.items(),
                                                 key=lambda i: keyorder.index(i[0])))
        plot.bar_profile_median(gene2median_degree_of_clustering.values(),
                                gene2median_degree_of_clustering.keys(),
                                gene2error_degree_of_clustering.values(), gene2CI, molecule_type="protein")

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
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def compute_degree_of_clustering(genes_list, repo, molecule_type):
    gene2median_degree_of_clustering = {}
    gene2error_degree_of_clustering = {}
    for gene in genes_list:
        image_set = ImageSet(repo, ['{0}/{1}/'.format(molecule_type, gene)])
        degree_of_clustering = np.array(image_set.compute_degree_of_clustering())

        median_degree_of_clustering = np.median(degree_of_clustering)
        gene2median_degree_of_clustering[gene] = math.log(median_degree_of_clustering)

        #Standard error computation
        gene2error_degree_of_clustering[gene] = np.std(np.log(degree_of_clustering))/math.sqrt(len(degree_of_clustering))

    return gene2median_degree_of_clustering, gene2error_degree_of_clustering


''' 
Figure 2E left panel: plots the log mRNA degree of clustering for original
Figure 2E right panel: plots the log protein degree of clustering normalized for original

Figure 2E left panel: plots the log mRNA degree of clustering normalized for prrc2c
Figure 4B right panel: plots the log protein degree of clustering normalized for prrc2c

Figure 2G : plots the log protein degree of clustering normalized for CHX
Figure S2D right panel: plots the log protein degree of clustering normalized for CHX
'''

configurations = [
    ["src/analysis/degree_of_clustering/config_original.json", "", "", ""],
    ["src/analysis/degree_of_clustering/config_chx.json", "arhgdia", "chx+", "Gene"],
    ["src/analysis/degree_of_clustering/config_prrc2c.json", "arhgdia", "prrc2c_depleted", "Timepoint"]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()

        ## mRNA
        genes_list = constants.dataset_config['MRNA_GENES']
        mrna_time_points = constants.dataset_config['TIMEPOINTS_MRNA']
        gene2median_degree_of_clustering, gene2error_degree_of_clustering = compute_degree_of_clustering(genes_list, repo,
                                                                                                         molecule_type="mrna")

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)

        plot.bar_profile_median(gene2median_degree_of_clustering.values(),
                                gene2median_degree_of_clustering.keys(),
                                gene2error_degree_of_clustering.values(), figname=tgt_fp)
        logger.info("Generated image at {}", tgt_fp)

        ## Proteins
        protein_time_points = constants.dataset_config['TIMEPOINTS_PROTEIN']
        protein_list = constants.dataset_config['PROTEINS']

        gene2median_degree_of_clustering, gene2error_degree_of_clustering = compute_degree_of_clustering(protein_list, repo,
                                                                                                         molecule_type="protein")
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        plot.bar_profile_median(gene2median_degree_of_clustering.values(),
                                gene2median_degree_of_clustering.keys(),
                                gene2error_degree_of_clustering.values(), figname=tgt_fp)
        logger.info("Generated image at {}", tgt_fp)

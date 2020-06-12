#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import math

from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


# Figure S2E : plots the log mRNA degree of clustering normalized by log(0.5)

logger.info("mRNA Degree of clustering for the original data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/degree_of_clustering/config_nocodazole.json"))

dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

base = math.log(0.5)
gene2median_mrna_degree_of_clustering = {}
gene2error_mrna_degree_of_clustering = {}
for gene in constants.analysis_config['MRNA_GENES']:
    image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
    # degree_of_clustering = np.array(image_set.compute_mtoc_dependent_degree_of_clustering())
    degree_of_clustering = np.array(image_set.compute_degree_of_clustering())

    median_degree_of_clustering = np.median(degree_of_clustering)
    gene2median_mrna_degree_of_clustering[gene] = math.log(median_degree_of_clustering) - base
    err = np.median(np.abs(median_degree_of_clustering - degree_of_clustering))
    gene2error_mrna_degree_of_clustering[gene] = math.log(median_degree_of_clustering + err) - math.log(
        median_degree_of_clustering) - base

# generate image
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
plot.bar_profile_median(gene2median_mrna_degree_of_clustering.values(), gene2median_mrna_degree_of_clustering.keys(),
                        gene2error_mrna_degree_of_clustering.values(), figname=tgt_fp)
logger.info("Generated image at {}", tgt_fp)


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import time

import numpy as np
from loguru import logger

import constants
import helpers
import plot
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def compute_relative_densities(analysis_repo, molecule_type, quadrants_num=4):
    densities = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    if molecule_type == 'mrna':
        timepoints = constants.dataset_config['TIMEPOINTS_MRNA']
        genes = constants.analysis_config['MRNA_GENES']
    else:
        timepoints = constants.dataset_config['TIMEPOINTS_PROTEIN']
        genes = constants.analysis_config['PROTEINS']

    for gene in genes:
        median_densities = []
        for timepoint in timepoints:
            image_set = ImageSet(analysis_repo, [molecule_type + "/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalised_quadrant_densities(quadrants_num=quadrants_num,
                                                                  peripheral_flag=False,
                                                                  stripes=stripes, stripes_flag=True)
            num_images = arr.shape[0] // (quadrants_num * stripes)
            aligned_densities = arr[:, 0].reshape(num_images, quadrants_num * stripes)
            median_densities_per_slice = np.nanmedian(aligned_densities, axis=0)
            median_densities.append(median_densities_per_slice)
        densities[gene] = median_densities

    return densities


if __name__ == '__main__':
    np.random.seed(int(round(time.time())))
    logger.info("Clustering analysis")
    conf_full_path = pathlib.Path(global_root_dir, "src/analysis/clustering/config_original.json")
    constants.init_config(analysis_config_js_path=conf_full_path)
    repo = helpers.open_repo()

    mrna_densities = compute_relative_densities(repo, 'mrna', quadrants_num=8)
    plot.plot_fine_grained_clusters('mrna', mrna_densities)
    prot_densities = compute_relative_densities(repo, 'protein', quadrants_num=8)
    plot.plot_fine_grained_clusters('protein', prot_densities)


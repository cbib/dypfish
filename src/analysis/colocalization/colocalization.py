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
from plot import plot_heatmap


def compute_relative_densities(analysis_repo, molecule_type, quadrants_num=4, peripheral_flag=False):
    densities = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    if peripheral_flag:
        stripes_flag = False
        stripes = 1
    else:
        stripes_flag = True

    if molecule_type == 'mrna':
        timepoints = constants.dataset_config['TIMEPOINTS_MRNA']
    else:
        timepoints = constants.dataset_config['TIMEPOINTS_PROTEIN']

    for gene in constants.analysis_config['PROTEINS']:
        median_densities, gene_clusters = [], []
        for timepoint in timepoints:
            image_set = ImageSet(analysis_repo, [molecule_type + "/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalised_quadrant_densities(quadrants_num=quadrants_num,
                                                                  peripheral_flag=peripheral_flag,
                                                                  stripes=stripes, stripes_flag=stripes_flag)
            num_images = arr.shape[0] // (quadrants_num * stripes)
            aligned_densities = arr[:, 0].reshape(num_images, quadrants_num * stripes)
            median_densities_per_slice = np.nanmedian(aligned_densities, axis=0)
            median_densities.append(median_densities_per_slice)
        densities[gene] = median_densities

    return densities


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/colocalization/config_original.json"],
    ["src/analysis/colocalization/config_original_periph.json"],
    ["src/analysis/colocalization/config_nocodazole_arhgdia.json"],
    ["src/analysis/colocalization/config_nocodazole_pard3.json"]
]

# Figure 5D Analysis Colocalization Score (CS) for original data (5 figures)
# Figure 6E Analysis Colocalization Score (CS) for nocodazole arhgdia data (3 figures)
# Figure 6E Analysis Colocalization Score (CS) for nocodazole pard3 data (3 figures)

if __name__ == '__main__':
    np.random.seed(int(round(time.time())))
    for conf in configurations:
        logger.info("Colocalization Score")
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = helpers.open_repo()
        peripheral_flag = "periph" in conf[0]
        stripes = constants.analysis_config["STRIPE_NUM"]
        mrna_densities = compute_relative_densities(repo, 'mrna', quadrants_num=8, peripheral_flag=peripheral_flag)
        prot_densities = compute_relative_densities(repo, 'protein', quadrants_num=8, peripheral_flag=peripheral_flag)

        css, results = [], {}
        for gene in constants.analysis_config['PROTEINS']:
            cs, p, ranking = helpers.calculate_colocalization_score(mrna_densities[gene], prot_densities[gene],
                                                                    constants.dataset_config['TIMEPOINTS_NUM_MRNA'],
                                                                    constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'],
                                                                    peripheral_flag, stripes=stripes, quadrants=8)
            css.append(cs)
            results[gene] = [cs, p]
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS'].format(gene=gene)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            if "original" in conf[0]:
                plot_heatmap(ranking, gene, tgt_fp)
                logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])
            else:
                plot_heatmap(ranking, gene, tgt_fp, size=2, xtickslabel=['3h', '5h'], ytickslabel=['3h', '5h'])
                logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS_HISTOGRAM']
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.bar_profile(css, tgt_fp, constants.analysis_config['PLOT_COLORS'])
        logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])
        logger.info("Colocalization scores and p-values are: {}", results)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import time
import constants
import plot
from plot import compute_heatmap
from helpers import calculate_colocalization_score, open_repo
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir

def compute_relative_density_per_quadrants_and_slices(analysis_repo, molecule_type, quadrants_num=4):
    protein_cs_dict = {}
    stripes = constants.analysis_config['STRIPE_NUM']
    if molecule_type == 'mrna':
        timepoints = constants.dataset_config['TIMEPOINTS_MRNA']
    else:
        timepoints = constants.dataset_config['TIMEPOINTS_PROTEIN']

    for gene in constants.analysis_config['PROTEINS']:
        mean_densities = []
        for timepoint in timepoints:
            image_set = ImageSet(analysis_repo, [molecule_type + "/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalised_quadrant_densities(quadrants_num=quadrants_num,
                                                                  peripheral_flag=False,
                                                                  stripes=stripes, stripes_flag=True)
            aligned_densities = arr[:, 0].reshape(image_set.__sizeof__(), quadrants_num * stripes)
            mean_densities_per_slice = np.nanmean(aligned_densities, axis=0)
            mean_densities.append(mean_densities_per_slice)
        protein_cs_dict[gene] = mean_densities

    return protein_cs_dict


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/colocalization/config_original.json", [], 10000],
    ["src/analysis/colocalization/config_nocodazole_arhgdia.json", ["arhgdia", "Nocodazole+"], 24],
    ["src/analysis/colocalization/config_nocodazole_pard3.json", ["pard3", "Nocodazole+"], 24]
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
        repo = open_repo()

        # Use annot=True if you want to add stats annotation in plots
        mrna_cs_dict = compute_relative_density_per_quadrants_and_slices(repo, 'mrna', quadrants_num=8)
        prot_cs_dict = compute_relative_density_per_quadrants_and_slices(repo, 'protein', quadrants_num=8)

        css, p_vals = [], {}
        for gene in constants.analysis_config['PROTEINS']:
            cs, p, ranking = calculate_colocalization_score(mrna_cs_dict[gene], prot_cs_dict[gene],
                                                              constants.dataset_config['TIMEPOINTS_NUM_MRNA'],
                                                              constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'],
                                                              permutation_num=conf[2])
            css.append(cs)
            p_vals[gene] = p
            print("gene: ", gene, " p-values (random permutation test): ", p)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS'].format(gene=gene)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            if len(conf[1]) == 0:
                compute_heatmap(ranking, gene, tgt_fp)
            else:
                compute_heatmap(ranking, gene, tgt_fp, size=2, xtickslabel=['3h', '5h'], ytickslabel=['3h', '5h'])

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_CS_HISTOGRAM']
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.bar_profile(css, tgt_fp, constants.analysis_config['PLOT_COLORS'])
        logger.info("Colocalization score p-values are: {}", p_vals)

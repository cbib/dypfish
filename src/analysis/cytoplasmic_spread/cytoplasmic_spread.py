#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import collections
import plot
import numpy as np
from image_set import ImageSet
from helpers import open_repo
from path import global_root_dir
import helpers
import itertools

def build_cytoplasmic_statistics(analysis_repo, statistics_type, molecule_type, genes, keyorder):
    gene2stat, gene2median, gene2error, gene2confidence_interval = {}, {}, {}, {}

    for gene in genes:
        logger.info("Running {} cytoplasmic spread analysis for {}", molecule_type, gene)
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)])
        if statistics_type == 'centrality':
            if molecule_type == 'mrna':
                gene2stat[gene] = image_set.compute_cytoplasmic_spots_centrality()
            else:
                gene2stat[gene] = image_set.compute_intensities_cytoplasmic_centrality()
        if statistics_type == 'spread':
            if molecule_type == 'mrna':
                gene2stat[gene] = image_set.compute_cytoplasmic_spots_spread()
            else:
                gene2stat[gene] = image_set.compute_intensities_cytoplasmic_spread()
        gene2median[gene] = np.median(gene2stat[gene])
        gene2error[gene] = helpers.sem(gene2stat[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2stat[gene])
        gene2confidence_interval[gene] = [lower, higher]

    gene2median = collections.OrderedDict(sorted(gene2median.items(),
                                                             key=lambda i: keyorder.index(i[0])))
    gene2error = collections.OrderedDict(sorted(gene2error.items(), key=lambda i: keyorder.index(i[0])))
    gene2confidence_interval = collections.OrderedDict(sorted(gene2confidence_interval.items(),
                                                              key=lambda i: keyorder.index(i[0])))
    return gene2median, gene2stat, gene2error, gene2confidence_interval


def plot_bar_profile_median_and_violin(statistics_type, molecule_type,
                                       medians, all_values, errors, confidence_intervals):
    if statistics_type == 'spread':
        tgt_image_name = constants.analysis_config['FIGURE_NAME_SPREAD'].format(molecule_type=molecule_type)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
    if statistics_type == 'centrality':
        tgt_image_name = constants.analysis_config['FIGURE_NAME_CENTRALITY'].format(molecule_type=molecule_type)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)

    if molecule_type == 'mrna':
        xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    else:
        xlabels = constants.analysis_config['PROTEINS_LABEL']

    # generate the bar profile plot
    plot.bar_profile_median(medians, errors.values(), molecule_type,
                            xlabels, tgt_fp, confidence_intervals,
                            annot=False, data_to_annot=None)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])

    # generate violin plot image
    if statistics_type == 'spread':
        tgt_image_name = constants.analysis_config['FIGURE_NAME_SPREAD_VIOLIN'].format(molecule_type=molecule_type)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
    if statistics_type == 'centrality':
        tgt_image_name = constants.analysis_config['FIGURE_NAME_CENTRALITY_VIOLIN'].format(molecule_type=molecule_type)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)

    plot.violin_profile(all_values, tgt_fp, xlabels, rotation=0, annot=False)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])


# Figure 6B top left panel : mRNA cytoplasmic spread arhgdia and arhgdia nocodazole
# Figure 6B top right panel : protein cytoplasmic spread arhgdia and arhgdia nocodazole
# Figure 6B bottom left panel : mRNA cytoplasmic spread pard3 and pard3 nocodazole
# Figure 6B bottom right panel : protein cytoplasmic spread pard3 and pard3 nocodazole
# Figure S4D : mRNA cytoplasmic spread arhgdia prrc2c
# Figure S4D : protein cytoplasmic spread arhgdia prrc2c
# Figure S6B left: mRNA cytoplasmic spread for cytod
# Figure S6B right : protein cytoplasmic spread for cytod


configurations = [
    ["src/analysis/cytoplasmic_spread/config_original.json", ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"], "Gene"],
    ["src/analysis/cytoplasmic_spread/config_nocodazole_arhgdia.json", ["arhgdia", "arhgdia_nocodazole"], "Gene"],
    ["src/analysis/cytoplasmic_spread/config_nocodazole_pard3.json", ["pard3", "pard3_nocodazole"], "Gene"],
    ["src/analysis/cytoplasmic_spread/config_prrc2c.json", ["arhgdia/control", "arhgdia/prrc2c_depleted"], "Timepoint"],
    ["src/analysis/cytoplasmic_spread/config_cytod.json", ["arhgdia_control", "arhgdia_cytod"], "Gene"]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        keyorder = conf[1]
        for molecule_type, statistics_type in itertools.product(['mrna', 'protein'],
                                                                ['centrality', 'spread']):
            if molecule_type == 'mrna':
                molecules = constants.analysis_config['MRNA_GENES']
            else:
                molecules = constants.analysis_config['PROTEINS']
            medians, spread, err, CI = build_cytoplasmic_statistics(repo, statistics_type, molecule_type,
                                                                    molecules, keyorder)
            plot_bar_profile_median_and_violin(statistics_type, molecule_type, medians, spread, err, CI)

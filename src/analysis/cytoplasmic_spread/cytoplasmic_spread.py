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


def build_cytoplasmic_spread_statistics(analysis_repo, molecule_type, genes, keyorder):
    gene2cyto_spread = {}
    gene2median_cyto_spread = {}
    gene2error = {}
    gene2confidence_interval = {}

    for gene in genes:
        logger.info("Running {} cytoplasmic spread analysis for {}", molecule_type, gene)
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)])
        if molecule_type == 'mrna':
            gene2cyto_spread[gene] = image_set.compute_cytoplasmic_spots_centrality()
        else:
            gene2cyto_spread[gene] = image_set.compute_intensities_cytoplasmic_centrality()
        gene2median_cyto_spread[gene] = np.median(gene2cyto_spread[gene])
        gene2error[gene] = helpers.sem(gene2cyto_spread[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2cyto_spread[gene])
        gene2confidence_interval[gene] = [lower, higher]

    # generate bar plot image

    gene2median_cyto_spread = collections.OrderedDict(sorted(gene2median_cyto_spread.items(),
                                                             key=lambda i: keyorder.index(i[0])))
    gene2error = collections.OrderedDict(sorted(gene2error.items(), key=lambda i: keyorder.index(i[0])))
    gene2confidence_interval = collections.OrderedDict(sorted(gene2confidence_interval.items(),
                                                              key=lambda i: keyorder.index(i[0])))
    return gene2median_cyto_spread, gene2cyto_spread, gene2error, gene2confidence_interval


def plot_bar_profile_median_and_violin(molecule_type, median_spread, spread, errors,
                                       confidence_intervals):
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    if molecule_type == 'mrna':
        xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    else:
        xlabels = constants.analysis_config['PROTEINS_LABEL']

    # generate the bar profile plot
    plot.bar_profile_median(median_spread, errors.values(), molecule_type,
                            xlabels, tgt_fp, confidence_intervals,
                            annot=False, data_to_annot=None)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.violin_profile(spread, tgt_fp, xlabels, rotation=0, annot=False)
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
        for molecule_type in ['mrna', 'protein']:
            if molecule_type == 'mrna':
                molecules = constants.analysis_config['MRNA_GENES']
            else:
                molecules = constants.analysis_config['PROTEINS']
            medians, spread, err, CI = build_cytoplasmic_spread_statistics(repo, molecule_type,
                                                                           molecules, keyorder)
            plot_bar_profile_median_and_violin(molecule_type, medians, spread, err, CI)

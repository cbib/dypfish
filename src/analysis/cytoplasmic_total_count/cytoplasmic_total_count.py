#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
import helpers
from image_set import ImageSet
from helpers import open_repo
from path import global_root_dir
import collections


def mrna_cytoplasmic_total_count(analysis_repo, keyorder):
    gene2image_set = {}
    gene2cyto_count = {}
    gene2median_cyto_count = {}
    gene2error = {}
    gene2confidence_interval = {}

    for gene in constants.analysis_config['MRNA_GENES']:
        logger.info("Running mrna cytoplasmic total count analysis for {}", gene)
        gene2image_set[gene] = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2cyto_count[gene] = gene2image_set[gene].compute_cytoplasmic_spots_counts()
        gene2median_cyto_count[gene] = np.median(gene2cyto_count[gene])
        gene2error[gene] = helpers.sem(gene2cyto_count[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2cyto_count[gene])
        gene2confidence_interval[gene] = [lower, higher]

    # generate bar plot image

    gene2median_cyto_count = collections.OrderedDict(sorted(gene2median_cyto_count.items(), key=lambda i: keyorder.index(i[0])))
    gene2error = collections.OrderedDict(sorted(gene2error.items(), key=lambda i: keyorder.index(i[0])))
    gene2confidence_interval = collections.OrderedDict(sorted(gene2confidence_interval.items(), key=lambda i: keyorder.index(i[0])))
    xlabels = constants.analysis_config['MRNA_GENES_LABEL']

    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile_median(gene2median_cyto_count,
                            gene2error.values(),
                            'mrna',
                            xlabels,
                            tgt_fp,
                            gene2confidence_interval.values(),
                            annot=False,
                            data_to_annot=gene2cyto_count
                            )

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type="mrna")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.violin_profile(gene2cyto_count, tgt_fp, xlabels, rotation=0, annot=False)


def intensities_cytoplasmic_total_count(analysis_repo, keyorder):
    gene2cyto_count = {}
    gene2median_cyto_count = {}
    gene2error = {}
    gene2confidence_interval = {}
    for gene in constants.analysis_config['PROTEINS']:
        logger.info("Running protein cytoplasmic total count analysis for {}", gene)
        imageset = ImageSet(analysis_repo, ['protein/%s/' % gene])
        gene2cyto_count[gene] = imageset.compute_cytoplasmic_intensities()
        gene2median_cyto_count[gene] = np.median(gene2cyto_count[gene])
        gene2error[gene] = helpers.sem(gene2cyto_count[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2cyto_count[gene])
        gene2confidence_interval[gene] = [lower, higher]

    # generate bar plot image

    gene2median_cyto_count = collections.OrderedDict(sorted(gene2median_cyto_count.items(), key=lambda i: keyorder.index(i[0])))
    gene2error = collections.OrderedDict(sorted(gene2error.items(), key=lambda i: keyorder.index(i[0])))
    gene2confidence_interval = collections.OrderedDict(sorted(gene2confidence_interval.items(), key=lambda i: keyorder.index(i[0])))
    xlabels = constants.analysis_config['PROTEINS_LABEL']

    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile_median(gene2median_cyto_count,
                            gene2error.values(),
                            'proteins',
                            xlabels,
                            tgt_fp,
                            gene2confidence_interval.values(),
                            annot=False,
                            data_to_annot=gene2cyto_count
                            )

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type="protein")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    plot.violin_profile(gene2cyto_count, tgt_fp, xlabels, rotation=0, annot=False)


''' 
Figure 4A bottom panel: arhgdia and arhgdia prrc2c cytoplasmic total count
Figure S6A top panel: arhgdia and arhgdia nocodazole cytoplasmic total count
Figure S6A bottom panel: pard3 and pard3 nocodazole cytoplasmic total count
Figure S6A bottom panel: arhgdia and arhgdia CytoD cytoplasmic total count
'''

configurations = [
    ["src/analysis/cytoplasmic_total_count/config_prrc2c.json", ["arhgdia/control", "arhgdia/prrc2c_depleted"], "Timepoint"],
    ["src/analysis/cytoplasmic_total_count/config_nocodazole_arhgdia.json", ["arhgdia", "arhgdia_nocodazole"], "Gene"],
    ["src/analysis/cytoplasmic_total_count/config_nocodazole_pard3.json", ["pard3", "pard3_nocodazole"], "Gene"],
    ["src/analysis/cytoplasmic_total_count/config_cytod.json", ["arhgdia_control", "arhgdia_cytod"], "Gene"]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        analysis_repo = open_repo()
        keyorder = conf[1]
        mrna_cytoplasmic_total_count(analysis_repo, keyorder)
        intensities_cytoplasmic_total_count(analysis_repo, keyorder)

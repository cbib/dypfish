#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import plot
from image_set import ImageSet
from helpers import open_repo
from path import global_root_dir


def mrna_cytoplasmic_total_count(analysis_repo):
    gene2image_set = {}
    gene2cyto_count = {}
    for gene in constants.analysis_config['MRNA_GENES']:
        logger.info("Running mrna cytoplasmic total count analysis for {}", gene)
        gene2image_set[gene] = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2cyto_count[gene] = gene2image_set[gene].compute_cytoplasmic_spots_counts()

    # generate image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


def intensities_cytoplasmic_total_count(analysis_repo):
    gene2cyto_count = {}
    for gene in constants.analysis_config['PROTEINS']:
        logger.info("Running mrna cytoplasmic total count analysis for {}", gene)
        imageset = ImageSet(analysis_repo, ['protein/%s/' % gene])
        gene2cyto_count[gene] = imageset.compute_cytoplasmic_intensities()

    # generate image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


''' 
Figure 4A bottom panel: arhgdia and arhgdia prrc2c cytoplasmic total count
Figure S6A top panel: arhgdia and arhgdia nocodazole cytoplasmic total count
Figure S6A bottom panel: pard3 and pard3 nocodazole cytoplasmic total count
Figure S6A bottom panel: arhgdia and arhgdia CytoD cytoplasmic total count
'''

configurations = [
    ["src/analysis/cytoplasmic_total_count/config_prrc2c.json", "arhgdia", "prrc2c_depleted", "Timepoint"],
    ["src/analysis/cytoplasmic_total_count/config_nocodazole_arhgdia.json", "arhgdia", "Nocodazole+", "Gene"],
    ["src/analysis/cytoplasmic_total_count/config_nocodazole_pard3.json", "pard3", "Nocodazole+", "Gene"],
    ["src/analysis/cytoplasmic_total_count/config_cytod.json", "arhgdia", "cytod+", "Gene"]
]

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        analysis_repo = open_repo()
        mrna_cytoplasmic_total_count(analysis_repo)
        intensities_cytoplasmic_total_count(analysis_repo)

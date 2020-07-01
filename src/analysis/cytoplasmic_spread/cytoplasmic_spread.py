#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger

import constants
import plot
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
from path import global_root_dir


configurations = [

    ["src/analysis/cytoplasmic_spread/config_nocodazole_arhgdia.json", "arhgdia", "Nocodazole+", "Gene"],
    #["src/analysis/peripheral_fraction_profile/config_chx.json", "arhgdia", "chx+", "Gene"],
    ["src/analysis/cytoplasmic_spread/config_nocodazole_pard3.json", "pard3", "Nocodazole+", "Gene"],
    ["src/analysis/cytoplasmic_spread/config_prrc2c.json", "arhgdia", "prrc2c_depleted", "Timepoint"]
    ["src/analysis/cytoplasmic_spread/config_cytod.json", "arhgdia", "cytod+", "Gene"]
]

# Figure 6B top left panel : mRNA cytoplasmic spread arhgdia and arhgdia nocodazole
# Figure 6B top right panel : protein cytoplasmic spread arhgdia and arhgdia nocodazole
# Figure 6B bottom left panel : mRNA cytoplasmic spread pard3 and pard3 nocodazole
# Figure 6B bottom right panel : protein cytoplasmic spread pard3 and pard3 nocodazole
# Figure S4D : mRNA cytoplasmic spread arhgdia prrc2c
# Figure S4D : protein cytoplasmic spread arhgdia prrc2c
# Figure S6B left: mRNA cytoplasmic spread for cytod
# Figure S6B right : protein cytoplasmic spread for cytod

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        analysis_repo = open_repo()

        logger.info("Running mrna cytoplasmic spread analysis")
        gene2image_set = {}
        gene2cyto_count = {}
        for gene in constants.analysis_config['MRNA_GENES']:
            gene2image_set[gene] = ImageSet(analysis_repo, ['mrna/%s/' % gene])
            gene2cyto_count[gene] = gene2image_set[gene].compute_spots_cytoplasmic_spread()

        # generate image
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
        logger.info("Generated image at {}", tgt_fp)

        logger.info("Running protein cytoplasmic spread analysis")
        gene2image_set = {}
        gene2cyto_count = {}
        for gene in constants.analysis_config['PROTEINS']:
            gene2image_set[gene] = ImageSet(analysis_repo, ['protein/%s/' % gene])
            gene2cyto_count[gene] = gene2image_set[gene].compute_intensities_cytoplasmic_spread()

        # generate image
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
        logger.info("Generated image at {}", tgt_fp)
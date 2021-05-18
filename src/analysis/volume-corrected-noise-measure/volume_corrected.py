#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import tqdm
from loguru import logger

import constants
import plot
from helpers import open_repo
from image_set import ImageSet
from path import global_root_dir

# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/volume-corrected-noise-measure/config_original.json"],
    ["src/analysis/volume-corrected-noise-measure/config_standard.json"]
]

if __name__ == '__main__':
    for conf in configurations:

        stat_annotations = False
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        molecule_type = 'mrna'
        plot_colors = constants.analysis_config['PLOT_COLORS']
        if "original" in conf[0]:
            # Figure S1.B The volume-corrected noise measure (Padovan-Merhar et al., 2015) across time for 6 mRNAs was compared.
            for i, gene in enumerate(tqdm.tqdm(constants.dataset_config['MRNA_GENES'], desc="Genes")):
                nms = []
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BARPLOT'].format(gene=gene)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                      tgt_image_name)
                for timepoint in tqdm.tqdm(constants.dataset_config['TIMEPOINTS_MRNA'], desc="Timepoint"):
                    imageset = ImageSet(repo, ["{0}/{1}/{2}/".format(molecule_type, gene, timepoint)])
                    nm = imageset.compute_surface_corrected_nm()
                    nms.append(nm)
                plot.bar_profile(nms, tgt_fp, plot_colors[i])
                logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])
        else:
            timepoint = '3h'
            nms = []
            for i, gene in enumerate(constants.analysis_config['MRNA_GENES']):
                imageset = ImageSet(repo, ['{0}/{1}/{2}/'.format(molecule_type, gene, timepoint)])
                nms.append(imageset.compute_volume_corrected_nm())
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_VOLUME']
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.bar_profile(nms, tgt_fp, plot_colors)
                logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])
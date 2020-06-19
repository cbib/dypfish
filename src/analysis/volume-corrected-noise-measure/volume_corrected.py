#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import tqdm
import constants
from plot import histogram_noise_measured, bar_profile_simple
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,
                                                           "src/analysis/volume-corrected-noise-measure/config_original.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

imageset = ImageSet(analysis_repo, ["mrna/arhgdia/3h/"])
nm_arhgdia = imageset.compute_surface_corrected_nm()
imageset = ImageSet(analysis_repo, ["mrna/arhgdia_cultured/"])
nm_arhgdia_cultured = imageset.compute_surface_corrected_nm()
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_SURFACE']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
histogram_noise_measured(nm_arhgdia, nm_arhgdia_cultured, tgt_fp)

# imageset = ImageSet(analysis_repo, ["mrna/arhgdia/3h/"])
# nm_arhgdia = imageset.compute_volume_corrected_nm()
# imageset = ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"])
# nm_arhgdia_cultured = imageset.compute_volume_corrected_nm()
# tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_VOLUME']
# tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
# histogram_noise_measured(nm_arhgdia, nm_arhgdia_cultured, tgt_fp)

plot_colors = constants.analysis_config['PLOT_COLORS']
for i, gene in enumerate(tqdm.tqdm(constants.dataset_config['MRNA_GENES'], desc="Genes")):
    nms = []
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BARPLOT'].format(gene=gene)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    for timepoint in tqdm.tqdm(constants.dataset_config['TIMEPOINTS_MRNA'], desc="Timepoint"):
        imageset = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
        nm = imageset.compute_surface_corrected_nm()
        nms.append(nm)
    bar_profile_simple(nms, tgt_fp, plot_colors[i])
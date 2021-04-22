#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import tqdm
import constants
from plot import bar_profile
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint

constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,
                                                           "src/analysis/volume-corrected-noise-measure/config_original.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
plot_colors = constants.analysis_config['PLOT_COLORS']



# Figure S1.D bottom left volume-corrected noise measure Standard vs Micropatterned
imageset = ImageSet(analysis_repo, ["mrna/arhgdia/3h/"])
nm_arhgdia = imageset.compute_volume_corrected_nm()
imageset = ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"])
nm_arhgdia_cultured = imageset.compute_volume_corrected_nm()
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_VOLUME']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
bar_profile([nm_arhgdia, nm_arhgdia_cultured], tgt_fp, plot_colors)

# Figure S1.B The volume-corrected noise measure (Padovan-Merhar et al., 2015) across time for 6 mRNAs was compared.
for i, gene in enumerate(tqdm.tqdm(constants.dataset_config['MRNA_GENES'], desc="Genes")):
    nms = []
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BARPLOT'].format(gene=gene)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    for timepoint in tqdm.tqdm(constants.dataset_config['TIMEPOINTS_MRNA'], desc="Timepoint"):
        imageset = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
        nm = imageset.compute_surface_corrected_nm()
        nms.append(nm)
    bar_profile(nms, tgt_fp, plot_colors[i])

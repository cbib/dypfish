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

# Figure S5B left: mrna cytoplasmic spread for cytod
constants.init_config(
    analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_spread/config_cytod.json"))

dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

logger.info("Running mrna cytoplasmic spread analysis for Cytod data")
gene2image_set = {}
gene2cyto_count = {}
for gene in constants.analysis_config['MRNA_GENES']:
    gene2image_set[gene] = ImageSet(analysis_repo, ['mrna/%s/' % gene])
    gene2cyto_count[gene] = gene2image_set[gene].compute_spots_cytoplasmic_spread()
# generate image
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
logger.info("Generated image at {}", tgt_fp)

# Figure 5B right : protein cytoplasmic spread for cytod
logger.info("Running protein cytoplasmic spread analysis for cytod data")
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

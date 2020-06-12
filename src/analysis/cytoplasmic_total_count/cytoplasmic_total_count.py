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


def mrna_cytoplasmic_total_count():
    gene2image_set = {}
    gene2cyto_count = {}
    for gene in constants.analysis_config['MRNA_GENES']:
        gene2image_set[gene] = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2cyto_count[gene] = gene2image_set[gene].compute_cytoplasmic_spots_counts()

    # generate image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


def intensities_cytoplasmic_total_count():
    gene2image_set = {}
    gene2cyto_count = {}
    for gene in constants.analysis_config['PROTEINS']:
        gene2image_set[gene] = ImageSet(analysis_repo, ['protein/%s/' % gene])
        gene2cyto_count[gene] = gene2image_set[gene].compute_cytoplasmic_intensities()

    # generate image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2cyto_count.values(), gene2cyto_count.keys(), figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


# Figure S5A upper panel: arhgdia and arhgdia nocodazole cytoplasmic total count
# this should be called as soon as possible
logger.info("Running cytoplasmic total count analysis for Nocodazole")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_total_count/config_nocodazole_arhgdia.json"))

dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])

analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_cytoplasmic_total_count()
intensities_cytoplasmic_total_count()

# Figure S5A : arhgdia and arhgdia nocodazole cytoplasmic total count
# this should be called as soon as possible
logger.info("Running cytoplasmic total count analysis for Nocodazole")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_total_count/config_nocodazole_pard3.json"))

dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])

analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_cytoplasmic_total_count()
intensities_cytoplasmic_total_count()

# Figure S5A : arhgdia and arhgdia CytoD cytoplasmic total count
logger.info("Running cytoplasmic total count analysis for CytoD")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_total_count/config_cytod.json"))

dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])

analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_cytoplasmic_total_count()
intensities_cytoplasmic_total_count()

# Figure ?? : arhgdia and arhgdia prrc2c cytoplasmic total count
logger.info("Running cytoplasmic total count analysis for prrc2c")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_total_count/config_prrc2c.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])

analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_cytoplasmic_total_count()
intensities_cytoplasmic_total_count()


# Figure ?? : arhgdia and arhgdia original cytoplasmic total count
logger.info("Running cytoplasmic total count analysis for prrc2c")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/cytoplasmic_total_count/config_prrc2c.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])

analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
mrna_cytoplasmic_total_count()
intensities_cytoplasmic_total_count()

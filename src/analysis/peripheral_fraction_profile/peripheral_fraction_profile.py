#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import pprint as pp

from loguru import logger

import constants
import plot
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def plot_mrna_peripheral_fraction_profiles():
    dataset_root_fp = pathlib.Path(
        constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

    gene2profile_mrna_periph_fraction = []
    for gene in constants.analysis_config['MRNA_GENES']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2profile_mrna_periph_fraction.append(
            np.average(image_set.compute_histogram_spots_peripheral_fraction(), axis=0))
    # generate image

    # normalized by gapdh profile
    gene2profile_mrna_periph_fraction = gene2profile_mrna_periph_fraction / np.matlib.repmat(
        gene2profile_mrna_periph_fraction[2], len(constants.dataset_config['MRNA_GENES']), 1)

    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.profile(gene2profile_mrna_periph_fraction, constants.dataset_config['MRNA_GENES'],
                 constants.analysis_config['NUM_CONTOURS'], figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)

# TODO : code redundancy below
def plot_mrna_histogram_peripheral_fraction():
    dataset_root_fp = pathlib.Path(
        constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

    gene2histogram_mrna_periph_fraction = {}
    for gene in constants.analysis_config['MRNA_GENES']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2histogram_mrna_periph_fraction[gene] = image_set.compute_histogram_spots_peripheral_counts()
    # generate image
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_HISTOGRAM_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2histogram_mrna_periph_fraction.values(), gene2histogram_mrna_periph_fraction.keys(),
                     figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)

def plot_mrna_histogram_peripheral_fraction_2D():
    dataset_root_fp = pathlib.Path(
        constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

    gene2histogram_mrna_periph_fraction = {}
    for gene in constants.analysis_config['MRNA_GENES']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene],force2D=True)
        gene2histogram_mrna_periph_fraction[gene] = image_set.compute_histogram_spots_peripheral_counts()
    # generate image
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_HISTOGRAM_FORMAT_2D'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2histogram_mrna_periph_fraction.values(), gene2histogram_mrna_periph_fraction.keys(),
                     figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


def plot_mrna_intensity_histogram_peripheral_fraction():
    dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

    gene2histogram_protein_periph_fraction = {}
    for gene in constants.analysis_config['PROTEINS']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2histogram_protein_periph_fraction[gene] = image_set.compute_histogram_intensities_peripheral_fractions()
    logger.debug("Intensities fractions {}", gene2histogram_protein_periph_fraction)
    # generate image
    fraction = fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2histogram_protein_periph_fraction.values(), gene2histogram_protein_periph_fraction.keys(),
                     figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)

def plot_protein_histogram_peripheral_fraction():
    dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
    gene2histogram_protein_periph_fraction = {}
    for gene in constants.analysis_config['PROTEINS']:
        image_set = ImageSet(analysis_repo, ['protein/%s/' % gene])
        gene2histogram_protein_periph_fraction[gene] = image_set.compute_histogram_intensities_peripheral_fractions()
    logger.debug("Intensities fractions {}", gene2histogram_protein_periph_fraction)
    # generate image
    fraction = fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.bar_profile(gene2histogram_protein_periph_fraction.values(), gene2histogram_protein_periph_fraction.keys(),
                     figname=tgt_fp)
    logger.info("Generated image at {}", tgt_fp)


# Figure S2B : mRNA peripheral fraction at 30% for beta-actin, arhgdia, gapdh, pard3, pkp4, rab13
# to obtain 10% plots, change the config_original file accordingly: PERIPHERAL_FRACTION_THRESHOLD = 10
logger.info("Peripheral fraction histogram for mRNA the original data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,"src/analysis/peripheral_fraction_profile/config_original.json"))
plot_mrna_histogram_peripheral_fraction()

# Figure 2B : Comparison of the enrichment of 5 mRNAs with respect to Gapdh mRNA in a peripheral cellular-
# region whose width varies from 0-100% of the radial distance from the plasma membrane to the nucleus.
logger.info("Peripheral fraction profile for the mRNA original data")
plot_mrna_peripheral_fraction_profiles()

# Figure 2G upper panel: protein peripheral fraction at 30% for CHX arhgdia and pard3
logger.info("Peripheral fraction for the protein CHX data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/peripheral_fraction_profile/config_chx.json"))
plot_protein_histogram_peripheral_fraction()

# Figure S2C : mrna treated as intensity data peripheral fraction at 30% for CHX arhgdia and pard3
# to obtain 10% plots, change the config_original file accordingly: PERIPHERAL_FRACTION_THRESHOLD = 10
logger.info("Peripheral fraction for the mRNA CHX data")
plot_mrna_intensity_histogram_peripheral_fraction()

# Figure 5C top panel: 2D data : peripheral fraction nocodazole arhgdia, pard3 at 30% -> barplot mrna and protein

constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,"src/analysis/peripheral_fraction_profile/config_nocodazole_arhgdia.json"))
plot_mrna_histogram_peripheral_fraction_2D()
plot_protein_histogram_peripheral_fraction()

# Figure 5C bottom panel: 2D data : peripheral fraction nocodazole arhgdia, pard3 at 30% -> barplot mrna and protein
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,"src/analysis/peripheral_fraction_profile/config_nocodazole_pard3.json"))
plot_mrna_histogram_peripheral_fraction_2D()
plot_protein_histogram_peripheral_fraction()

# Figure S5C : mRNA peripheral fraction at 30% for arhgdia, cytod
# to obtain 10% plots, change the config_original file accordingly: PERIPHERAL_FRACTION_THRESHOLD = 10
logger.info("Peripheral fraction for the cytod data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,"src/analysis/peripheral_fraction_profile/config_cytod.json"))
plot_mrna_histogram_peripheral_fraction()
plot_protein_histogram_peripheral_fraction()

# New Figure : mRNA peripheral fraction at 30% for arhgdia for control and prrc2c depleted
# to obtain 10% plots, change the config_original file accordingly: PERIPHERAL_FRACTION_THRESHOLD = 10
logger.info("Peripheral fraction for the prrc2c data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir,"src/analysis/peripheral_fraction_profile/config_prrc2c.json"))
plot_mrna_histogram_peripheral_fraction()
plot_protein_histogram_peripheral_fraction()
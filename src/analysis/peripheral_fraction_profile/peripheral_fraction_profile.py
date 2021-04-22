#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import pprint as pp

from loguru import logger

import constants
import plot
import pandas as pd
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
from helpers import open_repo
# this should be called as soon as possible
from path import global_root_dir
import helpers


def plot_mrna_peripheral_fraction_profiles(analysis_repo):
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
def plot_mrna_histogram_peripheral_fraction(analysis_repo, annot=False, force2D=False):
    gene2mrna_periph_fraction = {}
    gene2median_mrna_periph_fraction = {}
    gene2error = {}
    gene2confidence_interval = {}
    for gene in constants.analysis_config['MRNA_GENES']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene], force2D=force2D)
        gene2mrna_periph_fraction[gene] = image_set.compute_histogram_spots_peripheral_counts()
        gene2median_mrna_periph_fraction[gene] = np.median(gene2mrna_periph_fraction[gene])
        gene2error[gene] = helpers.sem(gene2mrna_periph_fraction[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2mrna_periph_fraction[gene])
        gene2confidence_interval[gene] = [lower, higher]
    # generate bar plot image
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_HISTOGRAM_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)

    xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    plot.bar_profile_median(gene2median_mrna_periph_fraction,
                            gene2error.values(),
                            'mrna',
                            xlabels,
                            tgt_fp,
                            gene2confidence_interval.values(),
                            annot=annot,
                            data_to_annot=gene2mrna_periph_fraction
                            )

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
    xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    plot.violin_profile(gene2mrna_periph_fraction, tgt_fp, xlabels, annot=annot)


def plot_mrna_intensity_histogram_peripheral_fraction(analysis_repo, annot=False):
    '''
    Specific call to compute_histogram_intensities_peripheral_fractions for IF data treated as FISH.
    '''
    gene2protein_periph_fraction = {}
    gene2median_protein_periph_fraction = {}
    gene2error = {}
    gene2confidence_interval = {}
    for gene in constants.analysis_config['PROTEINS']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2protein_periph_fraction[gene] = image_set.compute_histogram_intensities_peripheral_fractions()
        gene2median_protein_periph_fraction[gene] = np.median(gene2protein_periph_fraction[gene])
        gene2error[gene] = helpers.sem(gene2protein_periph_fraction[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2protein_periph_fraction[gene])
        gene2confidence_interval[gene] = [lower, higher]
    logger.debug("Intensities fractions {}", gene2protein_periph_fraction)
    # generate image
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    plot.bar_profile_median(gene2median_protein_periph_fraction,
                            gene2error.values(),
                            'mrna',
                            xlabels,
                            tgt_fp,
                            gene2confidence_interval.values(),
                            annot=annot,
                            data_to_annot=gene2protein_periph_fraction
                            )

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type="mrna", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
    xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    plot.violin_profile(gene2protein_periph_fraction, tgt_fp, xlabels, annot=annot)


def plot_protein_histogram_peripheral_fraction(analysis_repo, annot=False):
    gene2protein_periph_fraction = {}
    gene2median_protein_periph_fraction = {}
    gene2error = {}
    gene2confidence_interval = {}
    for gene in constants.analysis_config['PROTEINS']:
        image_set = ImageSet(analysis_repo, ['protein/%s/' % gene])
        gene2protein_periph_fraction[gene] = image_set.compute_histogram_intensities_peripheral_fractions()
        gene2median_protein_periph_fraction[gene] = np.median(gene2protein_periph_fraction[gene])
        gene2error[gene] = helpers.sem(gene2protein_periph_fraction[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2protein_periph_fraction[gene])
        gene2confidence_interval[gene] = [lower, higher]

    # generate image
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="protein", fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    xlabels = constants.analysis_config['PROTEINS_LABEL']
    plot.bar_profile_median(gene2median_protein_periph_fraction,
                            gene2error.values(),
                            'protein',
                            xlabels,
                            tgt_fp,
                            gene2confidence_interval.values(),
                            annot=annot,
                            data_to_annot=gene2protein_periph_fraction
                            )

    # generate violin plot image
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type="protein",
                                                                                   fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    xlabels = constants.analysis_config['PROTEINS_LABEL']
    plot.violin_profile(gene2protein_periph_fraction, tgt_fp, xlabels, annot=annot)


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/peripheral_fraction_profile/config_original.json", ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/peripheral_fraction_profile/config_chx.json", ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/peripheral_fraction_profile/config_cytod.json", ["arhgdia", "cytod+"]],
    ["src/analysis/peripheral_fraction_profile/config_prrc2c.json", ['arhgdia/control', "arhgdia/prrc2c_depleted"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_arhgdia.json", ["arhgdia", "Nocodazole+"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_pard3.json", ["pard3", "Nocodazole+"]]
]

# TODO recuperer commentaires indiquant nom de figure généré, le passer en argument des fonctions de plot
if __name__ == '__main__':
    for conf in configurations:
        stat_annotations = False
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        if "original" in conf[0]:
            logger.info("Peripheral fraction profile for the mRNA original data")
            plot_mrna_peripheral_fraction_profiles(repo)
            logger.info("Peripheral fraction histogram for mRNA the original data")
            plot_mrna_histogram_peripheral_fraction(repo, annot=stat_annotations)
        elif "chx" in conf[0]:
            logger.info("Peripheral fraction for the protein CHX data")
            plot_protein_histogram_peripheral_fraction(repo, annot=stat_annotations)
            logger.info("Peripheral fraction for the mRNA CHX data")
            plot_mrna_intensity_histogram_peripheral_fraction(repo, annot=stat_annotations)
        elif "nocodazole" in conf[0]:
            logger.info("2D Peripheral fraction for the mRNA nocodazole data")
            plot_mrna_histogram_peripheral_fraction(repo, annot=stat_annotations, force2D=True)
            logger.info("Peripheral fraction for the mRNA nocodazole data")
            plot_protein_histogram_peripheral_fraction(repo, annot=stat_annotations)
        else:
            logger.info("Peripheral fraction for the prrc2c data")
            plot_mrna_histogram_peripheral_fraction(repo, annot=stat_annotations)
            plot_protein_histogram_peripheral_fraction(repo, annot=stat_annotations)



#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

# This analysis will produce incoherent results if PERIPHERAL_FRACTION value provided in
# the analysis_config.json file goes beyond the cytoplasm and within the nucleus

import collections
import pathlib

import numpy as np
from loguru import logger

import constants
import helpers
import plot
from helpers import open_repo
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir


def build_mrna_peripheral_fraction_profiles(analysis_repo, normalisation_gene=None):
    genes= constants.analysis_config['MRNA_GENES']
    gene2m_fractions = {}
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        peripheral_fractions = image_set.compute_cytoplasmic_spots_fractions_per_periphery()
        gene2m_fractions[gene] = np.median(peripheral_fractions, axis=0)

    # normalized by gapdh profile
    if normalisation_gene:
        for gene in genes:
            gene2m_fractions[gene]  = gene2m_fractions[gene] / gene2m_fractions[normalisation_gene]

    fractions = collections.OrderedDict(sorted(gene2m_fractions.items(), key=lambda i: keyorder.index(i[0])))
    return fractions


def build_histogram_peripheral_fraction(analysis_repo, molecule_type, keyorder, force2D=False):
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    gene2fractions = {}
    gene2median = {}
    gene2error = {}
    gene2ci = {}
    if molecule_type == 'mrna':
        genes = constants.analysis_config['MRNA_GENES']
    else:
        genes = constants.analysis_config['PROTEINS']
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)], force2D=force2D)
        if molecule_type == 'mrna':
            gene2fractions[gene] = image_set.compute_cytoplasmic_spots_fractions_per_periphery()
        else:
            gene2fractions[gene] = image_set.compute_cytoplasmic_intensities_fractions_per_periphery()

        median_pfractions = np.median(gene2fractions[gene], axis=0)
        gene2median[gene] = median_pfractions[fraction]
        gene2error[gene] = helpers.sem(median_pfractions, factor=0)
        lower, higher = helpers.median_confidence_interval(median_pfractions)
        gene2ci[gene] = [lower, higher]
        gene2fractions[gene]=gene2fractions[gene][:, fraction]

    fractions = collections.OrderedDict(sorted(gene2fractions.items(), key=lambda i: keyorder.index(i[0])))
    gene2median = collections.OrderedDict(sorted(gene2median.items(), key=lambda i: keyorder.index(i[0])))
    gene2error = collections.OrderedDict(sorted(gene2error.items(), key=lambda i: keyorder.index(i[0])))
    gene2ci = collections.OrderedDict(sorted(gene2ci.items(), key=lambda i: keyorder.index(i[0])))

    return gene2median, fractions, gene2error, gene2ci


def plot_bar_profile_median_and_violin(molecule_type, medians, fractions, errors, CI, annotations):
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_HISTOGRAM_FORMAT'].format(molecule_type=molecule_type,
                                                                                      fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    if molecule_type=='mrna':
        xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    else:
        xlabels = constants.analysis_config['PROTEINS_LABEL']

    # generate the bar profile plot
    plot.bar_profile_median(medians, errors, 'mrna', xlabels, tgt_fp, confidence_interval=CI,
                            annot=annotations, data_to_annot=fractions)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])

    # generate the violin plot
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type=molecule_type,
                                                                                   fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.violin_profile(fractions, tgt_fp, xlabels, annot=annotations)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/peripheral_fraction_profile/config_original.json",
     ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/peripheral_fraction_profile/config_chx.json",
     ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/peripheral_fraction_profile/config_cytod.json",
     ["arhgdia", "cytod+"]],
    ["src/analysis/peripheral_fraction_profile/config_prrc2c.json",
     ['arhgdia/control', "arhgdia/prrc2c_depleted"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_arhgdia.json", ["arhgdia", "Nocodazole+"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_pard3.json", ["pard3", "Nocodazole+"]]
]

if __name__ == '__main__':
    for conf in configurations:
        stat_annotations = False
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        keyorder = conf[1]
        repo = open_repo()
        if "original" in conf[0]:
            logger.info("Peripheral fraction profile for the mRNA original data")
            medians = build_mrna_peripheral_fraction_profiles(repo)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.profile(medians, constants.dataset_config['MRNA_GENES'],
                         constants.analysis_config['NUM_CONTOURS'], figname=tgt_fp)
            logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

            logger.info("Peripheral fraction histogram for mRNA the original data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna',
                                                                              keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        elif "chx" in conf[0]:
            logger.info("Peripheral fraction histogram for the protein CHX data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein',
                                                                              keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

            logger.info("Peripheral fraction histogram for the mRNA CHX data")
            # this analysis is done in 2D
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein',
                                                                              keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        elif "nocodazole" in conf[0]:
            logger.info("2D Peripheral fraction histogram for the mRNA nocodazole data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna',
                                                                              force2D=True, keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions = fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

            logger.info("Peripheral fraction histogram for the mRNA nocodazole data")
            build_histogram_peripheral_fraction(repo, molecule_type='mrna')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        else:
            logger.info("Peripheral fraction histogram for the prrc2c data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna',
                                                                              keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein')
            plot_bar_profile_median_and_violin(molecule_type='protein', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)


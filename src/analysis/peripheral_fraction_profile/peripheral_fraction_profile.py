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


def build_mrna_peripheral_fraction_profiles(analysis_repo, normalisation_gene='None'):
    genes = constants.analysis_config['MRNA_GENES']
    gene2m_fractions = {}
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        peripheral_fractions = image_set.compute_cytoplasmic_spots_fractions_per_periphery()
        gene2m_fractions[gene] = np.median(peripheral_fractions, axis=0)

    # normalized by gapdh profile
    if normalisation_gene != 'None':
        for gene in genes:
            gene2m_fractions[gene] = gene2m_fractions[gene] / gene2m_fractions[normalisation_gene]

    fractions = collections.OrderedDict(sorted(gene2m_fractions.items(), key=lambda i: keyorder.index(i[0])))
    return fractions


def build_histogram_peripheral_fraction(analysis_repo, molecule_type, keyorder, force2D=False, mRNAFromContinuousSignal=False):
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    gene2fractions, medians, erros, CI = {}, {}, {}, {}
    if molecule_type == 'mrna':
        genes = constants.analysis_config['MRNA_GENES']
    else:
        genes = constants.analysis_config['PROTEINS']
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)], force2D=force2D)
        if molecule_type == 'mrna':
            if mRNAFromContinuousSignal:
                gene2fractions[gene] = image_set.compute_cytoplasmic_intensities_fractions_per_periphery()
            else:
                gene2fractions[gene] = image_set.compute_cytoplasmic_spots_fractions_per_periphery(force2D=force2D)
        else:
            gene2fractions[gene] = image_set.compute_cytoplasmic_intensities_fractions_per_periphery()

        median_pfractions = np.median(gene2fractions[gene], axis=0)
        medians[gene] = median_pfractions[fraction]
        erros[gene] = helpers.sem(median_pfractions, factor=0)
        lower, higher = helpers.median_confidence_interval(median_pfractions)
        CI[gene] = [lower, higher]
        gene2fractions[gene] = gene2fractions[gene][:, fraction]

    fractions = collections.OrderedDict(sorted(gene2fractions.items(), key=lambda i: keyorder.index(i[0])))
    medians = collections.OrderedDict(sorted(medians.items(), key=lambda i: keyorder.index(i[0])))
    erros = collections.OrderedDict(sorted(erros.items(), key=lambda i: keyorder.index(i[0])))
    CI = collections.OrderedDict(sorted(CI.items(), key=lambda i: keyorder.index(i[0])))

    return medians, fractions, erros, CI


def plot_bar_profile_median_and_violin(molecule_type, medians, fractions, errors, CI, annotations):
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    tgt_image_name = constants.analysis_config['FIGURE_NAME_HISTOGRAM_FORMAT'].format(molecule_type=molecule_type,
                                                                                      fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)

    if molecule_type == 'mrna':
        xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    else:
        xlabels = constants.analysis_config['PROTEINS_LABEL']

    # generate the bar profile plot
    plot.bar_profile_median(medians, errors, 'mrna', xlabels, tgt_fp, confidence_interval=CI, annot=annotations, data_to_annot=fractions)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])

    # generate the violin plot
    tgt_image_name = constants.analysis_config['FIGURE_NAME_VIOLIN_FORMAT'].format(molecule_type=molecule_type,
                                                                                   fraction=fraction)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plot.violin_profile(fractions, tgt_fp, xlabels, annot=annotations)
    logger.info("Generated plot at {}", str(tgt_fp).split("analysis/")[1])


# Figure 3.A Peripheral fraction profile for the mRNA original data normalized by Gapdh
# Figure 3.F Peripheral fraction histograms for protein CHX data (30%)
# Figure 7C top left panel : Peripheral fraction histograms for mRNA arhgdia and arhgdia nocodazole
# Figure 7C top right panel : Peripheral fraction histograms for protein arhgdia and arhgdia nocodazole
# Figure 7C bottom left panel : Peripheral fraction histograms for mRNA pard3 and pard3 nocodazole
# Figure 7C bottom right panel : Peripheral fraction histograms for protein pard3 and pard3 nocodazole
# Figure S2A Peripheral fraction profile for the mRNA original data
# Figure S2C top Peripheral fraction histograms for mRNA original data (10%)
# Figure S2C bottom Peripheral fraction histograms for mRNA original data (30%)
# Figure S2D top left Peripheral fraction histograms for mRNA CHX data (10%)
# Figure S2D top right Peripheral fraction histograms for mRNA  CHX data (30%)
# Figure S2D bottom left Peripheral fraction histograms for protein CHX data (10%)
# Figure S4C top left panel : Peripheral fraction histograms for mRNA prrc2c data (10%)
# Figure S4C top right panel : Peripheral fraction histograms for protein prrc2c data (10%)
# Figure S4C bottom left panel : Peripheral fraction histograms for mRNA prrc2c data (30%)
# Figure S4C bottom right panel : Peripheral fraction histograms for protein prrc2c data (30%)
# Figure S6C left panel : Peripheral fraction histograms for mRNA cytod data (30%)
# Figure S6C right panel : Peripheral fraction histograms for protein cytod data (30%)


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/peripheral_fraction_profile/config_original.json",
     ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/peripheral_fraction_profile/config_chx.json",
     ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/peripheral_fraction_profile/config_cytod.json",
     ["arhgdia_control", "arhgdia_cytod"]],
    ["src/analysis/peripheral_fraction_profile/config_prrc2c.json",
     ['arhgdia/control', "arhgdia/prrc2c_depleted"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_arhgdia.json",
     ["arhgdia", "arhgdia_nocodazole"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_pard3.json",
     ["pard3", "pard3_nocodazole"]]
]

if __name__ == '__main__':
    for conf in configurations:
        stat_annotations = True
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        keyorder = conf[1]
        repo = open_repo()
        if "original" in conf[0]:
            logger.info("Peripheral fraction profile for the mRNA original data")
            normalisation_gene = constants.analysis_config['NORMALISATION_GENE']
            medians = build_mrna_peripheral_fraction_profiles(repo, normalisation_gene=normalisation_gene)
            medians = build_mrna_peripheral_fraction_profiles(repo)

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.profile(medians, figname=tgt_fp)

            logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

            logger.info("Peripheral fraction histogram for mRNA the original data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna',
                                                                              keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        if "chx" in conf[0]:
            logger.info("Peripheral fraction histograms for CHX data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna', keyorder=keyorder, force2D=True, mRNAFromContinuousSignal=True)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)
            # this analysis is done in 2D
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein',
                                                                              keyorder=keyorder, force2D=True)
            plot_bar_profile_median_and_violin(molecule_type='protein', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        if "nocodazole" in conf[0] or 'cytod' in conf[0]:
            logger.info("2D Peripheral fraction histograms")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna',
                                                                              force2D=True, keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein', force2D=True, keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='protein', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

        if 'prrc2c' in conf[0]:
            logger.info("Peripheral fraction histograms for prcc2c data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna', keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein', keyorder=keyorder)
            plot_bar_profile_median_and_violin(molecule_type='protein', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)

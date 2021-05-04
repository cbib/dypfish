#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

# This analysis will produce incogerent results if PERIPHERAL_FRACTION goes within the nucleus

import pathlib
from loguru import logger
import constants
import plot
import numpy as np
from image_set import ImageSet
from helpers import open_repo
# this should be called as soon as possible
from path import global_root_dir
import helpers


def build_mrna_peripheral_fraction_profiles(analysis_repo):
    genes= constants.analysis_config['MRNA_GENES']
    gene2mean_fractions = {}
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        peripheral_fractions = image_set.compute_cytoplsamic_spots_fractions_per_periphery()
        gene2mean_fractions[gene] = np.mean(peripheral_fractions, axis=0)

    # normalized by gapdh profile
    # for gene in genes:
    #     gene2mean_fractions[gene]  = gene2mean_fractions[gene] / gene2mean_fractions['gapdh']

    # this is because of the format that the plotting function expects
    fractions = np.array([list(v) for v in gene2mean_fractions.values()])
    return fractions


def build_histogram_peripheral_fraction(analysis_repo, molecule_type, force2D=False):
    fraction = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
    gene2periph_fractions = {}
    gene2median_pfraction = {}
    gene2error = {}
    gene2ci = {}
    if molecule_type == 'mrna':
        genes = constants.analysis_config['MRNA_GENES']
    else:
        genes = constants.analysis_config['PTOTEINS']
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)], force2D=force2D)
        if molecule_type == 'mrna':
            gene2periph_fractions[gene] = image_set.compute_cytoplsamic_spots_fractions_per_periphery()
        else:
            gene2periph_fractions[gene] = image_set.compute_cytoplasmic_intensities_fractions_per_periphery()
        median_pfractions = np.median(gene2periph_fractions[gene], axis=0)
        gene2median_pfraction[gene] = median_pfractions[fraction]
        gene2error[gene] = helpers.sem(median_pfractions, factor=0)
        lower, higher = helpers.median_confidence_interval(median_pfractions)
        gene2ci[gene] = [lower, higher]

    fractions = np.array([list(v) for v in gene2periph_fractions.values()])
    return gene2median_pfraction, fractions, gene2error, gene2ci

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
    ["src/analysis/peripheral_fraction_profile/config_original.json"], # ['beta_actin', 'arhgdia', 'gapdh', 'pard3', 'pkp4', 'rab13']],
    ["src/analysis/peripheral_fraction_profile/config_chx.json"],# ['arhgdia', 'arhgdia_CHX', 'pard3', 'pard3_CHX']],
    ["src/analysis/peripheral_fraction_profile/config_cytod.json"], #["arhgdia", "cytod+"]],
    ["src/analysis/peripheral_fraction_profile/config_prrc2c.json"], #['arhgdia/control', "arhgdia/prrc2c_depleted"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_arhgdia.json"], #["arhgdia", "Nocodazole+"]],
    ["src/analysis/peripheral_fraction_profile/config_nocodazole_pard3.json"], #["pard3", "Nocodazole+"]]
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
            gene2mean_fraction = build_mrna_peripheral_fraction_profiles(repo)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.profile(gene2mean_fraction, constants.dataset_config['MRNA_GENES'],
                         constants.analysis_config['NUM_CONTOURS'], figname=tgt_fp)
            logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

            logger.info("Peripheral fraction histogram for mRNA the original data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err.values(), CI=CI, annotations=stat_annotations)
            exit()

        elif "chx" in conf[0]:
            logger.info("Peripheral fraction histogram for the protein CHX data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

            logger.info("Peripheral fraction histogram for the mRNA CHX data")
            # this analysis is done in 2D
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

        elif "nocodazole" in conf[0]:
            logger.info("2D Peripheral fraction histogram for the mRNA nocodazole data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna', force2D=True)
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions = fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

            logger.info("Peripheral fraction histogram for the mRNA nocodazole data")
            build_histogram_peripheral_fraction(repo, molecule_type='mrna')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

        else:
            logger.info("Peripheral fraction histogram for the prrc2c data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein')
            plot_bar_profile_median_and_violin(molecule_type='protein', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)


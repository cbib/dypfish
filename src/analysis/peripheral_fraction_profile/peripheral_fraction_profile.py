#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

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
    gene2profile_mrna_periph_fraction = []
    for gene in constants.analysis_config['MRNA_GENES']:
        image_set = ImageSet(analysis_repo, ['mrna/%s/' % gene])
        gene2profile_mrna_periph_fraction.append(
            np.average(image_set.compute_spots_fractions_per_periphery(), axis=0))

    # normalized by gapdh profile
    gene2profile_mrna_periph_fraction = gene2profile_mrna_periph_fraction / \
                                        np.matlib.repmat(gene2profile_mrna_periph_fraction[2],
                                                         len(constants.dataset_config['MRNA_GENES']), 1)
    return gene2profile_mrna_periph_fraction


def build_histogram_peripheral_fraction(analysis_repo, molecule_type, force2D=False):
    gene2periph_fraction = {}
    gene2median_periph_fraction = {}
    gene2error = {}
    gene2ci = {}
    if molecule_type == 'mrna':
        genes = constants.analysis_config['MRNA_GENES']
    else:
        genes = constants.analysis_config['PTOTEINS']
    for gene in genes:
        image_set = ImageSet(analysis_repo, ['{0}/{1}/'.format(molecule_type, gene)], force2D=force2D)
        if molecule_type == 'mrna':
            gene2periph_fraction[gene] = image_set.compute_histogram_spots_peripheral_counts()
        else:
            gene2periph_fraction[gene] = image_set.compute_intensities_fractions_from_periphery()
        gene2median_periph_fraction[gene] = np.median(gene2periph_fraction[gene])
        gene2error[gene] = helpers.sem(gene2periph_fraction[gene], factor=0)
        lower, higher = helpers.median_confidence_interval(gene2periph_fraction[gene])
        gene2ci[gene] = [lower, higher]

    return gene2median_periph_fraction, gene2periph_fraction, gene2error, gene2ci


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
    plot.bar_profile_median(medians, errors.values(), 'mrna', xlabels, tgt_fp, confidence_interval=CI,
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
            gene2profile_mrna_periph_fraction = build_mrna_peripheral_fraction_profiles(repo)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT'].format(molecule_type="mrna")
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.profile(gene2profile_mrna_periph_fraction, constants.dataset_config['MRNA_GENES'],
                         constants.analysis_config['NUM_CONTOURS'], figname=tgt_fp)
            logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])

            logger.info("Peripheral fraction histogram for mRNA the original data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='mrna')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

        elif "chx" in conf[0]:
            logger.info("Peripheral fraction histogram for the protein CHX data")
            medians, fractions, err, CI = build_histogram_peripheral_fraction(repo, molecule_type='protein')
            plot_bar_profile_median_and_violin(molecule_type='mrna', medians=medians, fractions=fractions,
                                               errors=err, CI=CI, annotations=stat_annotations)

            logger.info("Peripheral fraction histogram for the mRNA CHX data")
            # this analysis is done as if it is protein to be in 2D
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


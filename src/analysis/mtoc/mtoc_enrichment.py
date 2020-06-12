#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import pandas as pd
from loguru import logger

import constants
import plot
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint

"""
for each timepoint, compute given function compute median, low and up envelope
normalize data by mean of median timepoint result
"""

OUTLIERS_THRESHOLD = 3
SKIP_CALL = False

def compute_protein_intensity_per_quadrant():
    dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
    protein_intensities_per_quadrant = []
    for protein in constants.dataset_config['PROTEINS']:
        for timepoint in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(analysis_repo, ["protein/{0}/{1}/".format(protein, timepoint)])
            if image_set.__sizeof__() < 5:
                logger.debug("Image set is small for {}", protein)
            dict_protein = image_set.compute_normalised_quadrant_densities(
                quadrants_num=3, mtoc_quadrant_label='MTOC',
                quadrant_labels=["Non MTOC1", "Non MTOC2"])
            dict_protein["Gene"] = [protein for i in range(len(dict_protein["MTOC"]))]
            dict_protein["MTOC leading edge"] = image_set.mtoc_is_in_leading_edge()
            dict_protein["Timepoint"] = timepoint
            protein_intensities_per_quadrant.append(pd.DataFrame(dict_protein))

    return pd.concat(protein_intensities_per_quadrant)


def compute_mrna_counts_per_quadrant():
    dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
    mrna_spot_counts_per_quadrant = []
    for gene in constants.dataset_config['MRNA_GENES']:
        for timepoint in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
            dict_gene = image_set.compute_normalised_quadrant_densities(
                quadrants_num=4, mtoc_quadrant_label='MTOC',
                quadrant_labels=["Non MTOC1", "Non MTOC2", "Non MTOC3"])
            dict_gene["Gene"] = [gene for i in range(len(dict_gene["MTOC"]))]
            dict_gene["MTOC leading edge"] = image_set.mtoc_is_in_leading_edge()
            dict_gene["Timepoint"] = timepoint
            mrna_spot_counts_per_quadrant.append(pd.DataFrame(dict_gene))

    return pd.concat(mrna_spot_counts_per_quadrant)


# To run peripheral mode analysis
#   change periph_mode to True and also
#   change the analysis config file accordingly: SPOTS_MIN : 10 => 50
#   change the analysis config file accordingly: FIGURE_OUTPUT_PATH : output/ => output_periph
# periph_mode = False

#  Figure 3C : mRNA MTOC enrichment for original data
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_original.json"))

logger.info("MTOC enrichment for the mRNA original data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

mrna_df = compute_mrna_counts_per_quadrant()  # periph_mode=periph_mode)

# Additional Figure 3C mRNA MTOC ratio MTOC/Non MTOC for original data
plot.plot_hist_ratio(mrna_df, 'mrna', limit_threshold=OUTLIERS_THRESHOLD)


# Additional Figure 3.C left mRNA MTOC enrichment for original data
plot.plot_mtoc_enrichment(mrna_df, 'mrna', quadrant_labels=["Non MTOC1", "Non MTOC2", "Non MTOC3"], limit_threshold=OUTLIERS_THRESHOLD)
# Figure 3.D left mRNA MPI for original data
plot.plot_MPI(mrna_df, 'mrna', quadrant_labels=["Non MTOC1", "Non MTOC2", "Non MTOC3"])

# Figure 3C protein MTOC enrichment for original data
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_original.json"))

logger.info("MTOC enrichment for the protein original data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

prot_df = compute_protein_intensity_per_quadrant()#periph_mode=periph_mode)

# New Figure page 3 protein MTOC ratio MTOC/Non MTOC for original data
plot.plot_hist_ratio(prot_df, 'protein', limit_threshold=OUTLIERS_THRESHOLD)
# New Figure 3.C right protein MTOC enrichment for original data
plot.plot_mtoc_enrichment(prot_df, 'protein', quadrant_labels=["Non MTOC1", "Non MTOC2"], limit_threshold=5)
# Figure 3.D right protein MPI for original data
plot.plot_MPI(prot_df, 'protein', quadrant_labels=["Non MTOC1", "Non MTOC2"])
# Figure 3.E dynamic MPI for original data
plot.plot_boxplot_MPI(mrna_df, prot_df, constants.analysis_config['PROTEINS'])

#  Analysis MTOC enrichment for nocodazole arhgdia mRNA data
constants.init_config(
    analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_nocodazole_arhgdia.json"))

logger.info("MTOC enrichment for the nocodazole arhgdia FISH data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

mrna_df = compute_mrna_counts_per_quadrant()#periph_mode=periph_mode)

# New Figure page 5 mRNA MTOC ratio MTOC/Non MTOC for nocodazole data
plot.plot_hist_ratio(mrna_df, 'mrna')
# # # New Figure 5.D left mRNA MTOC enrichment for nocodazole data
plot.plot_mtoc_enrichment(mrna_df, 'mrna', limit_threshold=OUTLIERS_THRESHOLD)
# # Figure 5.E MPI for nocodazole arhdgia mRNA data
plot.plot_MPI(mrna_df, 'mrna')

# Analysis MTOC enrichment for nocodazole arhgdia IF data
logger.info("MTOC enrichment for the nocodazole arhgdia IF data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

prot_df = compute_protein_intensity_per_quadrant()#periph_mode=periph_mode)

# Figure protein MTOC ratio MTOC/Non MTOC for nocodazole data
plot.plot_hist_ratio(prot_df, 'protein')
# New Figure 5.D left protein MTOC enrichment for nocodazole data
plot.plot_mtoc_enrichment(prot_df, 'protein', limit_threshold=OUTLIERS_THERSHOLD)
# Figure XXXXXX MPI for nocodazole arhdgia protein data
plot.plot_MPI(prot_df, 'protein')

# MTOC analysis protein MTOC enrichment for nocodazole pard3 mRNA data
constants.init_config(
    analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_nocodazole_pard3.json"))

target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))
logger.info("MTOC enrichment for the nocodazole pard3 FISH data")

mrna_df = compute_mrna_counts_per_quadrant()#periph_mode=periph_mode)

# Figure mRNA MTOC ratio MTOC/Non MTOC for nocodazole data
plot.plot_hist_ratio(mrna_df, 'mrna')

# Analysis mRNA MTOC enrichment for nocodazole arhgdia data
plot.plot_mtoc_enrichment(mrna_df, 'mrna', limit_threshold=OUTLIERS_THERSHOLD)

# Figure mRNA MPI for nocodazole arhgdia data
plot.plot_MPI(mrna_df, 'mrna')

# MTOC analysis protein MTOC enrichment for nocodazole pard3 IF data
logger.info("MTOC enrichment for the nocodazole pard3 IF data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

prot_df = compute_protein_intensity_per_quadrant()#periph_mode=periph_mode)

# Figure protein MTOC ratio MTOC/Non MTOC for nocodazole data
plot.plot_hist_ratio(prot_df, 'protein')
# Figure protein MTOC enrichment for nocodazole arhgdia data
plot.plot_mtoc_enrichment(prot_df, 'protein', limit_threshold=OUTLIERS_THRESHOLD)
# # # Figure protein MPI for nocodazole arhgdia data
plot.plot_MPI(prot_df, 'protein')

# Figure XXXXX Analysis mRNA MTOC enrichment for prrc2c data
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_prrc2c.json"))

logger.info("MTOC enrichment for the prrc2c FISH data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

mrna_df = compute_mrna_counts_per_quadrant()#periph_mode=periph_mode)

# Figure mRNA MTOC ratio MTOC/Non MTOC for prrc2c data
plot.plot_hist_ratio(mrna_df, 'mrna', groupby=['Timepoint'])
# Analysis mRNA MTOC enrichment for prrc2c data
plot.plot_mtoc_enrichment(mrna_df, 'mrna', limit_threshold=OUTLIERS_THRESHOLD)

# Figure mRNA MPI for prcc2C data
plot.plot_MPI(mrna_df, 'mrna', groupby=['Gene', 'Timepoint'])

# Figure XXXXX Analysis protein MTOC enrichment for prrc2c data
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/mtoc/config_prrc2c.json"))

logger.info("MTOC enrichment for the prrc2c IF data")
target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir))

prot_df = compute_protein_intensity_per_quadrant()# periph_mode=periph_mode)

# Figure protein MTOC ratio MTOC/Non MTOC for prrc2c data
plot.plot_hist_ratio(prot_df, 'protein', groupby=['Timepoint'])
# # # Figure protein MTOC enrichment for prrc2c data
plot.plot_mtoc_enrichment(prot_df, 'protein', limit_threshold=OUTLIERS_THRESHOLD)
# # # Figure protein MPI for prcc2C data
plot.plot_MPI(prot_df, 'protein', groupby=['Gene', 'Timepoint'])

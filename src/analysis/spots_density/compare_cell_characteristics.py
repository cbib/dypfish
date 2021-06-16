#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import numpy as np
import pandas as pd
from loguru import logger

import constants
import helpers
import plot
from helpers import open_repo
from image_set import ImageSet
from path import global_root_dir


def compute_cells_area(gene, analysis_repo, gene_label, timepoint):
    dict_Cell_area = {"Gene": [], "value": []}
    image_set = ImageSet(analysis_repo, [f"{'mrna'}/{gene}/{timepoint}/"])
    [dict_Cell_area["value"].append(image.compute_cell_area()) for image in image_set.get_images()]
    [dict_Cell_area["Gene"].append(gene_label) for image in image_set.get_images()]
    return pd.DataFrame(dict_Cell_area)


def compute_nuclei_area(gene, analysis_repo, gene_label, timepoint):
    dict_nucleus_area = {"Gene": [], "value": []}
    image_set = ImageSet(analysis_repo, [f"{'mrna'}/{gene}/{timepoint}/"])
    [dict_nucleus_area["value"].append(image.compute_nucleus_area()) for image in image_set.get_images()]
    [dict_nucleus_area["Gene"].append(gene_label) for image in image_set.get_images()]
    return pd.DataFrame(dict_nucleus_area)


def compute_transcript_by_cell_area(analysis_repo, gene, timepoints):
    transcript_by_cell_area = {"total_transcript": [], "cell_area": []}
    for timepoint in timepoints:
        image_set = ImageSet(analysis_repo, [f"{'mrna'}/{gene}/{timepoint}/"], force2D=True)
        [transcript_by_cell_area["total_transcript"].append(image.compute_cytoplasmic_total_spots()) for image in image_set.get_images()]
        [transcript_by_cell_area["cell_area"].append(image.compute_cell_area()) for image in image_set.get_images()]
    return pd.DataFrame(transcript_by_cell_area)


# configurations contain the order in which the degree of clustering is plotted
configurations = [
    ["src/analysis/spots_density/config_standard.json", ["arhgdia", "arhgdia_cultured"]]
]

# Figure 2.C bottom left (nucleus area) Micropatterned vs Standard
# Figure 2.C bottom right (cell area) Micropatterned vs Standard
# Figure 2.C Top left(transcript total count by cell area) Standard cultured cell
# Figure 2.C Top right(transcript total count by cell area) Micropatterned cultured cell

if __name__ == '__main__':

    for conf in configurations:
        stat_annotations = False
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()

        genes = constants.analysis_config['MRNA_GENES']
        genes_label = constants.analysis_config['MRNA_GENES_LABEL']
        plot_colors = constants.analysis_config['PLOT_COLORS']
        timepoints = constants.analysis_config['TIMEPOINTS']
        my_pal = {"Micropatterned cells": "#66b2ff", "Standard cells": "#1a8cff"}

        logger.info("cell area comparison for arhgdia et arhgdia_cultured 3h data")
        df_cell_area = pd.DataFrame(columns=["Gene", "value"])
        df_nucleus_area = pd.DataFrame(columns=["Gene", "value"])
        for gene, timepoint, label in zip(genes, [constants.analysis_config['TIMEPOINTS'][0][1], constants.analysis_config['TIMEPOINTS'][1][1]], genes_label):
            df_cell_area = pd.concat([df_cell_area, compute_cells_area(gene, repo, label, timepoint)])
            df_nucleus_area = pd.concat([df_nucleus_area, compute_nuclei_area(gene, repo, label, timepoint)])

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BOXPLOT'].format(model="cell_area")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.sns_boxplot(df_cell_area, my_pal, tgt_fp)

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BOXPLOT'].format(model="nucleus_area")
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
        plot.sns_boxplot(df_nucleus_area, my_pal, tgt_fp)

        logger.info("transcript total count by cell area comparison for arhgdia et arhgdia_cultured data")
        for gene, timepoints, i in zip(genes, [constants.analysis_config['TIMEPOINTS'][0], constants.analysis_config['TIMEPOINTS'][1]], [0,1]):
            transcript_by_cell_area = compute_transcript_by_cell_area(repo, gene, timepoints)
            if gene == 'arhgdia':
                outliers = helpers.detect_outliers(transcript_by_cell_area['total_transcript'], threshold=1.8)
                transcript_by_cell_area = transcript_by_cell_area[~np.isin(transcript_by_cell_area["total_transcript"], outliers)]
                transcript_by_cell_area.dropna(inplace=True)
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT'].format(cell_type=genes[i])
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
            plot.sns_linear_regression(transcript_by_cell_area["cell_area"], transcript_by_cell_area["total_transcript"], plot_colors[i], tgt_fp, order=2)
            logger.info("Generated image at {}", tgt_fp)

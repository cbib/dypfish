#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from scipy import interpolate

import src.acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, plot_colors, cell_type_micropatterned, check_dir


def main():
    enable_logger()

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts

    with h5py.File(path.basic_file_path, "r") as file_handler:

        molecule_type = ['/mrna']

        #genes = ["rack1_control", "rack1_sirna"]
        plot_colors = ['#1E95bb', '#1ec5d4']
        genes = [ "arhgdia_control","arhgdia_cytod"]
        #plot_colors = ['#F16c1b', '#f1bc1b']
        #genes = ["pard3_control", "pard3_cytod"]
        #genes = ["rack1_sirna", "rack1", "arhgdia_cytod", "arhgdia_control"]
        #genes = ["rack1_sirna", "rack1_control", "arhgdia_cytod", "arhgdia_control", "pard3_cytod", "pard3_control"]

        cytoplasmic_transcripts = []
        for gene in genes:
            print(gene)
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data))

        plot_filename = molecule_type[0] + '_total_cytoplasmic_transcript_arhgdia_cytod.png'
        #plot_filename = molecule_type[0] + '_total_cytoplasmic_transcript_pard3_cytod.png'
        #plot_filename = molecule_type[0] + '_total_cytoplasmic_transcript_rack1_sirna.png'
        figname = check_dir(path.analysis_dir + 'analysis_cytoD/figures/') + plot_filename
        stan.plot_bar_profile(cytoplasmic_transcripts, genes, 30, "Cytoplasmic total transcript", figname,
                              plot_colors)




if __name__ == "__main__":
    main()

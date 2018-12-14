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

        # molecule_type = ['/protein']
        # genes = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        #
        # cytoplasmic_intensities = []
        # for gene in genes:
        #     image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
        #     cytoplasmic_intensities.append(adsc.compute_cytoplasmic_total(file_handler, image_list, path.path_data))
        # figname = path.analysis_dir + 'analysis_cytoplasmic_total_count/figures/' + molecule_type[
        #     0] + '_total_cytoplasmic_intensity.png'
        # stan.plot_bar_profile(cytoplasmic_intensities, genes, 1000000000,
        #                       "Cytoplasmic total intensity (arbitrary units)", figname, plot_colors)

    # # This part produces plot interpolation of cytoplasmic spread by timepoint
    # with h5py.File(path.basic_file_path, "r") as file_handler:
    #
    #     molecule_type = ['/mrna']
    #     genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    #     proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    #
    #     for i in range(len(genes)):
    #
    #         fig = plt.figure()
    #         ax = plt.axes()
    #         ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    #         ax.set_xlim(1, 8)
    #         ax.set_ylim(0, 3.50)
    #
    #         x_mrna = np.arange(2, 5, 0.01)
    #         mrna_data = np.zeros((3, 4))
    #         counter = 0
    #         timepoints = ["2h", "3h", "4h", "5h"]
    #         for timepoint in timepoints:
    #             print(genes[i], '_', timepoint)
    #             image_list = helps.preprocess_image_list3(file_handler, molecule_type, genes[i], [timepoint])
    #             cyt_spread = adsc.compute_cytoplasmic_total(file_handler, image_list, path.path_data)
    #             cyt_spread_median = np.median(cyt_spread)
    #             print(cyt_spread_median)
    #             err = np.median(np.abs(np.tile(np.median(cyt_spread), (1, len(cyt_spread))) - cyt_spread))
    #             upp_env = cyt_spread_median + err
    #             low_env = cyt_spread_median - err
    #             mrna_data[0, counter] = cyt_spread_median
    #             mrna_data[1, counter] = upp_env
    #             mrna_data[2, counter] = low_env
    #             counter += 1
    #         mrna_data = mrna_data / np.mean(mrna_data[0, :])
    #
    #         spl = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[0, :])
    #         spl_upp = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[1, :])
    #         spl_low = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[2, :])
    #         m_y_new = spl(x_mrna)
    #         m_y_new_upp = spl_upp(x_mrna)
    #         m_y_new_down = spl_low(x_mrna)
    #
    #         solid_mrna, = plt.plot(x_mrna, m_y_new, linestyle="-", color=plot_colors[i], linewidth=2)
    #         plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=plot_colors[i])
    #         plt.plot(x_mrna, m_y_new_down, linestyle="-", color=plot_colors[i])
    #         ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=plot_colors[i], alpha=0.5, interpolate=False)
    #
    #         if genes[i] in proteins:
    #
    #             counter = 0
    #             x_protein = np.arange(2, 7, 0.01)
    #             protein_data = np.zeros((3, 4))
    #             timepoints = ["2h", "3h", "5h", "7h"]
    #             for timepoint in timepoints:
    #                 print(genes[i], '_', timepoint)
    #                 image_list = helps.preprocess_image_list3(file_handler, ['/protein'], genes[i], [timepoint])
    #                 cyt_spread = adsc.compute_cytoplasmic_total(file_handler, image_list, path.path_data)
    #                 cyt_spread_median = np.median(cyt_spread)
    #                 print(cyt_spread_median)
    #                 err = np.median(np.abs(np.tile(np.median(cyt_spread), (1, len(cyt_spread))) - cyt_spread))
    #                 upp_env = cyt_spread_median + err
    #                 low_env = cyt_spread_median - err
    #                 protein_data[0, counter] = cyt_spread_median
    #                 protein_data[1, counter] = upp_env
    #                 protein_data[2, counter] = low_env
    #                 counter += 1
    #             protein_data = protein_data / np.mean(protein_data[0, :])
    #
    #             spl = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[0, :])
    #             spl_upp = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[1, :])
    #             spl_low = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[2, :])
    #             p_y_new = spl(x_protein)
    #             p_y_new_upp = spl_upp(x_protein)
    #             p_y_new_down = spl_low(x_protein)
    #
    #             dashed_protein, = plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=plot_colors[i],
    #                                        linewidth=2)
    #             plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=plot_colors[i])
    #             plt.plot(x_protein, p_y_new_down, linestyle="--", color=plot_colors[i])
    #             ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=plot_colors[i], alpha=0.25,
    #                             interpolate=False)
    #             plt.legend([dashed_protein, solid_mrna], ['Protein', 'Mrna'])
    #         else:
    #             plt.legend([solid_mrna], ['Mrna'])
    #
    #         ax.set_ylabel('Cytoplasmic total')
    #         ax.set_xlabel('Time(hrs)')
    #         ax.set_title(genes[i])
    #         plt.savefig(path.analysis_dir + 'analysis_cytoplasmic_total_count/figures/cyt_total_' + genes[i] + ".png",
    #                     format='png')
    #         plt.show()


if __name__ == "__main__":
    main()

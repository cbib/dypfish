#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues


import h5py
import src.plot as plot
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors, check_dir


def cytoplasmic_total_dynamic_profiles(file_handler, genes, proteins):
    data_generator = plot.data_extractor(genes, proteins, file_handler,
                                         adsc.compute_cytoplasmic_total, file_handler, path.path_data)
    for mrna_data, protein_data, i in data_generator:
        figpath = check_dir(path.analysis_dir + 'analysis_cytoplasmic_total_count/figures/') + 'cyt_total_' + genes[
            i] + ".png"
        plot.dynamic_profiles(mrna_data, protein_data, genes[i], plot_colors[i], 'Time(hrs)', 'Cytoplasmic total',
                              figpath)


def cytoplasmic_total_count(file_handler, molecule_type, genes):
    cytoplasmic_total = []
    for gene in genes:
        image_list = helps.preprocess_image_list2(file_handler, molecule_type, gene)
        cytoplasmic_total.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data))
    figname = check_dir(
        path.analysis_dir + 'analysis_cytoplasmic_total_count/figures/') + molecule_type + '_total_cytoplasmic_transcript.png'
    plot.bar_profile(cytoplasmic_total, genes, figname)


def main():
    enable_logger()

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts

    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(path.basic_file_path, "r") as file_handler:
        cytoplasmic_total_count(file_handler, '/mrna', genes)
        cytoplasmic_total_count(file_handler, '/protein', proteins)
        cytoplasmic_total_dynamic_profiles(file_handler, genes, proteins)


if __name__ == "__main__":
    main()

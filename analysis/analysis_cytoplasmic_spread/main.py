#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.plot as plot
from src import acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors, cell_type_micropatterned, check_dir


def cytoplasmic_spread(file_handler, genes, path_data, cell_type, molecule_type):  # TODO argument genes vs proteins ?
    # for each gene in genes list, compute cytoplasmic spread
    cyt_spreads = []
    plot_filename = check_dir(path.analysis_dir + 'analysis_cytoplasmic_spread/figures/' + cell_type + '/') + \
                    molecule_type + '_cytoplasmic_spread.png'
    for gene in genes:
        image_list = helps.preprocess_image_list2(file_handler, '/' + molecule_type, gene)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread(image_list, file_handler, path_data))
    plot.bar_profile(cyt_spreads, genes, 1.60, 'Cytoplasmic spread', plot_filename, plot_colors)


def cytoplasmic_spread_dynamic_profiles(file_handler, genes, proteins):
    # This part produces plot interpolation of cytoplasmic spread by timepoint
    data_generator = plot.data_extractor(genes, proteins, file_handler,
                                         adsc.compute_cytoplasmic_spread, file_handler, path.path_data)
    for mrna_data, protein_data, i in data_generator:
        figpath = check_dir(path.analysis_dir + '/analysis_cytoplasmic_spread/figures/') + 'cyt_spread_' + genes[
            i] + '.png'
        plot.dynamic_profiles(mrna_data, protein_data, genes[i], plot_colors[i], 'Time(hrs)', 'Cytoplasmic spread',
                              figpath)

'''compute cytoplasmic spread histogram and dynamic profiles'''
def main():
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()

    # Compute bar plot cytoplasmic spread
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(path.basic_file_path, "r") as file_handler:
        cytoplasmic_spread(file_handler, genes, path.path_data, cell_type_micropatterned, molecule_type='mrna')
        cytoplasmic_spread(file_handler, proteins, path.path_data, cell_type_micropatterned, molecule_type='protein')
        cytoplasmic_spread_dynamic_profiles(file_handler, genes, proteins)




if __name__ == "__main__":
    main()

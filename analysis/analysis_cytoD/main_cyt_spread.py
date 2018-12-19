#!/usr/bin/python
# encoding: UTF-8

import h5py
from src import acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, check_dir

def cytoplasmic_spread(file_handler, genes, path_data, plot_colors, molecule_type):
    cyt_spreads = []
    plot_filename = check_dir(path.analysis_dir + 'analysis_cytoD/figures/') + \
                    molecule_type + '_cytoplasmic_spread_arhgdia_ctrl_vs_cytod.png'
    for gene in genes:
        print(gene)
        image_list = helps.preprocess_image_list2(file_handler, '/' + molecule_type, gene)
        print(image_list)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread_2D(file_handler, image_list, path_data))
    stan.plot_bar_profile(cyt_spreads, genes, 1.60, 'Cytoplasmic spread',  plot_filename, plot_colors)

def main():
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()
    plot_colors = ['#1E95bb', '#1ec5d4']
    genes = ["arhgdia_control", "arhgdia_cytod"]
    cell_type = cell_type_micropatterned
    with h5py.File(path.basic_file_path, "r") as file_handler:
        cytoplasmic_spread(file_handler, genes, path.path_data, cell_type,plot_colors, molecule_type='mrna')

if __name__ == "__main__":
    main()

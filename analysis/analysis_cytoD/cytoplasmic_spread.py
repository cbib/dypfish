#!/usr/bin/python
# encoding: UTF-8

import h5py
import argparse
import tqdm

from src import acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, check_dir, plot_colors_cytoD, loadconfig
import src.plot as plot

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


def cytoplasmic_spread(file_handler, genes, path_data, plot_colors, molecule_type):
    cyt_spreads = []

    plot_filename = check_dir(path.analysis_dir + 'analysis_cytoD/figures/cyt_spread/') + molecule_type + '_cytoplasmic_spread_arhgdia_ctrl_vs_cytod.png'
    for gene in tqdm.tqdm(genes, desc="Processing mRNAs"):
        image_list = helps.preprocess_image_list2(file_handler, '/' + molecule_type, gene)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread_2D(file_handler, image_list, path_data))
    plot.bar_profile(cyt_spreads, genes, plot_filename)

def main():
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()
    configData = loadconfig(input_dir_name)
    genes = configData["GENES"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    check_dir(path.analysis_dir + 'analysis_cytoD/figures/')

    with h5py.File(path.data_dir+input_dir_name+"/"+basic_file_name, "r") as file_handler:
        cytoplasmic_spread(file_handler, genes, path.path_data, plot_colors_cytoD, molecule_type='mrna')
        cytoplasmic_spread(file_handler, genes, path.path_data, plot_colors_cytoD, molecule_type='protein')

if __name__ == "__main__":
    main()
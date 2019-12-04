#!/usr/bin/python
# encoding: UTF-8

import h5py
import argparse
import src.acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, check_dir, loadconfig, plot_colors_cytoD
import src.plot as plot

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def main():
    enable_logger()

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts

    configData = loadconfig(input_dir_name)
    genes = configData["GENES"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    check_dir(path.analysis_dir + 'analysis_cytoD/figures/')



    with h5py.File(path.data_dir+input_dir_name+"/"+basic_file_name, "r") as file_handler:
        molecule_type = ['/mrna']
        cytoplasmic_transcripts = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data))
        plot_filename = molecule_type[0] + '_total_cytoplasmic_transcript_arhgdia_cytod.png'
        figname = check_dir(path.analysis_dir + 'analysis_cytoD/figures/cyt_total/') + plot_filename

        plot.bar_profile(cytoplasmic_transcripts, genes, figname)
if __name__ == "__main__":
    main()

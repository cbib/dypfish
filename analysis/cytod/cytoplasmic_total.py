#!/usr/bin/python
# encoding: UTF-8

import h5py
import argparse
import tqdm
import src.acquisition_descriptors as adsc
import src.path as path
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
    molecule_types = configData["MOLECULE_TYPES"]
    colors=configData["COLORS"]
    check_dir(path.analysis_dir + 'cytod/figures/')


    with h5py.File(path.data_dir+input_dir_name+"/"+basic_file_name, "r") as file_handler:
        for molecule_type in molecule_types:
            cytoplasmic_transcripts = []
            for gene in tqdm.tqdm(genes,desc="Processing mRNAs"):
                image_list = helps.preprocess_image_list2(file_handler, molecule_type, gene)
                cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler))
            plot_filename = molecule_type + '_total_cytoplasmic_transcript_arhgdia_cytod.png'
            figname = check_dir(path.analysis_dir + 'cytod/figures/cyt_total/') + plot_filename
            plot.bar_profile(cytoplasmic_transcripts, genes, figname, colors)
if __name__ == "__main__":
    main()

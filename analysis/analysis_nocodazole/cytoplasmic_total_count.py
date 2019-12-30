#!/usr/bin/python
# encoding: UTF-8

import h5py
import logging
import argparse
import sys
import src.acquisition_descriptors as adsc
from src.path import path_data, analysis_dir
import src.path as path
import src.helpers as helps
from src.utils import enable_logger
from src.utils import check_dir, loadconfig
import src.plot as plot

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


if __name__ == "__main__":
    # Required descriptors: spots, IF, cell mask an height_map

    enable_logger()
    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    colors = configData["COLORS"]
    molecule_types = configData["MOLECULE_TYPES"]


    # Compute bar plot cytoplasmic total transcripts
    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name,  "r") as file_handler:
        for molecule_type in molecule_types:
            cytoplasmic_total = []
            for gene in mrnas:
                image_list = helps.preprocess_image_list2(file_handler, molecule_type, gene)
                cytoplasmic_total.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path_data))
            figname= check_dir(analysis_dir + 'analysis_nocodazole/figures/cytoplasmic_total/') + molecule_type +'_total_cytoplasmic_transcript.png'
            plot.bar_profile(cytoplasmic_total, mrnas, figname)





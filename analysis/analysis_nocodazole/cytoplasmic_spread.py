#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import argparse
import h5py
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import check_dir,loadconfig
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

def cytoplasmic_spread(file_handler, molecule_type,genes,image_width,image_height):
    cyt_spreads = []
    figname = check_dir(path.analysis_dir + '/analysis_nocodazole/figures/cytoplasmic_spread/') + molecule_type + '_cytoplasmic_spread.png'
    for gene in genes:
        image_list = helps.preprocess_image_list2(file_handler, molecule_type, gene)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread_2D(image_list, file_handler, image_width,image_height))
    plot.bar_profile(cyt_spreads, genes, figname)

if __name__ == "__main__":
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map

    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"]

    basic_file_name = configData["BASIC_FILE_NAME"]
    colors = configData["COLORS"]
    molecule_types = configData["MOLECULE_TYPES"]
    image_width = configData["IMAGE_WIDTH"]
    image_height = configData["IMAGE_HEIGHT"]


    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name, "r") as file_handler:
        for molecule_type in molecule_types:
            cytoplasmic_spread(file_handler,molecule_type,mrnas,image_width,image_height)

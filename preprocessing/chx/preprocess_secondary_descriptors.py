#!/usr/bin/python
# encoding: UTF-8

import h5py
import argparse
import src.path as path
import src.helpers as helps
import src.image_descriptors as idsc
from src.utils import enable_logger, loadconfig


''' 
2-This script is supposed to be ran to compute cell_mask_distance_map and spot_peripheral_distance descriptors in dypfish. It produces:
'''

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def main():

    enable_logger()

    configData = loadconfig(input_dir_name)
    max_cell_radius = configData["MAX_CELL_RADIUS"]
    simulation_number = configData["RIPLEY_K_SIMULATION_NUMBER"]
    contours_num=configData["NUM_CONTOURS"]
    image_width=configData["IMAGE_WIDTH"]
    image_height = configData["IMAGE_HEIGHT"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    hstar_file_name = configData["HSTAR_FILE_NAME"]
    size_coeff = configData["SIZE_COEFFICIENT"]


    with h5py.File(path.data_dir+input_dir_name+"/"+basic_file_name, "a") as input_file_handler, \
            h5py.File(path.data_dir+input_dir_name+"/"+secondary_file_name, "a") as secondary_file_handler, \
            h5py.File(path.data_dir+input_dir_name+"/"+hstar_file_name, "a") as hstar_file_handler:

        image_list = helps.preprocess_image(input_file_handler)
        for image in image_list:
            idsc.set_cell_mask_distance_map(input_file_handler, secondary_file_handler, image, contours_num,image_width, image_height, max_cell_radius)
            idsc.set_cell_area(input_file_handler, secondary_file_handler, image, size_coeff)
            idsc.set_nucleus_area(input_file_handler, secondary_file_handler, image, size_coeff)
            idsc.set_h_star_protein_2D(input_file_handler, hstar_file_handler, image,max_cell_radius, simulation_number)

if __name__ == "__main__":
    main()
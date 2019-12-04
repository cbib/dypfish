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
    # Required descriptors: cell_mask, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()
    configData = loadconfig("cytoD")
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    contours_num = configData["NUM_CONTOURS"]
    with h5py.File(path.basic_file_path, "a") as input_file_handler, \
            h5py.File(path.secondary_file_path, "a") as secondary_file_handler:

        image_list = helps.preprocess_image_list_1(input_file_handler, genes)
        #image_list = helps.preprocess_image(input_file_handler)
        for image in image_list:
            print(image)
            idsc.set_cell_mask_distance_map(input_file_handler, secondary_file_handler, image, contours_num)
            if 'mrna' in image:
                 idsc.set_spots_peripheral_distance_2D(input_file_handler, secondary_file_handler, image)
            idsc.set_cell_area(input_file_handler, secondary_file_handler, image)
            idsc.set_nucleus_area(input_file_handler, secondary_file_handler, image)

if __name__ == "__main__":
    main()
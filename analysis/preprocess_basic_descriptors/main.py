#!/usr/bin/python
# encoding: UTF-8

import matplotlib.pyplot as plt
import h5py
import math
import numpy as np
import sys
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
import src.plot as plot
import src.image_descriptors as idsc

from scipy import interpolate
from src.utils import enable_logger, plot_colors, check_dir

''' 
2-This script is supposed to be ran to compute cell_mask_distance_map and spot_peripheral_distance descriptors in dypfish. It produces:
'''




def main():
    # Required descriptors: cell_mask, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    enable_logger()

    with h5py.File(path.basic_file_path, "a") as input_file_handler, h5py.File(path.secondary_file_path, "a") as secondary_file_handler:
        image_list = helps.preprocess_image(input_file_handler)
        for image in image_list:
            print(image)
            idsc.set_cell_mask_distance_map(input_file_handler, secondary_file_handler, image)
            if 'mrna' in image:
                idsc.set_spots_peripheral_distance(input_file_handler, secondary_file_handler, image)
                idsc.set_spots_peripheral_distance_2D(input_file_handler, secondary_file_handler, image)

            idsc.set_cell_area(input_file_handler, secondary_file_handler, image)
            idsc.set_nucleus_area(input_file_handler, secondary_file_handler, image)


if __name__ == "__main__":
    main()




#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.path as path
import src.helpers as helps
import src.acquisition_descriptors as adsc
from src.utils import enable_logger, check_dir
from numpy import mean, median, absolute

def mean_absolute_deviation(data, axis=None):
    return mean(absolute(data - mean(data, axis)), axis)

def median_absolute_deviation(data, axis=None):
    return median(absolute(data - median(data, axis)), axis)

def get_spots_peripheral_distance_2D(spots,nucleus_mask,cell_mask,periph_dist_map):
    spots_peripheral_distance_2D = []
    for spot in spots:
        if nucleus_mask[spot[1], spot[0]] == 0:
            if cell_mask[spot[1], spot[0]] == 1 and periph_dist_map[spot[1], spot[0]] == 0:
                spots_peripheral_distance_2D.append(int(1))
            else:
                spots_peripheral_distance_2D.append(periph_dist_map[spot[1], spot[0]])
    print(spots_peripheral_distance_2D)
    return spots_peripheral_distance_2D



def main():
    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()

    ## Control reproducibility mrna distribution for arhgdia and arhgdia cultured
    with h5py.File(path.basic_file_path, "a") as file_handler, \
            h5py.File(path.secondary_file_path, "a") as sec_file_handler:
        check_dir(path.analysis_dir + "/analysis_stability/figures/")
        arhgdia = helps.build_image_list_2(file_handler, 'mrna', "arhgdia",["3h"])
        cyt_spreads=[]
        cyt_spreads.append(adsc.compute_cytoplasmic_spread_2D(file_handler, arhgdia, path.path_data))
        print(cyt_spreads)

if __name__ == "__main__":
    main()

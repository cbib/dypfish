#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import pandas as pd
from random import *
import h5py
import sys
from skimage import draw
from skimage import measure
from skimage import io
import matplotlib.pyplot as plt
from statsmodels import robust
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.image_descriptors as idsc
import src.acquisition_descriptors as adsc


from src.utils import enable_logger, check_dir

from numpy import mean, median, absolute

def mean_absolute_deviation(data, axis=None):

    return mean(absolute(data - mean(data, axis)), axis)

def median_absolute_deviation(data, axis=None):

    return median(absolute(data - median(data, axis)), axis)




def get_spots_peripheral_distance_2D(spots,nucleus_mask,cell_mask,periph_dist_map):

    # plt.imshow(periph_dist_map)
    # plt.show()

    spots_peripheral_distance_2D = []
    for spot in spots:

        if nucleus_mask[spot[1], spot[0]] == 0:
            if cell_mask[spot[1], spot[0]] == 1 and periph_dist_map[spot[1], spot[0]] == 0:
                spots_peripheral_distance_2D.append(int(1))
            else:
                spots_peripheral_distance_2D.append(periph_dist_map[spot[1], spot[0]])
    print(spots_peripheral_distance_2D)
    # xs = spots[:, 0]
    # ys = spots[:, 1]
    # plt.imshow(periph_dist_map, cmap='hot')
    # contours = measure.find_contours(nucleus_mask, 0.8)
    #
    # for n, contour in enumerate(contours):
    #     plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    # plt.scatter(xs, ys, color='blue', marker="o", facecolors='none', linewidths=0.5)
    # plt.show()
    return spots_peripheral_distance_2D



def main():
    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()

    try:
        arhgdia_df = pd.read_csv(path.analysis_dir + "analysis_stability/dataframe/arhgdia_3h_mrna.csv",index_col=0,dtype=np.float64)
    except IOError:
        print "Couldn't load file : ", path.analysis_dir + 'analysis_stability/dataframe/global_mtoc_file_all_mrna.csv'
        print "Maybe MTOC analysis hasn't been launched prior to this one"
        exit(1)
    try:
        arhgdia_cultured_df = pd.read_csv(path.analysis_dir + "analysis_stability/dataframe/arhgdia_cultured_3h_mrna.csv",index_col=0)
    except IOError:
        print "Couldn't load file : ", path.analysis_dir + 'analysis_stability/dataframe/arhgdia_3h_mrna.csv'
        print "Maybe MTOC analysis hasn't been launched prior to this one"
        exit(1)

    ## Control reproducibility mrna distribution for arhgdia and arhgdia cultured
    with h5py.File(path.basic_file_path, "a") as file_handler, \
            h5py.File(path.secondary_file_path, "a") as sec_file_handler:

        check_dir(path.analysis_dir + "/analysis_stability/figures/")

        # Build all image for an acquisition
        #arhgdia = helps.build_image_list(file_handler, 'mrna', "arhgdia")

        arhgdia = helps.build_image_list_2(file_handler, 'mrna', "arhgdia",["3h"])
        arhgdia_cultured = helps.build_image_list_2(file_handler, 'mrna', "arhgdia_cultured",["3h"])
        #arhgdia_noc = helps.build_image_list_2(file_handler, 'mrna', "arhgdia_nocodazole",["3h","5h"])
        cyt_spreads=[]
        cyt_spreads.append(adsc.compute_cytoplasmic_spread_2D(file_handler, arhgdia, path.path_data))
        print(cyt_spreads)






if __name__ == "__main__":
    main()

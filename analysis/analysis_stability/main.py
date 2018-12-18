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
from sklearn.metrics import median_absolute_error

import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.image_descriptors as idsc


from src.utils import enable_logger, check_dir

from numpy import mean, absolute

def mad(data, axis=None):

    return mean(absolute(data - mean(data, axis)), axis)




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

    total_mads_arhgdia = []
    total_mads_arhgdia_cultured=[]
    for j in range(500):
        print(j)
        mads = []


        for i in range(1,len(arhgdia_df.index.values)+1):

            arr=np.array(arhgdia_df.ix[:,:])
            arr = arr[np.random.choice(arr.shape[0], i, replace=True)]
            mean_arr=np.mean(arr,axis=0)
            rand_arr=arr[randint(0, len(arr)-1)]
            mse = mean_squared_error(mean_arr, rand_arr)

            if (i > 1):
                mads.append(mse)

        total_mads_arhgdia.append(mads)



        mads = []
        for i in range(1,len(arhgdia_cultured_df.index.values)+1):
            arr = np.array(arhgdia_cultured_df.ix[:,:])
            arr = arr[np.random.choice(arr.shape[0], i, replace=True)]

            mean_arr=np.mean(arr,axis=0)
            rand_arr=arr[randint(0, len(arr)-1)]

            mse = mean_squared_error(mean_arr, rand_arr)
            if (i > 1):
                mads.append(mse)
        total_mads_arhgdia_cultured.append(mads)

    plt.figure()
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ## the data
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    plt.plot(np.mean(total_mads_arhgdia, axis=0), color='blue')
    plt.plot(np.mean(total_mads_arhgdia_cultured, axis=0),color='black')
    plt.savefig(path.analysis_dir + "/analysis_stability/figures/stability_analysis_mtoc_max_quadrant_median_mean_squared_error.png")

    plt.show()



if __name__ == "__main__":
    main()

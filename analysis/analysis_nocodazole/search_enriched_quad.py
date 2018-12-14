#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import math
import pandas as pd
pd.set_option('display.max_rows', 500)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import interpolate
from scipy import stats
from skimage.draw import circle
from skimage import measure

import os
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.constants
from src.utils import check_dir



logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
np.set_printoptions(precision=4)
# if "log" not in globals():
# logger = Logger.init_logger('REFORMAT_%s'%(cfg.language_code), load_config())
logger.info("Running %s", sys.argv[0])


def box_plot_concentrations(k1, k2, title, figname):
    #k1 = sum(k1[:, 0], [])
    #print(k1)
    #k2 = sum(k2[:, 0], [])
    #k3 = sum(k3[:, 0], [])
    #k1_k2med = k1 - np.median(k2)
    fig, ax = plt.subplots(1, 1)
    plt.boxplot([np.matrix(k1), np.matrix(k2)], 'gD')
    plt.ylabel('Mean mrna concentration (total mrna spots) by quadrant')
    plt.title(title)
    y_lim = np.max(k1)
    plt.axis([0, 2.5, 0, y_lim + y_lim / 10])
    plt.axis([0, 2.5, -1, y_lim])

    MTOC_patch = mpatches.Patch(label='1 - Quadrants with MTOC')
    noMTOC_patch = mpatches.Patch(label='2 - Other quadrants')
    plt.legend(handles=[MTOC_patch, noMTOC_patch])
    ax.yaxis.grid(color='gray', linestyle='dashed')
    plt.savefig(figname)
    #plt.show()
    plt.close()

# function for setting the colors of the box plots pairs
def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')


def plot_bar_profile(data, genes, y_limit, ylabel, figname, colors):
    ## third technic
    fig = plt.figure()
    # ax = fig.add_subplot(111)
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ## the data
    N = len(genes)

    dataMedians = []
    dataStd = []
    y_lim_max = np.max(data) + 0.2
    y_lim_min = np.min(data) - 0.2
    print(y_lim_max)


    ## necessary variables
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars
    # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
    # colors =['#0A3950','#1E95BB','#A1BA6D','#F16C1B','#C02A18','#E9CB45']

    ## the bars
    rects1 = ax.bar(ind, data, width,
                    color=colors)

    # axes and labels
    ax.set_xlim(-width, len(ind) + width)

    ax.set_ylim(y_lim_min, y_lim_max)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    xTickMarks = ["" for i in range(0, N)]
    ax.set_xticks(ind)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.legend([gene for gene in genes], loc='upper right')
    ax.legend(rects1, genes, prop={'size': 8})
    plt.savefig(figname, format='svg')
    plt.show()


def plot_bar_profile_2_series(data,data_noc, genes,noc_genes, y_limit, ylabel, figname, colors):
    ## third technic
    fig = plt.figure()
    # ax = fig.add_subplot(111)
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ## the data
    N = len(genes)

    dataMedians = []
    dataStd = []
    y_lim_max = np.max(data) + 0.2
    y_lim_min = np.min(data) - 0.4
    print(y_lim_max)


    ## necessary variables
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars
    # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
    # colors =['#0A3950','#1E95BB','#A1BA6D','#F16C1B','#C02A18','#E9CB45']

    ## the bars
    rects1 = ax.bar(ind, data, width,
                    color=['grey','grey'])
    rects2 = ax.bar(ind+width, data_noc, width,
                    color=['#1E95BB','#F16C1B'])

    # axes and labels
    ax.set_xlim(-width, len(ind) + width)

    ax.set_ylim(y_lim_min, y_lim_max)
    ax.set_ylabel(ylabel)
    ax.set_title('')

    ax.set_xticks(ind + width)
    ax.set_xticklabels(('arhgdia', 'pard3'))
    ax.legend((rects1[0], rects2[0]), ('control', 'Nocodazole'))
    #plt.legend([gene for gene in genes], loc='upper right')
    #ax.legend(rects1, genes, prop={'size': 8})
    #ax.legend(rects2, noc_genes, prop={'size': 8})
    plt.savefig(figname, format='svg')
    plt.show()

if __name__ == "__main__":



    # Required descriptors: spots, IF, cell mask an height_map

    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    ''' 
    1-You need to create a password.txt file before running to connect via ssh
    '''
    basic_file_path = path.analysis_data_dir + 'basic.h5'
    secondary_file_path = path.analysis_data_dir + 'secondary.h5'
    mtoc_file_path = path.analysis_data_dir + 'mtoc.h5'



    with h5py.File(basic_file_path, "r") as file_handler,h5py.File(secondary_file_path, "r") as second_file_handler,h5py.File(mtoc_file_path, "r") as mtoc_file_handler:
        from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes
        molecule_type=['/mrna']
        colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B']
        mrnas = ["arhgdia", "arhgdia_nocodazole","pard3","pard3_nocodazole"]

        global_mean_mtoc=[]
        global_mean_mtoc_leading = []
        global_mean_non_mtoc=[]

        global_max_mtoc = []
        global_max_mtoc_leading = []
        global_max_non_mtoc = []

        global_mtoc=[]
        global_nmtoc=[]
        global_mtoc_leading = []

        global_mrna=[]
        global_image=[]
        global_index = []
        global_timepoint=[]
        global_mpis=[]
        global_p=[]
        global_degree=[]
        for mrna in mrnas:
            print(mrna)
            spots_mtoc_all = []
            spots_non_mtoc_all = []

            timepoints = ["3h", "5h"]

            for timepoint in timepoints:
                mean_spots_mtoc = []
                mean_spots_non_mtoc = []
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, mrna, [timepoint])


                for image in image_list:

                    #spot_by_quad= idsc.search_periph_mrna_quadrants(file_handler, second_file_handler, image)
                    spot_by_quad = idsc.search_mrna_quadrants(file_handler, image)
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    mtoc_spot = spot_by_quad[:, :, 1] == 1
                    #spots_mtoc_all.append(spot_by_quad[mtoc_spot][:,0])
                    non_mtoc_spot = spot_by_quad[:, :, 1] == 0

                    for i in range(90):

                        global_index.append(image.split("/")[4]+"_"+str(i+1))
                        global_image.append(image.split("/")[4])

                        global_mrna.append(mrna)
                        global_timepoint.append(timepoint)

                    global_mtoc.extend(spot_by_quad[mtoc_spot][:,0].flatten())
                    for i in range(0,270,3):
                        global_nmtoc.append(np.mean(spot_by_quad[non_mtoc_spot][:,0].flatten()[i:i+2]))

                    if mtoc_quad_j == 1:
                        global_mtoc_leading.extend(spot_by_quad[mtoc_spot][:,0].flatten())
                    else:
                        for i in range(90):
                            global_mtoc_leading.append(np.nan)


        df = pd.DataFrame({'Image': global_image,'Gene': global_mrna, 'timepoint': global_timepoint, 'Non MTOC': global_nmtoc,'MTOC': global_mtoc,
                            'MTOC leading edge': global_mtoc_leading},
                          index=global_index)
        df.to_csv(path.analysis_dir + 'analysis_nocodazole/df/' +'global_mtoc_file_mrna_all.csv')


        # #protein part
        molecule_type = ['/protein']
        # #proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        proteins = ["arhgdia", "arhgdia_nocodazole","pard3","pard3_nocodazole"]

        mean_intensities_mtoc = []
        mean_intensities_non_mtoc = []
        global_mean_mtoc = []
        global_mean_non_mtoc = []
        global_protein = []
        global_timepoint = []
        global_mean_mtoc = []
        global_mean_mtoc_leading = []
        global_mean_non_mtoc = []

        global_max_mtoc = []
        global_max_mtoc_leading = []
        global_max_non_mtoc = []

        global_mtoc = []
        global_nmtoc = []
        global_mtoc_leading = []

        global_image = []
        global_index = []
        global_timepoint = []
        global_mpis = []
        global_p = []
        global_degree = []
        for protein in proteins:


            #timepoints = ["2h", "3h", "5h", "7h"]
            timepoints = ["3h", "5h"]

            for timepoint in timepoints:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, protein, [timepoint])

                for image in image_list:

                    intensity_by_quad = idsc.search_protein_quadrants(file_handler, mtoc_file_handler,protein, image)
                    mtoc_intensity = intensity_by_quad[:, :, 1] == 1

                    non_mtoc_intensity = intensity_by_quad[:, :, 1] == 0


                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)

                    # needed for final boxplot
                    # global_image.append(image.split("/")[4])
                    # global_mrna.append(mrna)
                    # global_timepoint.append(timepoint)
                    #
                    # # global_mean_mtoc.append(np.mean(spot_by_quad[mtoc_spot][:,0]))
                    # # if mtoc_quad_j==1:
                    # #     global_mean_mtoc_leading.append(np.mean(spot_by_quad[mtoc_spot][:,0]))
                    # # else:
                    # #     global_mean_mtoc_leading.append(np.nan)
                    # # global_mean_non_mtoc.append(np.mean(spot_by_quad[non_mtoc_spot][:,0]))
                    #
                    # global_max_mtoc.append(np.max(spot_by_quad[mtoc_spot][:, 0]))
                    # if mtoc_quad_j == 1:
                    #     global_max_mtoc_leading.append(np.max(spot_by_quad[mtoc_spot][:, 0]))
                    # else:
                    #     global_max_mtoc_leading.append(np.nan)
                    # global_mean_non_mtoc.append(np.mean(spot_by_quad[non_mtoc_spot][:, 0]))

                    for i in range(90):
                        # needed for final boxplot
                        # global_degree.append(i)
                        global_index.append(image.split("/")[4] + "_" + str(i + 1))
                        global_image.append(image.split("/")[4])

                        global_protein.append(protein)
                        global_timepoint.append(timepoint)

                    global_mtoc.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                    # print(len(spot_by_quad[non_mtoc_spot][:,0].flatten()))
                    for i in range(0, 270, 3):
                        global_nmtoc.append(np.mean(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 2]))
                    # global_nmtoc.extend(spot_by_quad[non_mtoc_spot][:,0].flatten())

                    if mtoc_quad_j == 1:
                        # print(spot_by_quad[mtoc_spot][:,0].flatten())
                        global_mtoc_leading.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                    else:
                        for i in range(90):
                            global_mtoc_leading.append(np.nan)

        df = pd.DataFrame(
            {'Image': global_image, 'Gene': global_protein, 'timepoint': global_timepoint, 'Non MTOC': global_nmtoc,
             'MTOC': global_mtoc,
             'MTOC leading edge': global_mtoc_leading},
            index=global_index)
        df.to_csv(path.analysis_dir + 'analysis_nocodazole/df/' +'global_mtoc_file_protein_all.csv')







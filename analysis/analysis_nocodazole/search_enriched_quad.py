#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import pandas as pd
pd.set_option('display.max_rows', 500)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])

def box_plot_concentrations(k1, k2, title, figname):
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
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    y_lim_max = np.max(data) + 0.2
    y_lim_min = np.min(data) - 0.2
    print(y_lim_max)
    ind = np.arange(N)
    width = 0.35
    rects1 = ax.bar(ind, data, width,
                    color=colors)
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(y_lim_min, y_lim_max)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    ax.set_xticks(ind)
    plt.legend([gene for gene in genes], loc='upper right')
    ax.legend(rects1, genes, prop={'size': 8})
    plt.savefig(figname, format='svg')
    plt.show()


def plot_bar_profile_2_series(data,data_noc, genes,noc_genes, y_limit, ylabel, figname, colors):
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    y_lim_max = np.max(data) + 0.2
    y_lim_min = np.min(data) - 0.4
    print(y_lim_max)
    ind = np.arange(N)
    width = 0.35
    rects1 = ax.bar(ind, data, width,
                    color=['grey','grey'])
    rects2 = ax.bar(ind+width, data_noc, width,
                    color=['#1E95BB','#F16C1B'])
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(y_lim_min, y_lim_max)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('arhgdia', 'pard3'))
    ax.legend((rects1[0], rects2[0]), ('control', 'Nocodazole'))
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
        from pylab import setp
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
                    spot_by_quad = idsc.search_mrna_quadrants(file_handler, second_file_handler, image)
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    mtoc_spot = spot_by_quad[:, :, 1] == 1
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
        molecule_type = ['/protein']
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
            timepoints = ["3h", "5h"]
            for timepoint in timepoints:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, protein, [timepoint])
                for image in image_list:
                    intensity_by_quad = idsc.search_protein_quadrants(file_handler, second_file_handler, mtoc_file_handler,protein, image)
                    mtoc_intensity = intensity_by_quad[:, :, 1] == 1
                    non_mtoc_intensity = intensity_by_quad[:, :, 1] == 0
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    for i in range(90):
                        global_index.append(image.split("/")[4] + "_" + str(i + 1))
                        global_image.append(image.split("/")[4])
                        global_protein.append(protein)
                        global_timepoint.append(timepoint)
                    global_mtoc.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                    for i in range(0, 270, 3):
                        global_nmtoc.append(np.mean(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 2]))
                    if mtoc_quad_j == 1:
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
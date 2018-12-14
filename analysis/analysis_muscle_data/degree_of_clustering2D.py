#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import interpolate
import src.constants as cst
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])


def compute_minimal_distance(segment_summed):
    for i in range(15):
        if segment_summed[i]!=0:
            return i


def keep_cell_mask_spots(spots,cell_mask):
    #print(len(spots))
    new_spots_list=[]
    for spot in spots:
        if cell_mask[spot[1],spot[0]]==1:
            new_spots_list.append(spot)
    #print(len(new_spots_list))
    return new_spots_list




if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    molecule_type = ['/mrna']
    genes=['actn2','gapdh']
    timepoints=['mature']
    colors = ['#0A3950', '#1E95BB', '#A1BA6D']


    #compute degree of clusteirng
    with h5py.File(basic_file_path, "a") as file_handler:
        base = math.log(0.5)
        mrna_median = []
        mrna_err = []
        for gene in genes:
            image_list=[]
            h_star_l = []
            for timepoint in timepoints:
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:


                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    # if not "/mrna/actn2/mature/03" in image:
                    #      continue
                    cell_mask=idsc.get_cell_mask(file_handler, image)
                    nucleus_mask=idsc.get_nucleus_mask(file_handler,image)
                    nucleus_centroid=idsc.get_multiple_nucleus_centroid(file_handler,image)
                    spots=idsc.get_spots(file_handler,image)
                    x_nuc=[]
                    for nuc_n in range(len(nucleus_centroid)):
                        x_nuc.append(nucleus_centroid[nuc_n][0])
                    x_nuc=np.sort(x_nuc)
                    cell_mask[:, 0:x_nuc[0]] = 0
                    cell_mask[:, x_nuc[1]::] = 0
                    #radius_dist=int((x_nuc[1]-x_nuc[0])/2)
                    radius_dist=50
                    spots = keep_cell_mask_spots(spots, cell_mask)
                    spots = np.array(spots).reshape((len(spots), 3))

                    # plt.imshow(cell_mask)
                    # xs = spots[:, 0]
                    # ys = spots[:, 1]
                    # from skimage import measure
                    # contours = measure.find_contours(nucleus_mask, 0.8)
                    # for n, contour in enumerate(contours):
                    #     plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
                    # plt.scatter(xs, ys, color='blue', marker="o", facecolors='none', linewidths=0.5)
                    #
                    # plt.show()
                    # sys.exit()
                    h_star=helps.clustering_index_point_process_2d(spots,cell_mask,radius_dist)
                    #print(h_star)
                    #descriptor = image + '/h_star'
                    # output_file_handler[descriptor][:] = h_star
                    #file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)
                    d = np.array(h_star[h_star > 1] - 1).sum()
                    h_star_l.append(d)
            mrna_median.append(math.log(np.median(h_star_l)) - base)
            # mrna_median.append(numpy.median(dof))
            err = np.median(np.abs(np.tile(np.median(h_star_l), (1, len(h_star_l))) - h_star_l))
            error_median = math.log(np.median(h_star_l) + err)
            error_median = error_median - math.log(np.median(h_star_l)) - base
            mrna_err.append(error_median)
        fig = plt.figure()
        ax = plt.axes()
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        N = len(genes)
        ind = np.arange(N)
        width = 0.35
        # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
        colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']

        ## the bars
        rects1 = ax.bar(ind, mrna_median, width,
                        color=colors,
                        yerr=mrna_err,
                        error_kw=dict(elinewidth=1, ecolor='black'))
        # axes and labels
        ax.set_xlim(-width, len(ind) + width)
        ax.set_ylim(0, 20)
        ax.set_ylabel('Degree of clustering(log)')
        ax.set_title('Mrna degree of clustering')
        xTickMarks = ["" for i in range(0, 2)]
        ax.set_xticks(ind)
        xtickNames = ax.set_xticklabels(xTickMarks)
        plt.legend([gene for gene in genes], loc='upper right')
        ax.legend(rects1, genes)
        figname = path.analysis_dir + 'analysis_muscle_data/figures/mrna_degree_of_clustering_mature.svg'
        fig.savefig(figname)
        plt.show()
        print(h_star_l)

    molecule_type = ['/mrna']
    genes = ['actn2']
    timepoints = ['immature']
    colors = ['#0A3950']

    # compute degree of clusteirng
    with h5py.File(basic_file_path, "a") as file_handler:
        base = math.log(0.5)
        mrna_median = []
        mrna_err = []
        for gene in genes:
            image_list = []
            h_star_l = []
            for timepoint in timepoints:
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    nucleus_mask=idsc.get_nucleus_mask(file_handler,image)
                    cell_mask=idsc.get_cell_mask(file_handler, image)
                    #z_lines_mask = idsc.get_z_lines_mask(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                    x_nuc = []
                    for nuc_n in range(len(nucleus_centroid)):
                        x_nuc.append(nucleus_centroid[nuc_n][0])
                    x_nuc = np.sort(x_nuc)
                    cell_mask[:, 0:x_nuc[0]] = 0
                    cell_mask[:, x_nuc[1]::] = 0
                    radius_dist = int((x_nuc[1] - x_nuc[0]) / 2)
                    radius_dist = 50
                    spots = keep_cell_mask_spots(spots, cell_mask)
                    spots = np.array(spots).reshape((len(spots), 3))
                    # plt.imshow(cell_mask)
                    # xs = spots[:, 0]
                    # ys = spots[:, 1]
                    # from skimage import measure
                    # contours = measure.find_contours(nucleus_mask, 0.8)
                    # for n, contour in enumerate(contours):
                    #     plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
                    # plt.scatter(xs, ys, color='blue', marker="o", facecolors='none', linewidths=0.5)
                    #
                    # plt.show()


                    h_star = helps.clustering_index_point_process_2d(spots, cell_mask,radius_dist)
                    descriptor = image + '/h_star'
                    # output_file_handler[descriptor][:] = h_star
                    #file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)
                    d = np.array(h_star[h_star > 1] - 1).sum()
                    h_star_l.append(d)
            mrna_median.append(math.log(np.median(h_star_l)) - base)
            # mrna_median.append(numpy.median(dof))
            err = np.median(np.abs(np.tile(np.median(h_star_l), (1, len(h_star_l))) - h_star_l))
            error_median = math.log(np.median(h_star_l) + err)
            error_median = error_median - math.log(np.median(h_star_l)) - base
            mrna_err.append(error_median)






        fig = plt.figure()
        ax = plt.axes()
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        N = len(genes)
        ind = np.arange(N)
        width = 0.35
        # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
        #colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']

        ## the bars
        rects1 = ax.bar(ind, mrna_median, width,
                        color=colors,
                        yerr=mrna_err,
                        error_kw=dict(elinewidth=1, ecolor='black'))
        # axes and labels
        ax.set_xlim(-width, len(ind) + width)
        ax.set_ylim(0, 20)
        ax.set_ylabel('Degree of clustering(log)')
        ax.set_title('Mrna degree of clustering')
        xTickMarks = ["" for i in range(0, 1)]
        ax.set_xticks(ind)
        xtickNames = ax.set_xticklabels(xTickMarks)
        plt.legend([gene for gene in genes], loc='upper right')
        ax.legend(rects1, genes)
        figname = path.analysis_dir + 'analysis_muscle_data/figures/mrna_degree_of_clustering_immature.svg'
        fig.savefig(figname)
        plt.show()
        print(h_star_l)


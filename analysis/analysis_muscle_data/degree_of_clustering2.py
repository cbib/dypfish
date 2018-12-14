#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
from numpy import matlib
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


def keep_cell_mask_spots(spots,cell_mask,z_lines):
    new_spots_list=[]
    for spot in spots:
        #print(spot[2])
        if cell_mask[spot[1],spot[0]]==1:

            new_spots_list.append(spot)
    return new_spots_list


def get_quantized_grid(q, Qx, Qy):

    tmp_x = np.matrix(np.arange(Qx))
    tmp_y = np.matrix(np.arange(Qy))
    #print(tmp.transpose())

    #    qxs = np.kron(tmp,Q)
    qxs = matlib.repmat(tmp_x.transpose(), 1, Qx)
    qys = matlib.repmat(tmp_y, Qy, 1)
    #print(qxs)
    qxs = np.kron(qxs, np.ones((q, q)))
    print(qxs)
    plt.imshow(qxs)
    plt.show()
    qys = np.kron(qys, np.ones((q, q)))
    # print qys
    return qxs, qys

def get_variance(test):
    print(test)

    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot / N)

def compute_cell_mask_between_nucleus_centroid(cell_mask,nucleus_centroid,nuc_dist,cell_masks,nucs_dist):
    #x_nuc = []
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    #print(nucleus_centroid)
    im_mask=[]
    nucs=[]
    for nuc_n in range(len(nucleus_centroid)-1):
        cell_mask_copy=cell_mask.copy()
        #print(nucleus_centroid[nuc_n][0])
        #print(nucleus_centroid[nuc_n+1][0])
        nuc_dist.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])
        nucs.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])

        cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
        cell_mask_copy[:, nucleus_centroid[nuc_n+1][0]::] = 0
        im_mask.append(cell_mask_copy)
        # plt.imshow(cell_mask_copy)
        # plt.show()
    nucs_dist.append(nucs)
    cell_masks.append(im_mask)
    return nuc_dist,nucs_dist,cell_masks


def search_best_centroid(nucleus_centroid):
    x_nuc = []
    print(nucleus_centroid)
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    x_nuc.append(nucleus_centroid[0][0])
    x_nuc.append(nucleus_centroid[len(nucleus_centroid)-1][0])


    #for nuc_n in range(len(nucleus_centroid)):
    #    x_nuc.append(nucleus_centroid[nuc_n][0])
    print(x_nuc)
    return x_nuc


    # x_nuc = np.sort(x_nuc)
    # print(x_nuc)
    # left_nuc = x_nuc[0]
    # right_nuc = x_nuc[1]
    # cell_width = right_nuc - left_nuc

def show_descriptors(cell_mask,spots,nucleus_mask):

    plt.imshow(cell_mask)
    xs = spots[:, 0]
    ys = spots[:, 1]
    plt.grid(True)
    from skimage import measure
    # contours = measure.find_contours(nucleus_mask, 0.8)
    # for n, contour in enumerate(contours):
    #     plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)

    plt.show()

def compute_cell_mask_3d(cell_mask,z_lines):
    x_dim=cell_mask.shape[0]
    y_dim = cell_mask.shape[1]
    cell_mask_3d = np.zeros((x_dim, y_dim, len(z_lines)))

    for slice in range(len(z_lines)):
        cell_mask_3d[:, :, slice] = cell_mask

    return cell_mask_3d

def reject_outliers(data):
    index_cpt=0
    indexes=[]
    tmp_list=[]
    for i in data:
        for j in i:
            tmp_list.append(j)


    u=np.mean(tmp_list)
    for i in data:
        print(i)
        for j in i:
            if u- 200 < j < u + 200:
                indexes.append(index_cpt)
        index_cpt += 1

    #filtered = [e for e in data if (np.mean(data)-200 < e < np.mean(data)+200)]
    #print(data[filtered])
    return indexes

    #return data[abs(np.mean(data)-200 < data  < np.mean(data)+200)]

def compute_degree_of_clustering(spots_reduced,mask,z_lines):
    cell_mask_3d = compute_cell_mask_3d(mask, z_lines)
    h_star = helps.clustering_index_point_process(spots_reduced, cell_mask_3d)
    d = np.array(h_star[h_star > 1] - 1).sum()
    return d


def reduce_z_line_mask(z_lines,spots):

    cpt_z = 1
    z_lines_idx=[]
    for z_line_mask in z_lines:
        #print(cpt_z)
        spots_reduced = spots[spots[:, 2] == cpt_z]
        xs = spots_reduced[:, 0]
        ys = spots_reduced[:, 1]
        if len(spots_reduced) > 30:
            z_lines_idx.append(cpt_z)
            print("zline: "+ str(cpt_z))
            print("number of spots: " + str(len(spots_reduced)))

            # fig = plt.figures(figsize=(40, 10))
            # plt.imshow(z_line_mask)
            # plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)
            # plt.show()
        cpt_z += 1
    return z_lines_idx

if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    muscle_rebuild_file_path = path.analysis_data_dir + 'secondary_muscle_data.h5'

    molecule_type = ['/mrna']
    genes=['actn2','gapdh']
    timepoints=['mature']
    colors = ['#0A3950', '#1E95BB', '#A1BA6D']


    #compute quadrat

    all_median_profiles=[]

    with h5py.File(basic_file_path, "a") as file_handler, h5py.File(muscle_rebuild_file_path, "a") as muscle_file_handler:
        base = math.log(0.5)
        mrna_median = []
        mrna_err = []


        for gene in genes:
            image_list=[]
            h_star_l = []
            for timepoint in timepoints:
                total_image=0
                image_counter=0
                total_profile = []
                z_lines=[]
                for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(sec_image)
                    print(im)
                    image_number=im.split("_")[0]
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                    print(image)
                    spots = idsc.get_spots(muscle_file_handler, sec_image)
                    z_lines=idsc.get_z_lines_masks(file_handler, image)
                    cell_mask=idsc.get_cell_mask(muscle_file_handler, sec_image)

                    z_lines_idx = reduce_z_line_mask(z_lines, spots)
                    spots_reduced = spots[z_lines_idx[0] <= spots[:, 2]]
                    spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
                    print(spots_reduced)
                    h_star = helps.clustering_index_point_process_2d(spots_reduced, cell_mask, 100)
                    # print(h_star)
                    d = np.array(h_star[h_star > 1] - 1).sum()
                    h_star_l.append(d)


                    # slice_rejected = len(z_lines) / 5
                    # for i in range(slice_rejected, len(z_lines) - slice_rejected):
                    #     spots_by_slice=spots[spots[:,2]==i]
                    #     if spots_by_slice!=[]:
                    #         h_star = helps.clustering_index_point_process_2d(spots_by_slice, cell_mask, 100)
                    #         #print(h_star)
                    #         d = np.array(h_star[h_star > 1] - 1).sum()
                    #         h_star_l.append(d)


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
        figname = path.analysis_dir + 'analysis_muscle_data/figures/mrna_degree_of_clustering_mature_by_slice2D_radius'+str(cst.MAX_CELL_RADIUS)+'_bad_slice_removed.svg'
        fig.savefig(figname)
        #plt.show()
        print(h_star_l)



    molecule_type = ['/mrna']
    genes = ['actn2']
    timepoints = ['immature']
    colors = ['#0A3950', '#1E95BB', '#A1BA6D']

    # compute quadrat

    all_median_profiles = []

    with h5py.File(basic_file_path, "a") as file_handler, h5py.File(muscle_rebuild_file_path,
                                                                    "a") as muscle_file_handler:
        base = math.log(0.5)
        mrna_median = []
        mrna_err = []

        for gene in genes:
            image_list = []
            h_star_l = []
            for timepoint in timepoints:
                total_image = 0
                image_counter = 0
                total_profile = []
                z_lines = []
                for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(sec_image)
                    print(im)
                    image_number = im.split("_")[0]
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                    print(image)
                    spots = idsc.get_spots(muscle_file_handler, sec_image)
                    z_lines = idsc.get_z_lines_masks(file_handler, image)
                    cell_mask = idsc.get_cell_mask(muscle_file_handler, sec_image)

                    z_lines_idx = reduce_z_line_mask(z_lines, spots)
                    spots_reduced = spots[z_lines_idx[0] <= spots[:, 2]]
                    spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
                    print(spots_reduced)
                    h_star = helps.clustering_index_point_process_2d(spots_reduced, cell_mask, 100)
                    # print(h_star)
                    d = np.array(h_star[h_star > 1] - 1).sum()
                    h_star_l.append(d)

                    # slice_rejected=len(z_lines)/5
                    # for i in range(slice_rejected, len(z_lines) - slice_rejected):
                    #     spots_by_slice = spots[spots[:, 2] == i]
                    #     if len(spots_by_slice) != 0:
                    #         h_star = helps.clustering_index_point_process_2d(spots_by_slice, cell_mask, 100)
                    #         # print(h_star)
                    #         d = np.array(h_star[h_star > 1] - 1).sum()
                    #         h_star_l.append(d)

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
        figname = path.analysis_dir + 'analysis_muscle_data/figures/mrna_degree_of_clustering_immature_by_slice2D_radius'+str(cst.MAX_CELL_RADIUS)+'_bad_slice_removed.svg'
        fig.savefig(figname)
        #plt.show()
        print(h_star_l)
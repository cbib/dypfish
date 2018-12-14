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

def compute_cell_mask_between_nucleus_centroid(cell_mask,nucleus_centroid,nuc_dist,cell_masks,nucs_dist,nucs_pos,nucleus_mask):
    #x_nuc = []
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    #print(nucleus_centroid)
    im_mask=[]
    nucs=[]
    nuc_pos=[]


    for nuc_n in range(len(nucleus_centroid)-1):

        cell_mask_copy=cell_mask.copy()
        #print(nucleus_centroid[nuc_n][0])
        #print(nucleus_centroid[nuc_n+1][0])
        nuc_pos.append([nucleus_centroid[nuc_n][0],nucleus_centroid[nuc_n+1][0]])
        nuc_dist.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])
        nucs.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])

        cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
        cell_mask_copy[:, nucleus_centroid[nuc_n+1][0]::] = 0
        im_mask.append(cell_mask_copy)
        #plt.imshow(cell_mask_copy)
        #plt.show()
    nucs_dist.append(nucs)
    cell_masks.append(im_mask)
    nucs_pos.append(nuc_pos)

    return nuc_dist,nucs_dist,cell_masks,nucs_pos


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

if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    muscle_rebuild_file_path = path.analysis_data_dir + 'secondary_muscle_data.h5'

    molecule_type = ['/mrna']
    genes=['actn2']#,'gapdh']
    timepoints=['immature']
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
                nuc_dist = []
                nucs_pos = []
                cell_masks=[]
                nucs_dist=[]
                # compute new cell mask and cell mask width
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    nuc_dist, nucs_dist,cell_masks, nucs_pos = compute_cell_mask_between_nucleus_centroid(cell_mask, nucleus_centroid,nuc_dist,cell_masks,nucs_dist,nucs_pos,nucleus_mask)
                    total_image+=1
                sys.exit()
                hx, hy, _ =plt.hist(nuc_dist)
                bin_max = np.where(hx == hx.max())[0]

                hist_mod=hy[bin_max]
                print(hist_mod)
                if len(hist_mod) >1:
                    hist_mod=np.mean(hist_mod)
                image_counter=0
                total_profile = []

                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    z_lines=idsc.get_z_lines_masks(file_handler, image)
                    cell_masks_im = cell_masks[image_counter]
                    nuc_dist_im = nucs_dist[image_counter]
                    nuc_pos_im = nucs_pos[image_counter]
                    print(nuc_pos_im)
                    print(nuc_dist_im)




                    # cpt_z=1
                    # fig = plt.figures(figsize=(40, 10))
                    # for z_line_mask in z_lines:
                    #     plt.imshow(z_line_mask)
                    #     from skimage import measure
                    #
                    #     contours = measure.find_contours(nucleus_mask, 0.8)
                    #     for n, contour in enumerate(contours):
                    #         plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
                    #     spots_reduced=spots[spots[:, 2]==cpt_z]
                    #     print(spots_reduced)
                    #     xs = spots_reduced[:, 0]
                    #     ys = spots_reduced[:, 1]
                    #     if spots_reduced != []:
                    #         plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)
                    #     plt.show()
                    #     cpt_z+=1


                    mask_count=0
                    for i in range(len(cell_masks_im)):
                        mask = cell_masks_im[i]
                        nuc_d = nuc_dist_im[i]
                        nuc_pos = nuc_pos_im[i]
                        print(nuc_pos)


                        if hist_mod - 250 < nuc_d < hist_mod + 250:
                            #print(spots)

                            spots_reduced = keep_cell_mask_spots(spots, mask, z_lines)
                            spots_reduced = np.array(spots_reduced).reshape((len(spots_reduced), 3))
                            #print(spots_reduced)
                            #show_descriptors(mask,spots_reduced,nucleus_mask)
                            #print(mask.shape)
                            #spots_reduced=spots_reduced[(spots_reduced[:,2] > 15) & (spots_reduced[:,2] < (len(z_lines)-15))]
                            #print(spots_reduced)
                            #show_descriptors(mask, spots_reduced, nucleus_mask)
                            print("new dimension: ",nuc_pos[0]-50,str(nuc_pos[1]+50))
                            mask_c=mask.copy()
                            mask_reduced=mask_c[:,nuc_pos[0]-50:nuc_pos[1]+50]
                            spots_reduced[:,0]-=nuc_pos[0]-50
                            #print(spots_reduced)

                            #show_descriptors(mask_reduced,spots_reduced,nucleus_mask)
                            cell_mask_id='cell_mask'
                            new_im= molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im +'_'+str(mask_count)
                            descriptor = new_im +'/'+cell_mask_id
                            muscle_file_handler.create_dataset(descriptor, data=mask_reduced, dtype=np.int)
                            spots_id='spots'
                            descriptor = new_im+ '/'+spots_id
                            muscle_file_handler.create_dataset(descriptor, data=spots_reduced, dtype=np.float64, maxshape=(None, 3))

                            #compute spots minimal distance to z-lines
                            spots_count = 0
                            mask_count+=1
                    image_counter+=1
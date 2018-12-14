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


def keep_cell_mask_spots(spots,cell_mask):
    #print(len(spots))
    new_spots_list=[]
    for spot in spots:
        if cell_mask[spot[1],spot[0]]==1:
            new_spots_list.append(spot)
    #print(len(new_spots_list))
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

def compute_cell_mask_between_nucleus_centroid(cell_mask,nucleus_centroid,nuc_dist,cell_masks):
    #x_nuc = []
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    print(nucleus_centroid)
    im_mask=[]
    for nuc_n in range(len(nucleus_centroid)-1):
        cell_mask_copy=cell_mask.copy()
        #print(nucleus_centroid[nuc_n][0])
        #print(nucleus_centroid[nuc_n+1][0])
        nuc_dist.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])
        cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
        cell_mask_copy[:, nucleus_centroid[nuc_n+1][0]::] = 0
        im_mask.append(cell_mask_copy)
        # plt.imshow(cell_mask_copy)
        # plt.show()
    cell_masks.append(im_mask)
    return nuc_dist,cell_masks


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
    contours = measure.find_contours(nucleus_mask, 0.8)
    for n, contour in enumerate(contours):
        plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)

    plt.show()



if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    molecule_type = ['/mrna']
    genes=['actn2','gapdh']
    timepoints=['mature']
    colors = ['#0A3950', '#1E95BB', '#A1BA6D']


    #compute quadrat


    with h5py.File(basic_file_path, "a") as file_handler:
        for gene in genes:
            image_list=[]
            h_star_l = []
            for timepoint in timepoints:
                total_image=0
                nuc_dist = []
                cell_masks=[]
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    #show_descriptors(cell_mask,spots,nucleus_mask)
                    nuc_dist, cell_masks = compute_cell_mask_between_nucleus_centroid(cell_mask, nucleus_centroid,
                                                                                      nuc_dist,cell_masks)

                    total_image+=1
                #sys.exit()
                print(total_image)
                print(nuc_dist)
                print(len(cell_masks))
                image_counter=0
                # for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                #     image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                #     nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                #     nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                #     spots = idsc.get_spots(file_handler, image)
                #     cell_masks_im=cell_masks[image_counter]
                #     image_counter+=1
                #     for cell_mask in cell_masks_im:
                #         spots_reduced = keep_cell_mask_spots(spots, cell_mask)
                #         spots_reduced = np.array(spots_reduced).reshape((len(spots_reduced), 3))
                #         #show_descriptors(cell_mask,spots_reduced,nucleus_mask)
                #


                fig=plt.figure(figsize=(40,10))
                fig.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)

                image_cpt = 0

                #grid_medi
                nuc_dist=[]
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    # if not "/mrna/actn2/mature/03" in image:
                    #      continue
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    nucleus_centroid = idsc.get_multiple_nucleus_centroid(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    #z_lines_mask = idsc.get_z_lines_mask(file_handler, image)
                    # cytoplasm_mask=z_lines_mask[z_lines_mask==1]=2
                    # cytoplasm_mask=z_lines_mask[z_lines_mask == 0] = 1
                    # cytoplasm_mask=z_lines_mask[z_lines_mask == 2] = 0
                    # cytoplasm_mask=z_lines_mask[cell_mask==0]=0
                    spots = idsc.get_spots(file_handler, image)
                    x_nuc = []
                    for nuc_n in range(len(nucleus_centroid)):
                        x_nuc.append(nucleus_centroid[nuc_n][0])
                    x_nuc = np.sort(x_nuc)
                    left_nuc=x_nuc[0]
                    right_nuc=x_nuc[1]
                    cell_width = right_nuc - left_nuc
                    if cell_width<300:
                        left_nuc = x_nuc[1]
                        right_nuc = x_nuc[2]
                        cell_width = right_nuc - left_nuc

                    cell_mask[:, 0:left_nuc] = 0
                    cell_mask[:, right_nuc::] = 0
                    spots = keep_cell_mask_spots(spots, cell_mask)
                    spots = np.array(spots).reshape((len(spots), 3))


                    #cell_area = cell_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)

                    #show_descriptors(cell_mask,spots,nucleus_mask)

                    #quadrat_edge=int(math.ceil(math.sqrt(2)*(cell_mask.sum()/len(spots))))
                    band_width=20.0
                    quadrat_edge=cell_width/band_width
                    grid_1d=np.zeros((int(band_width)))
                    for spot in spots:
                        x = int(np.floor((spot[0]-left_nuc)/quadrat_edge))-1
                        grid_1d[x]+=1
                    grid=[]
                    for val in grid_1d:
                        if val !=0:
                            grid.append(val)
                    grid_mat=np.matrix(grid).reshape((1,len(grid)))
                    ax1 = plt.subplot2grid((total_image, 10), (image_cpt, 0), colspan=10)
                    #plt.imshow(grid_mat)
                    major_ticks = np.arange(0, len(grid), 1)
                    ax1.set_xticks(major_ticks)
                    ax1.grid(which='major', alpha=0.2)
                    plt.imshow(cell_mask)
                    #plt.show()
                    #plt.imshow(nucleus_mask)



                    #test=[]
                    # for val in grid.flatten():
                    #     if val !=-1:
                    #         test.append(val)
                    #var=get_variance(grid_1d.flatten())

                    #print(var/np.mean(grid_1d.flatten()))


                    #qxs, qys = get_quantized_grid(quadrat_edge, cell_mask.shape[1],cell_mask.shape[0])

                    #z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                    #cytoplasm_area = cell_area - z_lines_area
                    image_cpt+=1

                fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                plt.show()
                # plt.hist(nuc_dist)
                # plt.show()
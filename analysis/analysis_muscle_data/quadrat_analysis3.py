#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import matplotlib.patches as pct
from matplotlib.patches import Ellipse
from matplotlib import gridspec

from numpy import matlib
import math
import pandas as pd
from scipy import interpolate
import src.constants as cst
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.plot as plot
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
    ##print(len(spots))
    new_spots_list=[]
    for spot in spots:
        if cell_mask[spot[1],spot[0]]==1:
            new_spots_list.append(spot)
    ##print(len(new_spots_list))
    return new_spots_list


def get_quantized_grid(q, Qx, Qy):

    tmp_x = np.matrix(np.arange(Qx))
    tmp_y = np.matrix(np.arange(Qy))
    ##print(tmp.transpose())

    #    qxs = np.kron(tmp,Q)
    qxs = matlib.repmat(tmp_x.transpose(), 1, Qx)
    qys = matlib.repmat(tmp_y, Qy, 1)
    ##print(qxs)
    qxs = np.kron(qxs, np.ones((q, q)))
    ##print(qxs)
    plt.imshow(qxs)
    plt.show()
    qys = np.kron(qys, np.ones((q, q)))
    # #print qys
    return qxs, qys

def get_variance(test):
    #print(test)

    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot / N)

def compute_squared_diff_sum(test):
    #print(test)

    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot)

def compute_cell_mask_between_nucleus_centroid(cell_mask,nucleus_centroid,nuc_dist,cell_masks):
    #x_nuc = []
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    ##print(nucleus_centroid)
    im_mask=[]
    for nuc_n in range(len(nucleus_centroid)-1):
        cell_mask_copy=cell_mask.copy()
        ##print(nucleus_centroid[nuc_n][0])
        ##print(nucleus_centroid[nuc_n+1][0])
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
    #print(nucleus_centroid)
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    x_nuc.append(nucleus_centroid[0][0])
    x_nuc.append(nucleus_centroid[len(nucleus_centroid)-1][0])


    #for nuc_n in range(len(nucleus_centroid)):
    #    x_nuc.append(nucleus_centroid[nuc_n][0])
    ##print(x_nuc)
    return x_nuc


    # x_nuc = np.sort(x_nuc)
    # #print(x_nuc)
    # left_nuc = x_nuc[0]
    # right_nuc = x_nuc[1]
    # cell_width = right_nuc - left_nuc

def show_descriptors(cell_mask,spots,nucleus_mask):
    fig1 = plt.figure()
    plt.imshow(cell_mask)
    xs = spots[:, 0]
    ys = spots[:, 1]
    plt.grid(True)
    from skimage import measure
    contours = measure.find_contours(nucleus_mask, 0.8)
    for n, contour in enumerate(contours):
        ##print(contour[:, 1]-50)
        plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)

    plt.show()
    plt.close()

def reduce_z_line_mask(z_lines,spots):

    cpt_z = 1
    z_lines_idx=[]
    for z_line_mask in z_lines:
        ##print(cpt_z)
        spots_reduced = spots[spots[:, 2] == cpt_z]
        ##print(len(spots_reduced))

        xs = spots_reduced[:, 0]
        ys = spots_reduced[:, 1]
        if len(spots_reduced) > 25:
            z_lines_idx.append(cpt_z)
            ##print("zline: "+ str(cpt_z))
            ##print("number of spots: " + str(len(spots_reduced)))

            # fig = plt.figures(figsize=(40, 10))
            # plt.imshow(z_line_mask)
            # plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)
            # plt.show()
        cpt_z += 1
    return z_lines_idx


def get_cell_mask_width(cell_mask):
    ##print(cell_mask.shape)
    ##print(cell_mask)
    inx=[]
    for i in range(cell_mask.shape[1]):
        pix_val=cell_mask[int(cell_mask.shape[0]/2),i]
        if pix_val==1:
            inx.append(i)
    ##print(inx)

def get_cell_area(cell_mask):
    cell_area=np.sum(cell_mask==1)
    ##print(cell_area)
    return cell_area


def create_circle():
	circle= plt.Circle((0,0), radius= 5)
	return circle

def show_shape(patch):
	ax=plt.gca()
	ax.add_patch(patch)
	plt.axis('scaled')
	plt.show()


if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    muscle_rebuild_file_path = path.analysis_data_dir + 'secondary_muscle_data.h5'
    colors = ['#0A3950', '#A1BA6D','#1E95BB']
    molecule_type = ['/mrna']
    genes=['actn2-mature','gapdh-mature','actn2-immature']
    #timepoints=['mature']

    #compute quadrat
    with h5py.File(basic_file_path, "a") as file_handler, h5py.File(muscle_rebuild_file_path, "a") as muscle_file_handler:
        vmr_values_sig = []
        for gene in genes:
            print(gene)
            image_list=[]
            h_star_l = []
            [gene,timepoint] = gene.split("-")
            #timepoint=gene.split("-")[1]
            total_image=0
            nuc_dist = []
            cell_masks=[]
            image_counter=0
            image_cpt = 0
            for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                total_image += 1
            w = 40
            h = 2*total_image

            for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                #print(im)
                sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                #print(sec_image)
                image_number = im.split("_")[0]
                image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                #print(image)

                spots_reduced = idsc.get_spots(muscle_file_handler, sec_image)
                #nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                nucleus_mask = idsc.get_nucleus_mask(muscle_file_handler, sec_image)
                cell_mask = idsc.get_cell_mask(muscle_file_handler, sec_image)
                #show_descriptors(cell_mask,spots_reduced,nucleus_mask)


            fig = plt.figure()
            fig.set_size_inches(w, h)


            #fig, axes = plt.subplots(nrows=total_image*2, ncols=1)
            #axes = axes.flatten()
            vmr_values = []
            gs = gridspec.GridSpec(total_image * 6,1)

            for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                image_number = im.split("_")[0]
                image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                spots_reduced = idsc.get_spots(muscle_file_handler, sec_image)
                #all_spots = idsc.get_spots(file_handler, image)
                z_lines = idsc.get_z_lines_masks(file_handler, image)
                sums = np.zeros((len(z_lines), 1))
                for n in range(len(z_lines)):
                    slice = z_lines[n]
                    sums[n] = slice.sum()
                z_level=sums.argmax()
                # plt.imshow(z_lines[z_level])
                # plt.show()
                # for z_line in z_lines:
                #     plt.imshow(z_line)
                #     plt.show()
                nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                cell_mask = idsc.get_cell_mask(muscle_file_handler, sec_image)
                #show_descriptors(cell_mask,spots_reduced,nucleus_mask)
                z_lines_idx=reduce_z_line_mask(z_lines, spots_reduced)
                spots_reduced=spots_reduced[z_lines_idx[0] <= spots_reduced[:,2] ]
                spots_reduced=spots_reduced[spots_reduced[:,2] <= z_lines_idx[len(z_lines_idx)-1]]
                cell_area = get_cell_area(cell_mask)
                #print(cell_area)
                spot_surfacic_density = len(spots_reduced) / float(cell_area)
                cell_width=cell_mask.shape[1] -240
                #print(cell_width)
                band_n = 5.0
                #band_n=np.ceil(cell_width/20)
                #print(band_n)
                quadrat_edge = cell_width / band_n
                #quadrat_edge = cell_area / band_n

                #print(quadrat_edge)
                #print(spot_surfacic_density)
                grid_1d = np.zeros((int(band_n)))
                for spot in spots_reduced:
                    if spot[0] > 120 and spot[0]< cell_mask.shape[1]- 120:
                        x = int(np.floor((spot[0]-120)/ quadrat_edge))
                        grid_1d[x] += 1
                grid = []
                for val in grid_1d:
                    grid.append(val)
                #print(grid)
                #plt.show()
                grid_mat = np.matrix(grid).reshape((1, len(grid)))
                #print(grid_mat)
                grid_mat/=quadrat_edge
                #print(grid_mat)
                grid_mat/=spot_surfacic_density
                #print(grid_mat)
                #sys.exit()
                grid_list=list(np.array(grid_mat)[0])
                var = compute_squared_diff_sum(grid_list)
                # 30.144 is the value of chi-square for degree of freedom  19
                vmr_values.append((var / np.mean(grid_list))-30.144)

                #gs.update(hspace=0.5)
                ax0 = plt.subplot(gs[image_cpt])

                if image_cpt == 0:

                    im = ax0.imshow(grid_mat, cmap='coolwarm')
                    #clim = im.properties()['clim']
                    #cbarobj = plt.colorbar(im)
                    clim=im.set_clim(0, 100)
                    ax0.imshow(grid_mat,  cmap='coolwarm',clim=clim)
                else:
                    ax0.imshow(grid_mat,  cmap='coolwarm',clim=clim)

                circle1 = Ellipse((0, 0), width=0.6, height=0.8, clip_on=True, color='lightgreen')
                circle2 = Ellipse((band_n-1, 0), width=0.6, height=0.8, clip_on=True, color='lightgreen')
                ax0.add_artist(circle1)
                ax0.add_artist(circle2)
                #ax0.set_ylim(())
                #ax0.set_xlim((0, 19))
                ax0.set_yticks([])
                #ax0.set_xticks([])

                image_cpt+=1

                ax2 = plt.subplot(gs[image_cpt])
                ax2.imshow(z_lines[z_level], aspect='auto',extent=(0,band_n-1,0,140))
                #ax2.imshow(z_lines[z_level], aspect='auto')

                ax2.set_xlim((0,band_n-1))
                ax2.set_yticks([])
                image_cpt+=1

                ax1 = plt.subplot(gs[image_cpt,image_cpt+2:],sharex=ax0)
                plt.setp(ax0.get_xticklabels(), visible=False)
                #plt.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=1.0)
                #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.1, hspace=0.5)

                #ax1.plot(grid/np.mean(grid))
                x_mrna = np.arange(0, band_n, 0.5)
                #grid=grid/np.mean(grid)
                #print(np.array(grid_mat).flatten())
                fact=np.max(np.array(grid_mat))/10
                #print((np.array(grid_mat).flatten()/fact).astype(int))
                data=(np.array(grid_mat).flatten()/fact).astype(int)+1
                spl = interpolate.UnivariateSpline(np.arange(0, band_n, 1), data)

                m_y_new = spl(x_mrna)
                # new_val = np.zeros((25, 1))
                # for val in np.array(grid_mat).flatten():
                #     print(int(val / 10))
                #     new_val[int(val / 10), :] += 1
                ax1.plot(x_mrna,m_y_new)

                #print((np.array(grid_mat).flatten()/10).astype(int))
                #ax1.plot((np.array(grid_mat).flatten()/10).astype(int))
                #print(new_val.flatten())
                #ax1.plot(np.array(grid_mat).flatten())
                #ax1.plot(grid)



                #sys.exit()


                #ax1.plot(grid)
                circle1 = Ellipse((0, 5.5),width=1, height=11,clip_on=True,color='lightgreen')
                circle2 = Ellipse((band_n-1, 5.5),width=1, height=11,clip_on=True,color='lightgreen')

                #circle2 = plt.Circle((99, 10),radius=5, clip_on=False, color='blue')
                ax1.add_artist(circle1)
                ax1.add_artist(circle2)
                ax1.set_ylim((0,11))

                ax1.set_xlim((0, band_n-1))
                ax1.set_yticks([])
                #ax1.set_xticks((0,int(band_n)))


                image_cpt+=4
            #plt.show()
            print(vmr_values)
            vmr_values_sig.append(np.sum(np.array(vmr_values)>0)/float(len(vmr_values)))
            plt.savefig("quadrat_analysis_"+gene+"_"+timepoint+str(int(band_n))+"_zline.png")
            plt.close()

        figname="clustering_index.png"
        plot.bar_profile(vmr_values_sig, genes, 1.0, 'index of clustering', figname,colors)

    # genes = ['actn2']
    # timepoints = ['immature']
    # # compute quadrat
    #
    # with h5py.File(basic_file_path, "a") as file_handler, h5py.File(muscle_rebuild_file_path,
    #                                                                 "a") as muscle_file_handler:
    #     for gene in genes:
    #         image_list = []
    #         h_star_l = []
    #         #print(gene)
    #         for timepoint in timepoints:
    #             total_image = 0
    #             nuc_dist = []
    #             cell_masks = []
    #             image_counter = 0
    #             fig = plt.figures(figsize=(40, 10))
    #             fig.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    #             image_cpt = 0
    #             for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
    #                 total_image += 1
    #
    #             for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
    #                 sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
    #                 #print(sec_image)
    #                 # #print(im)
    #                 image_number = im.split("_")[0]
    #                 image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
    #                 # #print(image)
    #                 spots_reduced = idsc.get_spots(muscle_file_handler, sec_image)
    #                 # #print(spots_reduced)
    #                 # all_spots = idsc.get_spots(file_handler, image)
    #                 # #print(all_spots)
    #                 z_lines = idsc.get_z_lines_masks(file_handler, image)
    #                 nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
    #                 cell_mask = idsc.get_cell_mask(muscle_file_handler, sec_image)
    #                 z_lines_idx = reduce_z_line_mask(z_lines, spots_reduced)
    #                 spots_reduced = spots_reduced[z_lines_idx[0] <= spots_reduced[:, 2]]
    #                 spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
    #                 cell_area = get_cell_area(cell_mask)
    #                 spot_volumic_density = len(spots_reduced) / float(cell_area)
    #                 # #print(spots_reduced)
    #                 # show_descriptors(cell_mask,spots_reduced,nucleus_mask)
    #                 # cell_width=get_cell_mask_width(cell_mask)
    #                 cell_width = cell_mask.shape[1] - 240
    #                 band_width = 20.0
    #
    #                 quadrat_edge = cell_width / band_width
    #                 #print(quadrat_edge)
    #                 grid_1d = np.zeros((int(band_width)))
    #                 for spot in spots_reduced:
    #                     if spot[0] > 120 and spot[0]< cell_mask.shape[1]- 120:
    #
    #                         x = int(np.floor((spot[0] - 120) / quadrat_edge)) - 1
    #                         # #print(x)
    #                         grid_1d[x] += 1
    #                 grid = []
    #                 for val in grid_1d:
    #                     if val != 0:
    #                         grid.append(val)
    #                     else:
    #                         grid.append(1.0)
    #                 grid_mat = np.matrix(grid).reshape((1, len(grid)))
    #                 grid_mat /= quadrat_edge
    #                 grid_mat/=spot_volumic_density
    #                 ax1 = plt.subplot2grid((total_image, 10), (image_cpt, 0), colspan=10)
    #                 plt.imshow(grid_mat)
    #                 plt.clim(0, 100)
    #                 major_ticks = np.arange(0, len(grid), 1)
    #
    #                 ax1.set_xticks(major_ticks)
    #                 ax1.grid(which='major', alpha=0.2)
    #                 image_cpt += 1
    #
    #             fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    #             plt.show()
    #
    #
    #
    #
    #
    #

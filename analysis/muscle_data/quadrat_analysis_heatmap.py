#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import logging
import sys
import math
import argparse
import numpy as np
import h5py
import src.image_descriptors as idsc
import src.path as path
import src.plot as plot
from src.utils import check_dir, loadconfig

import matplotlib.pyplot as plt
from matplotlib import gridspec
from numpy import matlib
from skimage import measure

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def compute_minimal_distance(segment_summed):
    for i in range(15):
        if segment_summed[i]!=0:
            return i

def keep_cell_mask_spots(spots,cell_mask):
    new_spots_list=[]
    for spot in spots:
        if cell_mask[spot[1],spot[0]]==1:
            new_spots_list.append(spot)
    return new_spots_list

def get_quantized_grid(q, Qx, Qy):
    tmp_x = np.matrix(np.arange(Qx))
    tmp_y = np.matrix(np.arange(Qy))
    qxs = matlib.repmat(tmp_x.transpose(), 1, Qx)
    qys = matlib.repmat(tmp_y, Qy, 1)
    qxs = np.kron(qxs, np.ones((q, q)))
    plt.imshow(qxs)
    plt.show()
    qys = np.kron(qys, np.ones((q, q)))
    return qxs, qys

def get_variance(test):
    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot / N)

def compute_squared_diff_sum(test):
    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot)

def compute_cell_mask_between_nucleus_centroid(cell_mask,nucleus_centroid,nuc_dist,cell_masks):
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    im_mask=[]
    for nuc_n in range(len(nucleus_centroid)-1):
        cell_mask_copy=cell_mask.copy()
        nuc_dist.append(nucleus_centroid[nuc_n+1][0]-nucleus_centroid[nuc_n][0])
        cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
        cell_mask_copy[:, nucleus_centroid[nuc_n+1][0]::] = 0
        im_mask.append(cell_mask_copy)
    cell_masks.append(im_mask)
    return nuc_dist,cell_masks

def search_best_centroid(nucleus_centroid):
    x_nuc = []
    nucleus_centroid=np.sort(nucleus_centroid,axis=0)
    x_nuc.append(nucleus_centroid[0][0])
    x_nuc.append(nucleus_centroid[len(nucleus_centroid)-1][0])
    return x_nuc

def show_descriptors(cell_mask,spots,nucleus_mask):
    fig1 = plt.figure()
    plt.imshow(cell_mask)
    xs = spots[:, 0]
    ys = spots[:, 1]
    plt.grid(True)

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
        spots_reduced = spots[spots[:, 2] == cpt_z]
        xs = spots_reduced[:, 0]
        ys = spots_reduced[:, 1]
        if len(spots_reduced) > 25:
            z_lines_idx.append(cpt_z)
        cpt_z += 1
    return z_lines_idx

def get_cell_mask_width(cell_mask):
    inx=[]
    for i in range(cell_mask.shape[1]):
        pix_val=cell_mask[int(cell_mask.shape[0]/2),i]
        if pix_val==1:
            inx.append(i)

def get_cell_area(cell_mask):
    cell_area=np.sum(cell_mask==1)
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

    enable_logger()
    configData = loadconfig(input_dir_name)
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    colors = configData["COLORS"]
    molecule_type = ['/mrna']
    genes=['actn2-mature','gapdh-mature','actn2-immature']

    #compute quadrat
    with h5py.File(path.data_dir + input_dir_name + '/' + basic_file_name, "a") as file_handler, h5py.File(
            path.data_dir + input_dir_name + '/' + secondary_file_name, "a") as secondary_file_handler:
        vmr_values_sig = []
        for gene in genes:
            image_list=[]
            h_star_l = []
            [gene,timepoint] = gene.split("-")
            total_image=0
            nuc_dist = []
            cell_masks=[]
            image_counter=0
            image_cpt = 0
            for im in secondary_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                total_image += 1
            w = 40
            h = 2*total_image
            for im in secondary_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                image_number = im.split("_")[0]
                image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                spots_reduced = idsc.get_spots(secondary_file_handler, sec_image)
                nucleus_mask = idsc.get_nucleus_mask(secondary_file_handler, sec_image)
                cell_mask = idsc.get_cell_mask(secondary_file_handler, sec_image)
            vmr_values = []
            gs = gridspec.GridSpec(total_image * 6,1)

            for im in secondary_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                fig = plt.figure()
                fig.set_size_inches(w, 10)
                sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                image_number = im.split("_")[0]
                image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                spots_reduced = idsc.get_spots(secondary_file_handler, sec_image)
                z_lines = idsc.get_z_lines_masks(file_handler, image)
                sums = np.zeros((len(z_lines), 1))
                for n in range(len(z_lines)):
                    slice = z_lines[n]
                    sums[n] = slice.sum()
                z_level=sums.argmax()
                nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                cell_mask = idsc.get_cell_mask(secondary_file_handler, sec_image)
                z_lines_idx=reduce_z_line_mask(z_lines, spots_reduced)
                spots_reduced=spots_reduced[z_lines_idx[0] <= spots_reduced[:,2] ]
                spots_reduced=spots_reduced[spots_reduced[:,2] <= z_lines_idx[len(z_lines_idx)-1]]
                cell_area = get_cell_area(cell_mask)
                spot_surfacic_density = len(spots_reduced) / float(cell_area)
                cell_width=cell_mask.shape[1] -240
                band_n = 4.0
                quadrat_edge = cell_width / band_n
                grid_1d = np.zeros((int(band_n)))

                for spot in spots_reduced:
                    if spot[0] > 120 and spot[0]< cell_mask.shape[1]- 120:
                        x = int(np.floor((spot[0]-120)/ quadrat_edge))
                        grid_1d[x] += 1
                grid = []
                for val in grid_1d:
                    grid.append(val)
                grid_mat = np.matrix(grid).reshape((1, len(grid)))
                grid_mat/=quadrat_edge
                grid_mat/=spot_surfacic_density
                grid_list=list(np.array(grid_mat)[0])
                var = compute_squared_diff_sum(grid_list)
                vmr_values.append((var / np.mean(grid_list))-30.144)
                ax0 = plt.subplot()
                ax0.set_yticks([])
                major_ticks = np.arange(0,  int(band_n)+1, 1)
                ax0.tick_params(axis='both', which='major', labelsize=20)
                ax0.set_xticks(major_ticks)
                plt.imshow(grid_mat,  cmap='coolwarm',aspect=band_n/3)
                plt.ylim((0,0.4))
                image_cpt+=1
                plt.savefig(check_dir(path.analysis_dir+"/muscle_data/heatmap/"+"bucket_"+ str(int(band_n)))+"/mrna_distribution_" + gene + "_" + timepoint +"_image_"+str(image_cpt)+ ".png")
                plt.close()
            vmr_values_sig.append(np.sum(np.array(vmr_values)>0)/float(len(vmr_values)))
        figname="clustering_index.png"
        plot.bar_profile(vmr_values_sig, genes,figname)
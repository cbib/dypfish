#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt
from numpy import matlib
import math
import pandas as pd
import src.image_descriptors as idsc
import src.path as path
import src.plot as plot
from src.utils import enable_logger, check_dir, loadconfig

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def compute_minimal_distance(segment_summed):
    for i in range(15):
        if segment_summed[i] != 0:
            return i

def keep_cell_mask_spots(spots, cell_mask, z_lines):
    new_spots_list = []
    for spot in spots:
        if cell_mask[spot[1], spot[0]] == 1:
            new_spots_list.append(spot)
    return new_spots_list

def get_quantized_grid(q, Qx, Qy):
    tmp_x = np.matrix(np.arange(Qx))
    tmp_y = np.matrix(np.arange(Qy))
    qxs = matlib.repmat(tmp_x.transpose(), 1, Qx)
    qys = matlib.repmat(tmp_y, Qy, 1)
    qxs = np.kron(qxs, np.ones((q, q)))
    print(qxs)
    plt.imshow(qxs)
    plt.show()
    qys = np.kron(qys, np.ones((q, q)))
    return qxs, qys

def get_variance(test):
    print(test)
    N = len(test) - 1
    var = test - np.mean(test)
    tot = 0
    for elem in var:
        tot += math.pow(elem, 2)
    return (tot / N)

def compute_cell_mask_between_nucleus_centroid(cell_mask, nucleus_centroid, nuc_dist, cell_masks, nucs_dist):
    nucleus_centroid = np.sort(nucleus_centroid, axis=0)
    im_mask = []
    nucs = []
    for nuc_n in range(len(nucleus_centroid) - 1):
        cell_mask_copy = cell_mask.copy()
        nuc_dist.append(nucleus_centroid[nuc_n + 1][0] - nucleus_centroid[nuc_n][0])
        nucs.append(nucleus_centroid[nuc_n + 1][0] - nucleus_centroid[nuc_n][0])

        cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
        cell_mask_copy[:, nucleus_centroid[nuc_n + 1][0]::] = 0
        im_mask.append(cell_mask_copy)
    nucs_dist.append(nucs)
    cell_masks.append(im_mask)
    return nuc_dist, nucs_dist, cell_masks

def search_best_centroid(nucleus_centroid):
    x_nuc = []
    nucleus_centroid = np.sort(nucleus_centroid, axis=0)
    x_nuc.append(nucleus_centroid[0][0])
    x_nuc.append(nucleus_centroid[len(nucleus_centroid) - 1][0])
    return x_nuc

def show_descriptors(cell_mask, spots, nucleus_mask):
    plt.imshow(cell_mask)
    xs = spots[:, 0]
    ys = spots[:, 1]
    plt.grid(True)
    plt.scatter(xs, ys, color='white', marker=".", facecolors='none', linewidths=0.5)
    plt.show()

def compute_cell_mask_3d(cell_mask, z_lines):
    x_dim = cell_mask.shape[0]
    y_dim = cell_mask.shape[1]
    cell_mask_3d = np.zeros((x_dim, y_dim, len(z_lines)))
    for slice in range(len(z_lines)):
        cell_mask_3d[:, :, slice] = cell_mask
    return cell_mask_3d

def reject_outliers(data):
    index_cpt = 0
    indexes = []
    tmp_list = []
    for i in data:
        for j in i:
            tmp_list.append(j)
    u = np.mean(tmp_list)
    for i in data:
        print(i)
        for j in i:
            if u - 200 < j < u + 200:
                indexes.append(index_cpt)
        index_cpt += 1
    return indexes

# def compute_degree_of_clustering(spots_reduced, mask, z_lines):
#     cell_mask_3d = compute_cell_mask_3d(mask, z_lines)
#     h_star = helps.clustering_index_point_process(spots_reduced, cell_mask_3d)
#     d = np.array(h_star[h_star > 1] - 1).sum()
#     return d

def reduce_z_line_mask(z_lines, spots):
    cpt_z = 1
    z_lines_idx = []
    for z_line_mask in z_lines:
        spots_reduced = spots[spots[:, 2] == cpt_z]
        if len(spots_reduced) > 25:
            z_lines_idx.append(cpt_z)
        cpt_z += 1
    return z_lines_idx

def compute_zline_distance(file_handler, muscle_file_handler, molecule_type, genes, timepoints, z_line_spacing):
    all_median_profiles = []
    for gene in genes:
        for timepoint in timepoints:
            image_counter = 0
            total_profile = []
            print(molecule_type[0] + '/' + gene + '/' + timepoint)
            for im in muscle_file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                sec_image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                image_number = im.split("_")[0]
                image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + image_number
                print("computing " + image + " .....")
                spots = idsc.get_spots(muscle_file_handler, sec_image)
                z_lines = idsc.get_z_lines_masks(file_handler, image)
                image_spots_min_distance = np.zeros(len(spots))
                z_lines_idx = reduce_z_line_mask(z_lines, spots)
                spots_reduced = spots[z_lines_idx[0] <= spots[:, 2]]
                spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
                spots_count = 0
                for spot in spots_reduced:
                    z_line_mask = z_lines[int(spot[2])]
                    if z_line_mask[spot[1], spot[0]] == 0:
                        total_segment = np.zeros((360, z_line_spacing))
                        for degree in range(360):
                            z_line_segment = np.zeros(z_line_spacing)
                            line = np.zeros((z_line_spacing, 2))
                            angle = degree * 2 * math.pi / 360
                            x_slope, y_slope = math.sin(angle), math.cos(angle)
                            for point in range(z_line_spacing):
                                x = int(round(spot[0] + point * x_slope))
                                y = int(round(spot[1] + point * y_slope))
                                line[point, 0] = x
                                line[point, 1] = y
                                if (x >= 0 and x < z_line_mask.shape[1] and y >= 0 and y < z_line_mask.shape[0]):
                                    z_line_segment[point] = z_line_mask[y, x]
                                    total_segment[degree, point] = z_line_mask[y, x]
                                else:
                                    z_line_segment[point] = 0
                                    total_segment[degree, point] = 0
                        distance = compute_minimal_distance(np.sum(total_segment, axis=0))
                        image_spots_min_distance[spots_count] = distance
                    else:
                        image_spots_min_distance[spots_count] = 0
                    spots_count += 1
                z_line_distance_profile = np.zeros(z_line_spacing)
                for i in range(z_line_spacing):
                    z_line_distance_profile[i] = float(len(np.where(image_spots_min_distance == i)[0])) / float(
                        len(image_spots_min_distance))
                total_profile.append(z_line_distance_profile)
                image_counter += 1
            total_profile = np.array(total_profile).reshape((image_counter, z_line_spacing))
            all_median_profiles.append(np.median(total_profile, axis=0))
    return all_median_profiles

if __name__ == "__main__":
    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()

    configData = loadconfig(input_dir_name)

    # timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    # timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    # timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    # timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    # mime_type = configData["PNG_IMAGES_MIME_TYPE"]
    # molecule_types = configData["MOLECULE_TYPES"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    z_line_spacing = configData["Z_LINE_SPACING"]
    image_height = configData["IMAGE_HEIGHT"]

    molecule_type = ['/mrna']
    # genes = ['actn2', 'gapdh']
    # timepoints = ['mature']
    # colors = ['#0A3950', '#1E95BB', '#A1BA6D']
    # with h5py.File(path.data_dir + input_dir_name + '/' + basic_file_name, "a") as file_handler,h5py.File(path.data_dir + input_dir_name + '/' + secondary_file_name, "a") as muscle_file_handler:
    #
    #     all_median_profiles = compute_zline_distance(file_handler, muscle_file_handler, molecule_type, genes, timepoints, z_line_spacing)
    #     df = pd.DataFrame(all_median_profiles)
    #     df.to_csv(check_dir(path.analysis_dir + 'muscle_data/dataframe/') + "all_median_profiles_by_slice_mature_remove_bad_slice.csv")


    genes = ['actn2']
    timepoints = ['immature']
    with h5py.File(path.data_dir + input_dir_name + '/' + basic_file_name, "a") as file_handler, h5py.File(path.data_dir + input_dir_name + '/' + secondary_file_name,"a") as muscle_file_handler:
        all_median_profiles = compute_zline_distance(file_handler, muscle_file_handler, molecule_type, genes,timepoints,z_line_spacing)
        df = pd.DataFrame(all_median_profiles)
        df.to_csv(check_dir(path.analysis_dir + 'muscle_data/dataframe/') + "all_median_profiles_by_slice_immature_remove_bad_slice.csv")
    all_median_profiles = []
    plot_name = check_dir(path.analysis_dir + 'muscle_data/figures/') + 'z_line_distance' + str(z_line_spacing) + 'contours_mature_remove_bad_slice.png'
    figure_title = 'z line spots distance profile'
    genes = ["actn2 mature", "gapdh mature", "actn2 immature"]
    df = pd.read_csv(path.analysis_dir + "muscle_data/dataframe/all_median_profiles_by_slice_mature_remove_bad_slice.csv",index_col=0)
    df_im = pd.read_csv(path.analysis_dir + "muscle_data/dataframe/all_median_profiles_by_slice_immature_remove_bad_slice.csv",index_col=0)
    all_median_profiles.append(df.ix[0].values)
    all_median_profiles.append(df.ix[1].values)
    all_median_profiles.append(df_im.ix[0].values)
    plot.profile(all_median_profiles, genes, z_line_spacing, plot_name, figure_title, colors, True)

#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import math

import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
from skimage.draw import circle
from skimage import measure
import pandas as pd

import os
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.constants as cst

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


def compute_eight_quadrant_max(file_handler, image, degree_max):
    import matplotlib.pyplot as plt
    import sys
    print(image)
    mtoc_position = idsc.get_mtoc_position(file_handler, image)
    height_map = idsc.get_height_map(file_handler, image)
    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
    cell_mask = idsc.get_cell_mask(file_handler, image)
    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
    height_map = height_map.astype(float)
    height_map[cell_mask == 0] = 0
    height_map[nucleus_mask == 1] = 0
    # print(degree_max)

    # the quadrant of MTOC is defined by two lines 45 degrees to the right
    right_point = helps.rotate_point(nucleus_centroid, mtoc_position, int(degree_max))

    #print(right_point)
    s = helps.slope_from_points(nucleus_centroid, right_point)
    #print(s)

    # plt.imshow(cell_mask)
    # plt.scatter(mtoc_position[0], mtoc_position[1], color='black', marker="d", linewidths=3)
    # plt.scatter(right_point[0], right_point[1], color='white', marker="d", linewidths=3)
    # plt.show()
    corr = np.arctan(s)  # angle wrt to x axis
    xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
    rotated_xx, rotated_yy = helps.rotate_meshgrid(xx, yy, -corr)
    sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (math.pi)))
    sliceno = sliceno.astype(int)

    # plt.imshow(sliceno)
    # plt.scatter(mtoc_position[0], mtoc_position[1], color='black', marker="d", linewidths=3)
    # plt.show()
    quadrant_mask = sliceno + cell_mask
    quadrant_mask[quadrant_mask == 9] = 8
    quadrant_mask[cell_mask == 0] = 0

    # plt.imshow(quadrant_mask)
    # plt.scatter(mtoc_position[0], mtoc_position[1], color='black', marker="d", linewidths=3)
    #
    # plt.show()
    mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]

    return quadrant_mask, mtoc_quad


# Calculates the quadrant mask for the MTOC
def compute_quadrant_max(file_handler, image,degree_max):
    print(image)
    mtoc_position = idsc.get_mtoc_position(file_handler, image)
    height_map = idsc.get_height_map(file_handler, image)
    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
    cell_mask = idsc.get_cell_mask(file_handler, image)
    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
    spots = idsc.get_spots(file_handler, image)
    spot_by_quad = np.zeros((90, 4, 2))
    spot_by_quad_mean = np.zeros((90, 2, 2))
    height_map = height_map.astype(float)
    height_map[cell_mask == 0] = 0
    height_map[nucleus_mask == 1] = 0
    print(degree_max)


    # the quadrant of MTOC is defined by two lines 45 degrees to the right
    right_point = helps.rotate_point(nucleus_centroid, mtoc_position, int(degree_max)-1)
    s = helps.slope_from_points(nucleus_centroid, right_point)
    corr = np.arctan(s) # angle wrt to x axis
    xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
    rotated_xx, rotated_yy = helps.rotate_meshgrid(xx, yy, -corr)
    sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
    sliceno = sliceno.astype(int)
    quadrant_mask = sliceno + cell_mask
    quadrant_mask[quadrant_mask == 5]=4
    quadrant_mask[cell_mask == 0] = 0
    mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
    #mtoc_quads.append(mtoc_quad)



    # # show quadrant mask with mtoc and spots
    # xs = spots[:, 0]
    # ys = spots[:, 1]
    # fig, ax = plt.subplots()
    # # plt.scatter(xs, ys, color='black', marker="*", facecolors='none', linewidths=0.5)
    # contours = measure.find_contours(nucleus_mask, 0.8)
    # for n, contour in enumerate(contours):
    #     plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    # rr_m, cc_m = circle(mtoc_position[1], mtoc_position[0], 5, quadrant_mask.shape)
    # quadrant_mask[rr_m, cc_m] = 0
    # rr, cc = circle(nucleus_centroid[1], nucleus_centroid[0], 5, quadrant_mask.shape)
    # quadrant_mask[rr, cc] = 0
    # ax.imshow(quadrant_mask)
    # plt.show()

    return quadrant_mask,mtoc_quad

def reindex_quadrant_mask(quad_mask,mtoc_quad):
    print(mtoc_quad)
    if mtoc_quad==1:
        return quad_mask

    elif mtoc_quad==2:

        quad_mask[quad_mask == 1] = 5
        quad_mask[quad_mask == 2] = 1
        quad_mask[quad_mask == 3] = 2
        quad_mask[quad_mask == 4] = 3
        quad_mask[quad_mask == 5] = 4
        #plt.imshow(quad_mask)
        #plt.show()
    elif mtoc_quad==3:
        quad_mask[quad_mask == 3] = 5
        quad_mask[quad_mask == 4] = 6
        quad_mask[quad_mask == 1] = 7
        quad_mask[quad_mask == 2] = 8
        quad_mask[quad_mask == 5] = 1
        quad_mask[quad_mask == 6] = 2
        quad_mask[quad_mask == 7] = 3
        quad_mask[quad_mask == 8] = 4

        #plt.imshow(quad_mask)
        #plt.show()
    else:
        quad_mask[quad_mask == 4] = 5
        quad_mask[quad_mask == 1] = 6
        quad_mask[quad_mask == 2] = 7
        quad_mask[quad_mask == 3] = 8
        quad_mask[quad_mask == 5] = 1
        quad_mask[quad_mask == 6] = 2
        quad_mask[quad_mask == 7] = 3
        quad_mask[quad_mask == 8] = 4
        #plt.imshow(quad_mask)
        #plt.show()
    return quad_mask

def map_index(x,mtoc_quad):
    print(mtoc_quad)

    #for i in range(1, 9):
    if x >= mtoc_quad:
        new_index=(x-mtoc_quad+1)
    else:
        new_index=(x+8-mtoc_quad+1)
    return new_index




def reindex_quadrant_mask3(quad_mask, mtoc_quad):
     print(mtoc_quad)
     import pandas as pd
     df=pd.DataFrame(quad_mask)
     df=df.applymap(lambda x: x-mtoc_quad+1 if x >= mtoc_quad else (x+8-mtoc_quad+1 if x >0 else 0))
     quad_mask=np.array(df)

     return quad_mask



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
    # basic_md5_path=path.analysis_data_dir + 'basic.md5'
    # secondary_md5_path=path.analysis_data_dir + 'secondary.md5'
    # mtoc_md5_path = path.analysis_data_dir + 'mtoc.md5'
    # pwd_path = path.analysis_data_dir + 'password.txt'
    path_data = path.raw_data_dir

    # f = open(pwd_path, 'r')
    # loginpwd = f.readline()
    # login = loginpwd.split(":")[0]
    # pwd = loginpwd.split(":")[1]
    # helps.check_data(basic_file_path, 'basic', login, pwd)
    # helps.check_data(mtoc_file_path, 'mtoc', login, pwd)
    # q = 32; # pixels per voxel width
    # Q = 512 / q; # voxels per image width
    # qxs,qys=helps.get_quantized_grid(q,Q)

    df = pd.read_csv('global_mtoc_file_all_mrna_nocodazole')
    df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()

    degree_max_mrna={}
    for gene, line in df_sorted.groupby(['Gene','timepoint']):
        key=gene[0]+"_"+gene[1]
        degree_max_mrna[key]=line['Unnamed: 0'].values

    df = pd.read_csv('global_mtoc_file_all_protein_nocodazole')
    df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()

    degree_max_protein = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_protein[key] = line['Unnamed: 0'].values

    # Compute global TIS
    with h5py.File(basic_file_path, "r") as file_handler,h5py.File(secondary_file_path, "r") as second_file_handler,h5py.File(mtoc_file_path, "r") as mtoc_file_handler:
        molecule_type=['/mrna']
        #mrnas = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        mrnas = ["arhgdia_nocodazole","pard3_nocodazole"]
        #mrnas = ["pard3"]
        mrna_timepoints = ["3h","5h"]
        #h_mrna = np.empty((4, 4), dtype=object)
        gene_count=0
        for mrna in mrnas:
            print(mrna)
            time_count = 0
            for timepoint in mrna_timepoints:
                key=mrna+"_"+timepoint
                image_count = 0

                h_array = np.zeros((len(degree_max_mrna[key]), cst.STRIPE_NUM*8))
                for i in range(len(degree_max_mrna[key])):
                    print(degree_max_mrna[key][i])
                    image = degree_max_mrna[key][i].split("_")[0]
                    degree = degree_max_mrna[key][i].split("_")[1]
                    image = "/mrna/" + mrna + "/" + timepoint + "/" + image

                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask = reindex_quadrant_mask3(quad_mask, mtoc_quad)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    height_map = idsc.get_height_map(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    mtoc_position = idsc.get_mtoc_position(file_handler, image)
                    #quad_mask = idsc.set_quadrants(file_handler, image)
                    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
                    #mtoc_quad = quad_mask[mtoc_position[1], mtoc_position[0]]
                    #print(mtoc_quad)
                    #sys.exit()
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    #h = np.zeros((Q, Q))

                    cell_mask_dist_map[(cell_mask == 1) & (cell_mask_dist_map == 0)] = 1
                    cell_mask_dist_map[(nucleus_mask == 1)] = 0

                    for i in range(len(spots)):
                        if idsc.is_in_cytoplasm(file_handler, image, [spots[i, 1], spots[i, 0]]):
                            # round_value=np.round(cell_mask_dist_map[spots[i, 1], spots[i, 0]])

                            quad=quad_mask[spots[i, 1], spots[i, 0]]

                            value = cell_mask_dist_map[spots[i, 1], spots[i, 0]]
                            # plt.imshow(cell_mask_dist_map)
                            # plt.show()
                            if value == 100:
                                value = 99
                            dist = value
                            value = np.floor(value / (100.0 / cst.STRIPE_NUM))
                            value=(int(value)*8)+int(quad-1)

                            # if value == 40:
                            #     value = 39
                            #print(np.sum(cell_mask[(cell_mask_dist_map==dist) & (quad_mask==quad)]))* math.pow((1 /  cst.SIZE_COEFFICIENT), 2)
                            h_array[image_count, int(value)] += (1.0 / len(spots))/float(np.sum(cell_mask[(cell_mask_dist_map==dist) & (quad_mask==quad)])* math.pow((1 / cst.SIZE_COEFFICIENT), 2))

                    image_count += 1
                mrna_tp_df = pd.DataFrame(h_array)
                mrna_tp_df.to_csv(path.analysis_dir+"analysis_nocodazole/df/"+mrna + '_' + timepoint + "_mrna.csv")
                # h_mrna[gene_count,time_count]=np.mean(h_array,axis=0)
                time_count += 1

            gene_count += 1

        proteins = ["arhgdia_nocodazole","pard3_nocodazole"]

        #proteins = ["pard3"]
        prot_timepoints = ["3h", "5h"]
        #h_prot = np.empty((len(mrnas), len(prot_timepoints)), dtype=object)
        gene_count = 0
        for protein in proteins:
            print(protein)
            time_count = 0
            for timepoint in prot_timepoints:
                key = protein + "_" + timepoint

                image_count = 0
                #h_array = np.zeros((len(image_list), 10))
                h_array = np.zeros((len(degree_max_protein[key]), 8*cst.STRIPE_NUM))

                for i in range(len(degree_max_protein[key])):
                    print(degree_max_protein[key][i])
                    image = degree_max_protein[key][i].split("_")[0]
                    degree = degree_max_protein[key][i].split("_")[1]
                    image = "/protein/" + protein + "/" + timepoint + "/" + image

                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask = reindex_quadrant_mask3(quad_mask, mtoc_quad)


                    image_number = image.split("/")[4]
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    height_map = idsc.get_height_map(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    mtoc_position = idsc.get_mtoc_position(file_handler, image)
                    #quad_mask = idsc.set_quadrants(file_handler, image)
                    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
                   # mtoc_quad = quad_mask[mtoc_position[1], mtoc_position[0]]
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    IF = helps.get_IF_image_z_summed(protein, 'protein', timepoint, image_number,path_data)
                    count=0
                    IF[(cell_mask == 0)] = 0
                    IF[(nucleus_mask == 1)] = 0

                    cell_mask_dist_map[(cell_mask==1) & (cell_mask_dist_map==0)]=1
                    #cell_mask_dist_map[(nucleus_mask == 1)] = 0
                    # if "/protein/pard3/2h/9" in image:
                    #     plt.imshow(cell_mask_dist_map)
                    #     plt.show()
                    #     plt.imshow(quad_mask)
                    #     plt.show()

                    for i in range(1,99):
                        #print(np.sum(IF))
                        #plt.imshow(IF)
                        #plt.show()

                        for j in range (1,9):

                            value = int(np.floor(i / (100.0 / cst.STRIPE_NUM)))
                            value = (int(value) * 8) + int(j - 1)
                            #print(float(np.sum(IF[(cell_mask_dist_map>i) & (cell_mask_dist_map < i+9)]))/float(np.sum(IF)))
                            # IF[(cell_mask_dist_map == i) & (quad_mask == j)]=10000000000
                            # plt.imshow(IF)
                            # plt.show()

                            # print(i)
                            # print(j)
                            #
                            # print(np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == 1)]))
                            # print(np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == 2)]))
                            # print(np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == 3)]))
                            # print(np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == 4)]))
                            # print(np.sum(cell_mask[(cell_mask_dist_map==i) & (quad_mask==j)]))
                            if np.sum(cell_mask[(cell_mask_dist_map==i) & (quad_mask==j)])==0:
                                h_array[image_count, value] =0.0
                            else:
                                h_array[image_count,value]=(float(np.sum(IF[(cell_mask_dist_map == i) & (quad_mask==j)]))/float(np.sum(IF)))/np.sum(cell_mask[(cell_mask_dist_map==i) & (quad_mask==j)])

                        count+=1
                    image_count += 1
                    # print(spots_peripheral_distance/len(spots))
                    # h_array.append(spots_peripheral_distance/len(spots))
                #print(h_array)`
                #print(np.mean(h_array,axis=0))
                prot_tp_df = pd.DataFrame(h_array)
                prot_tp_df.to_csv(path.analysis_dir+"analysis_nocodazole/df/"+protein + '_' + timepoint + "_protein.csv")
                #h_prot[gene_count, time_count] = np.mean(h_array,axis=0)
                time_count += 1

            gene_count += 1
        #print(h_prot)







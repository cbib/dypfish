#!/usr/bin/python
# encoding: UTF-8
# author: benjamin Dartigues
import numpy as np
import h5py
import math
import pandas as pd
import sys
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
import src.constants as cst
import matplotlib.pyplot as plt

from src.utils import enable_logger, check_dir


def reindex_quadrant_mask(quad_mask, mtoc_quad):
    import pandas as pd
    df = pd.DataFrame(quad_mask)
    df = df.applymap(lambda x: x - mtoc_quad + 1 if x >= mtoc_quad else (x + 8 - mtoc_quad + 1 if x > 0 else 0))
    quad_mask = np.array(df)

    return quad_mask

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

def main():
    # Required descriptors: spots, IF, cell mask an height_map
    # you need to run before MTOC analysis script called search enriched quad

    enable_logger()

    try:
        df = pd.read_csv(path.analysis_dir + 'analysis_MTOC/dataframe/' + 'periph_global_mtoc_file_all_mrna.csv')
    except IOError:
        print "Couldn't load file : ", path.analysis_dir + 'analysis_MTOC/dataframe/' + 'periph_global_mtoc_file_all_mrna.csv'
        print "Maybe MTOC analysis hasn't been launched prior to this one"
        exit(1)
    # TODO check sort type
    df_sorted = df.sort_values('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],
                                                               as_index=False).first()  #  TODO check sort type
    degree_max_mrna = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_mrna[key] = line['Unnamed: 0'].values

    df = pd.read_csv(path.analysis_dir + 'analysis_MTOC/dataframe/' + 'periph_global_mtoc_file_all_protein.csv')
    df_sorted = df.sort_values('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()

    degree_max_protein = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_protein[key] = line['Unnamed: 0'].values
    limit=30
    # Compute global TIS
    with h5py.File(path.basic_file_path, "r") as file_handler, \
            h5py.File(path.secondary_file_path, "r") as second_file_handler:
        mrnas = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        mrna_timepoints = ["2h", "3h", "4h", "5h"]
        gene_count = 0
        for mrna in mrnas:
            time_count = 0
            for timepoint in mrna_timepoints:
                key = mrna + "_" + timepoint
                image_count = 0
                h_array = np.zeros((len(degree_max_mrna[key]),cst.STRIPE_NUM*8 ))
                for i in range(len(degree_max_mrna[key])):
                    image = degree_max_mrna[key][i].split("_")[0]
                    degree = degree_max_mrna[key][i].split("_")[1]
                    image = "/mrna/" + mrna + "/" + timepoint + "/" + image
                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask = reindex_quadrant_mask(quad_mask, mtoc_quad)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    cell_mask_dist_map[(cell_mask == 1) & (cell_mask_dist_map == 0)] = 1
                    cell_mask_dist_map[(nucleus_mask == 1)] = 0
                    peripheral_binary_mask = (cell_mask_dist_map > 0) & (cell_mask_dist_map <= limit).astype(int)
                    #plt.imshow(peripheral_binary_mask)
                    #plt.show()
                    spots_in_periph=0
                    for i in range(len(spots)):
                        #if idsc.is_in_cytoplasm(file_handler, image, [spots[i, 1], spots[i, 0]]):
                        if peripheral_binary_mask[spots[i, 1], spots[i, 0]]:
                            spots_in_periph+=1

                    for i in range(len(spots)):
                        #if idsc.is_in_cytoplasm(file_handler, image, [spots[i, 1], spots[i, 0]]):
                        if peripheral_binary_mask[spots[i, 1], spots[i, 0]]:

                            quad = quad_mask[spots[i, 1], spots[i, 0]]
                            value = cell_mask_dist_map[spots[i, 1], spots[i, 0]]

                            if value == limit:
                                value = limit-1
                            slice_area = np.floor(value / (float(limit)/cst.STRIPE_NUM))
                            value = (int(slice_area) * 8) + int(quad)
                            print("low_limit",(slice_area*(limit/cst.STRIPE_NUM))+1)
                            print("high limit",(slice_area *(limit/cst.STRIPE_NUM))+(limit/cst.STRIPE_NUM))
                            spots_relative=(1.0 / spots_in_periph)
                            surface_relative=float(np.sum(cell_mask[((cell_mask_dist_map >= (slice_area*(limit/cst.STRIPE_NUM))+1) & (cell_mask_dist_map < (slice_area *(limit/cst.STRIPE_NUM))+(limit/cst.STRIPE_NUM))) & (quad_mask == quad)]) * math.pow((1 / cst.SIZE_COEFFICIENT), 2))/float(np.sum(cell_mask[cell_mask==1])* math.pow((1 / cst.SIZE_COEFFICIENT), 2))
                            #print(spots_relative)
                            #print(surface_relative)
                            if surface_relative==0.0:
                                h_array[image_count, int(value) - 1] +=0
                            else:
                                h_array[image_count, int(value)-1] += (spots_relative/surface_relative)

                            #h_array[image_count, int(value)-1] += (spots_relative/float(np.sum(cell_mask[((cell_mask_dist_map >= (slice_area*(limit/cst.STRIPE_NUM))+1) & (cell_mask_dist_map < (slice_area *(limit/cst.STRIPE_NUM))+(limit/cst.STRIPE_NUM))) & (quad_mask == quad)]) * math.pow((1 / cst.SIZE_COEFFICIENT), 2)))

                    image_count += 1
                mrna_tp_df = pd.DataFrame(h_array)
                mrna_tp_df.to_csv(check_dir(path.analysis_dir + "analysis_temporal_interactions/dataframe/") +
                                  mrna + '_' + timepoint+ '_' +str(limit) +"_mrna.csv")
                time_count += 1
            gene_count += 1

        proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        prot_timepoints = ["2h", "3h", "5h", "7h"]
        gene_count = 0
        for protein in proteins:
            print(protein)
            time_count = 0
            for timepoint in prot_timepoints:
                key = protein + "_" + timepoint

                image_count = 0
                h_array = np.zeros((len(degree_max_protein[key]),cst.STRIPE_NUM*8 ))

                for i in range(len(degree_max_protein[key])):
                    #print(degree_max_protein[key][i])
                    image = degree_max_protein[key][i].split("_")[0]
                    degree = degree_max_protein[key][i].split("_")[1]
                    image = "/protein/" + protein + "/" + timepoint + "/" + image

                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask = reindex_quadrant_mask(quad_mask, mtoc_quad)

                    image_number = image.split("/")[4]
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    IF = helps.get_IF_image_z_summed(protein, 'protein', timepoint, image_number, path.path_data)
                    count = 0
                    IF[(cell_mask == 0)] = 0
                    IF[(nucleus_mask == 1)] = 0
                    cell_mask_dist_map[(cell_mask == 1) & (cell_mask_dist_map == 0)] = 1
                    cpt_stripes=0

                    for i in range((limit/cst.STRIPE_NUM), limit+1, (limit/cst.STRIPE_NUM)):
                        print("i",i)
                        for j in range(1, 9):
                            if np.sum(cell_mask[((cell_mask_dist_map >= 1+(i-(limit/cst.STRIPE_NUM))) & (cell_mask_dist_map < i)) & (quad_mask == j)]) == 0:
                                h_array[image_count, cpt_stripes] = 0.0
                            else:
                                IF_relative=float(np.sum(IF[((cell_mask_dist_map >= 1+(i-(limit/cst.STRIPE_NUM))) & (cell_mask_dist_map < i)) & (quad_mask == j)])) / np.sum(IF[cell_mask==1])/np.sum(IF[cell_mask==1])
                                surface_relative=np.sum(cell_mask[((cell_mask_dist_map >= 1+(i-(limit/cst.STRIPE_NUM))) & (cell_mask_dist_map < i)) & (quad_mask == j)]) * math.pow((1 / cst.SIZE_COEFFICIENT), 2) / np.sum(cell_mask[cell_mask==1])* math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                                h_array[image_count, cpt_stripes] =IF_relative/surface_relative
                            cpt_stripes+=1
                        count += 1
                    image_count += 1
                prot_tp_df = pd.DataFrame(h_array)
                prot_tp_df.to_csv(
                    path.analysis_dir + "analysis_temporal_interactions/dataframe/" + protein + '_' + timepoint+ '_' +str(limit) +"_protein.csv")
                time_count += 1

            gene_count += 1


# Calculates the quadrant mask for the MTOC
def compute_quadrant_max(file_handler, image, degree_max):
    print(image)
    mtoc_position = idsc.get_mtoc_position(file_handler, image)
    height_map = idsc.get_height_map(file_handler, image)
    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
    cell_mask = idsc.get_cell_mask(file_handler, image)
    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
    height_map = height_map.astype(float)
    height_map[cell_mask == 0] = 0
    height_map[nucleus_mask == 1] = 0
    print(degree_max)

    # the quadrant of MTOC is defined by two lines 45 degrees to the right
    right_point = helps.rotate_point(nucleus_centroid, mtoc_position, int(degree_max) - 1)
    s = helps.slope_from_points(nucleus_centroid, right_point)
    corr = np.arctan(s)  # angle wrt to x axis
    xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
    rotated_xx, rotated_yy = helps.rotate_meshgrid(xx, yy, -corr)
    sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
    sliceno = sliceno.astype(int)
    quadrant_mask = sliceno + cell_mask
    quadrant_mask[quadrant_mask == 5] = 4
    quadrant_mask[cell_mask == 0] = 0
    mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
    return quadrant_mask, mtoc_quad


def reindex_quadrant_mask(quad_mask, mtoc_quad):
    print(mtoc_quad)
    if mtoc_quad == 1:
        return quad_mask

    elif mtoc_quad == 2:

        quad_mask[quad_mask == 1] = 5
        quad_mask[quad_mask == 2] = 1
        quad_mask[quad_mask == 3] = 2
        quad_mask[quad_mask == 4] = 3
        quad_mask[quad_mask == 5] = 4

    elif mtoc_quad == 3:
        quad_mask[quad_mask == 3] = 5
        quad_mask[quad_mask == 4] = 6
        quad_mask[quad_mask == 1] = 7
        quad_mask[quad_mask == 2] = 8
        quad_mask[quad_mask == 5] = 1
        quad_mask[quad_mask == 6] = 2
        quad_mask[quad_mask == 7] = 3
        quad_mask[quad_mask == 8] = 4

    else:
        quad_mask[quad_mask == 4] = 5
        quad_mask[quad_mask == 1] = 6
        quad_mask[quad_mask == 2] = 7
        quad_mask[quad_mask == 3] = 8
        quad_mask[quad_mask == 5] = 1
        quad_mask[quad_mask == 6] = 2
        quad_mask[quad_mask == 7] = 3
        quad_mask[quad_mask == 8] = 4

    return quad_mask


if __name__ == "__main__":
    main()

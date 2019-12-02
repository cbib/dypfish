#!/usr/bin/python
# encoding: UTF-8
# author: benjamin Dartigues


import numpy as np
import h5py
import math
import pandas as pd
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
import src.constants as cst
from src.utils import enable_logger, check_dir, loadconfig

def compute_eight_quadrant_max(file_handler, image, degree_max):
    print(image)
    mtoc_position = idsc.get_mtoc_position(file_handler, image)
    height_map = idsc.get_height_map(file_handler, image)
    nucleus_centroid = idsc.get_nucleus_centroid(file_handler, image)
    cell_mask = idsc.get_cell_mask(file_handler, image)
    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
    height_map = height_map.astype(float)
    height_map[cell_mask == 0] = 0
    height_map[nucleus_mask == 1] = 0

    # the quadrant of MTOC is defined by two lines 45 degrees to the right
    right_point = helps.rotate_point(nucleus_centroid, mtoc_position, int(degree_max))
    s = helps.slope_from_points(nucleus_centroid, right_point)
    corr = np.arctan(s)  # angle wrt to x axis
    xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
    rotated_xx, rotated_yy = helps.rotate_meshgrid(xx, yy, -corr)
    sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (math.pi)))
    sliceno = sliceno.astype(int)
    quadrant_mask = sliceno + cell_mask
    quadrant_mask[quadrant_mask == 9] = 8
    quadrant_mask[cell_mask == 0] = 0
    mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]

    return quadrant_mask, mtoc_quad


def map_index(x,mtoc_quad):
    print(mtoc_quad)
    if x >= mtoc_quad:
        new_index=(x-mtoc_quad+1)
    else:
        new_index=(x+8-mtoc_quad+1)
    return new_index

def reindex_quadrant_mask(quad_mask, mtoc_quad):
     import pandas as pd
     df=pd.DataFrame(quad_mask)
     df=df.applymap(lambda x: x-mtoc_quad+1 if x >= mtoc_quad else (x+8-mtoc_quad+1 if x >0 else 0))
     quad_mask=np.array(df)
     return quad_mask

def main():
    enable_logger()
    try:
        df = pd.read_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_mrna.csv')
    except IOError:
        print "Couldn't load file : ", path.analysis_dir + 'analysis_MTOC/dataframe/' + 'global_mtoc_file_all_mrna.csv'
        print "Maybe MTOC analysis hasn't been launched prior to this one"
        exit(1)
    df_sorted = df.sort_values('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],
                                                               as_index=False).first()
    degree_max_mrna = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_mrna[key] = line['Unnamed: 0'].values
    df = pd.read_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_protein.csv')
    df_sorted = df.sort_values('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
    degree_max_protein = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_protein[key] = line['Unnamed: 0'].values

    configData = loadconfig("original")
    genes = configData["GENES"][0:4]
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    stripe_n= configData["STRIPE_NUM"]

    # Compute global TIS
    with h5py.File(path.basic_file_path, "r") as file_handler, \
            h5py.File(path.secondary_file_path, "r") as second_file_handler:
        gene_count = 0
        for mrna in genes:
            time_count = 0
            for timepoint in timepoints_mrna:
                key = mrna + "_" + timepoint
                image_count = 0
                h_array = np.zeros((len(degree_max_mrna[key]),stripe_n*8 ))
                for i in range(len(degree_max_mrna[key])):
                    image = degree_max_mrna[key][i].split("_")[0]
                    degree = degree_max_mrna[key][i].split("_")[1]
                    image = "/mrna/" + mrna + "/" + timepoint + "/" + image
                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask=reindex_quadrant_mask(quad_mask, mtoc_quad)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    cell_mask_dist_map[(cell_mask == 1) & (cell_mask_dist_map == 0)] = 1
                    cell_mask_dist_map[(nucleus_mask == 1)] = 0

                    for i in range(len(spots)):
                        if idsc.is_in_cytoplasm(file_handler, image, [spots[i, 1], spots[i, 0]]):
                            quad = quad_mask[spots[i, 1], spots[i, 0]]
                            value = cell_mask_dist_map[spots[i, 1], spots[i, 0]]
                            if value == 100:
                                value = 99
                            dist = value
                            slice_area = np.floor(value / (100.0/stripe_n))
                            value = (int(slice_area) * 8) + int(quad-1)
                            h_array[image_count, int(value)] += (1.0 / len(spots)) / float(np.sum(cell_mask[(cell_mask_dist_map == dist) & (quad_mask == quad)]) * math.pow((1 / cst.SIZE_COEFFICIENT), 2))
                    image_count += 1
                mrna_tp_df = pd.DataFrame(h_array)
                mrna_tp_df.to_csv(check_dir(path.analysis_dir + "analysis_temporal_interactions/dataframe/") +
                                  mrna + '_' + timepoint + "_mrna.csv")
                time_count += 1
            gene_count += 1
        gene_count = 0
        for protein in proteins:
            print(protein)
            time_count = 0
            for timepoint in timepoints_protein:
                key = protein + "_" + timepoint
                image_count = 0
                h_array = np.zeros((len(degree_max_protein[key]), 8* stripe_n))
                for i in range(len(degree_max_protein[key])):
                    image = degree_max_protein[key][i].split("_")[0]
                    degree = degree_max_protein[key][i].split("_")[1]
                    image = "/protein/" + protein + "/" + timepoint + "/" + image
                    quad_mask, mtoc_quad = compute_eight_quadrant_max(file_handler, image, degree)
                    quad_mask=reindex_quadrant_mask(quad_mask, mtoc_quad)
                    image_number = image.split("/")[4]
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    cell_mask_dist_map = idsc.get_cell_mask_distance_map(second_file_handler, image)
                    IF = helps.get_IF_image_z_summed(protein, 'protein', timepoint, image_number, path.path_data)
                    count = 0
                    IF[(cell_mask == 0)] = 0
                    IF[(nucleus_mask == 1)] = 0
                    cell_mask_dist_map[(cell_mask == 1) & (cell_mask_dist_map == 0)] = 1
                    for i in range(1, 99):
                        for j in range(1, 9):
                            value = int(np.floor(i / (100.0/stripe_n)))
                            value = (int(value) * 8) + int(j - 1)
                            if np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == j)]) == 0:
                                h_array[image_count, value] = 0.0
                            else:
                                h_array[image_count, value] = (float(
                                    np.sum(IF[(cell_mask_dist_map == i) & (quad_mask == j)])) / float(
                                    np.sum(IF))) / np.sum(cell_mask[(cell_mask_dist_map == i) & (quad_mask == j)])
                        count += 1
                    image_count += 1
                prot_tp_df = pd.DataFrame(h_array)
                prot_tp_df.to_csv(
                    path.analysis_dir + "analysis_temporal_interactions/dataframe/" + protein + '_' + timepoint + "_protein.csv")
                time_count += 1
            gene_count += 1

if __name__ == "__main__":
    main()
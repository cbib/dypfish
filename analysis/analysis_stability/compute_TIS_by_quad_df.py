#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import h5py
import math
import pandas as pd
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
import src.constants as cst
from src.utils import enable_logger, check_dir


def main():
    # Required descriptors: spots, IF, cell mask an height_map
    # you need to run before MTOC analysis script called search enriched quad
    enable_logger()

    try:
        df = pd.read_csv(path.analysis_dir + 'analysis_stability/dataframe/' + 'global_mtoc_file_all_mrna.csv')
    except IOError:
        print "Couldn't load file : ", path.analysis_dir + 'analysis_stability/dataframe/' + 'global_mtoc_file_all_mrna.csv'
        print "Maybe MTOC analysis hasn't been launched prior to this one"
        exit(1)
    df_sorted = df.sort_values('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],
                                                               as_index=False).first()
    degree_max_mrna = {}
    for gene, line in df_sorted.groupby(['Gene', 'timepoint']):
        key = gene[0] + "_" + gene[1]
        degree_max_mrna[key] = line['Unnamed: 0'].values

    # Compute global TIS
    with h5py.File(path.basic_file_path, "r") as file_handler, \
            h5py.File(path.secondary_file_path, "r") as second_file_handler:
        mrnas = ["arhgdia", "arhgdia_cultured"]
        mrna_timepoints = ["3h"]
        gene_count = 0
        for mrna in mrnas:
            time_count = 0
            for timepoint in mrna_timepoints:
                key = mrna + "_" + timepoint
                image_count = 0
                h_array = np.zeros((len(degree_max_mrna[key]), 40))
                for i in range(len(degree_max_mrna[key])):
                    image = degree_max_mrna[key][i].split("_")[0]
                    degree = degree_max_mrna[key][i].split("_")[1]
                    image = "/mrna/" + mrna + "/" + timepoint + "/" + image
                    quad_mask, mtoc_quad = compute_quadrant_max(file_handler, image, degree)
                    reindex_quadrant_mask(quad_mask, mtoc_quad)
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
                            value = np.floor(value / 10.0)
                            value = (int(value) * 4) + int(quad - 1)
                            h_array[image_count, int(value)] += (1.0 / len(spots)) / float(
                                np.sum(cell_mask[(cell_mask_dist_map == dist) & (quad_mask == quad)]) * math.pow(
                                    (1 / cst.SIZE_COEFFICIENT), 2))
                    image_count += 1
                mrna_tp_df = pd.DataFrame(h_array)
                mrna_tp_df.to_csv(check_dir(path.analysis_dir + "analysis_stability/dataframe/") +
                                  mrna + '_' + timepoint + "_mrna.csv")
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
    corr = np.arctan(s)
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

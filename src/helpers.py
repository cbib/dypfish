#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import math
import pathlib
from typing import List

from colormap import rgb2hex
from numpy import matlib
import numpy as np
from scipy.stats import *
import pandas as pd
import itertools
import constants
from scipy import stats
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint
from mpi_calculator import DensityStats


def checkpoint_decorator(path, dtype):
    def real_checkpoint_decorator(decorated_method):
        def wrapper(instance, *args, **kwargs):
            tgt_path = instance._path + path
            if not instance._repository.is_present(tgt_path):
                # logger.debug("Computing value for {} using {}", tgt_path, decorated_method)
                value = decorated_method(instance, *args, **kwargs)
                instance._repository.save_descriptor(tgt_path, value=value, dtype=dtype)
            else:
                # logger.debug("Reloading value for {} for {}", tgt_path, decorated_method)
                value = np.array(instance._repository.get(tgt_path)).astype(dtype)
            return value

        return wrapper

    return real_checkpoint_decorator


def volume_coeff():
    return ((1 / constants.dataset_config['SIZE_COEFFICIENT']) ** 2) * 0.3


def surface_coeff():
    return ((1 / constants.dataset_config['SIZE_COEFFICIENT']) ** 2)


def rotate_meshgrid(xx, yy, radians=0):
    """
    Rotate a meshgrid counter clockwise by an angle
    """
    R = None
    R = np.array([[np.cos(radians), np.sin(radians)],
                  [-np.sin(radians), np.cos(radians)]])

    return np.einsum('ji, mni -> jmn', R, np.dstack([xx, yy]))


def rotate_point(self, center_point, point, angle):
    """
    counter-clockwise rotate a point around the center_point; angle is in degrees.
    regular cartesian coordinates are used
    """

    angle = math.radians(angle)
    temp_point = [point[0] - center_point[0], point[1] - center_point[1]]
    temp_point = (temp_point[0] * math.cos(angle) - temp_point[1] * math.sin(angle),
                  temp_point[0] * math.sin(angle) + temp_point[1] * math.cos(angle))
    temp_point = np.array([int(round(temp_point[0] + center_point[0])), int(round(temp_point[1] + center_point[1]))])
    return temp_point


def slope_from_points(self, point1, point2):
    if (point2[1] != point1[1]):
        return (point2[0] - point1[0]) / (point2[1] - point1[1])
    else:
        return 0


def find_nearest(array, value):
    """
    Returns the index in the array where the value is the closest to the one in the argument
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def create_dir_if_needed_for_filepath(file_path):
    try:
        pathlib.Path(file_path).parent.mkdir(parents=True)
    except FileExistsError:
        return


def unit_circle(size, r) -> np.ndarray:
    """
    Create a matrix size x size with a circular mask of radius r in the center
    :return: an array of int
    """
    assert (size > r * 2), "Can't create a circular mask of radius greater or equal to the square's side"
    m = np.zeros((size, size))

    # specify circle centre ij
    ci, cj = int(size / 2), int(size / 2)

    # Create index arrays to m
    i, j = np.meshgrid(np.arange(m.shape[0]), np.arange(m.shape[1]))

    # calculate distance of all points to centre
    dist = np.sqrt((i - ci) ** 2 + (j - cj) ** 2)

    # Assign value of 1 to those points where dist<cr:
    m[np.where(dist < r)] = 1

    return m.astype(int)


def compute_statistics_random_h_star(h_sim: np.ndarray, max_cell_radius=None, simulation_number=None) -> (
        List[int], List[int], List[int]):
    """
    Build related statistics derived from Ripley's K function, normalize K
    """
    simulation_number = simulation_number or constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]
    max_cell_radius = max_cell_radius or constants.analysis_config["MAX_CELL_RADIUS"]

    h_sim = np.power((h_sim * 3) / (4 * math.pi), 1. / 3) - matlib.repmat(np.arange(1, max_cell_radius + 1),
                                                                          simulation_number, 1)
    h_sim_sorted = np.sort(h_sim)
    h_sim_sorted = np.sort(h_sim_sorted[:, ::-1], axis=0)
    synth95 = h_sim_sorted[int(np.floor(
        0.95 * simulation_number))]  # TODO : difference with V0 : floor since if the numbers are high we get simulation_sumber here
    synth50 = h_sim_sorted[int(np.floor(0.5 * simulation_number))]
    synth5 = h_sim_sorted[int(np.floor(0.05 * simulation_number))]

    return synth5, synth50, synth95


def compute_h_star(h: np.ndarray, synth5: List[int], synth50: List[int], synth95: List[int],
                   max_cell_radius=None) -> np.ndarray:
    """
    Compute delta between .95 percentile and .5 percentile; between .5 percentile and .05 percentile
    Fill the h_star array accordingly
    """
    max_cell_radius = max_cell_radius or constants.analysis_config["MAX_CELL_RADIUS"]
    delta1 = synth95 - synth50
    delta2 = synth50 - synth5
    idx_equal_median = np.where(h == synth50)[0]
    h_star = np.zeros(max_cell_radius)
    h_star[idx_equal_median] = 0
    idx_greater_median = np.where(h > synth50)[0]
    h_star[idx_greater_median] = (h[idx_greater_median] - synth50[idx_greater_median]) / delta1[idx_greater_median]
    idx_less_median = np.where(h < synth50)[0]
    h_star[idx_less_median] = -(synth50[idx_less_median] - h[idx_less_median]) / delta2[idx_less_median]
    h_star[h_star == - np.inf] = 0
    h_star[h_star == np.inf] = 0
    return h_star


def compute_statistics_random_h_star_2d(h_sim: np.ndarray, max_cell_radius=None, simulation_number=None) -> (
        List[int], List[int], List[int]):
    """
    Build related statistics derived from Ripley's K function, normalize K
    """
    simulation_number = simulation_number or constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]
    max_cell_radius = max_cell_radius or constants.analysis_config["MAX_CELL_RADIUS"]

    h_sim = np.sqrt((h_sim / math.pi)) - matlib.repmat(np.arange(1, max_cell_radius + 1), simulation_number, 1)
    h_sim_sorted = np.sort(h_sim)
    # TODO this line below was in VO
    h_sim_sorted = np.sort(h_sim_sorted[:, :], axis=0)
    # TODO this line below was in V1
    # h_sim_sorted = np.sort(h_sim_sorted[:, ::-1], axis=0)
    synth95 = h_sim_sorted[int(np.floor(
        0.95 * simulation_number))]  # TODO : difference with V0 : floor since if the numbers are high we get simulation_sumber here
    synth50 = h_sim_sorted[int(np.floor(0.5 * simulation_number))]
    synth5 = h_sim_sorted[int(np.floor(0.05 * simulation_number))]

    return synth5, synth50, synth95


def compute_h_star_2d(h: np.ndarray, synth5: List[int], synth50: List[int], synth95: List[int],
                   max_cell_radius=None) -> np.ndarray:
    """
    Compute delta between .95 percentile and .5 percentile; between .5 percentile and .05 percentile
    Fill the h_star array accordingly
    """
    max_cell_radius = max_cell_radius or constants.analysis_config["MAX_CELL_RADIUS"]
    #delta1 = synth95 - synth50
    #delta2 = synth50 - synth5
    idx_equal_median = np.where(h == synth50)[0]
    h_star = np.zeros(max_cell_radius)
    h_star[idx_equal_median] = 0
    idx_greater_median = np.where(h > synth50)[0]
    h_star[idx_greater_median] = (h[idx_greater_median] - synth50[idx_greater_median])
    idx_less_median = np.where(h < synth50)[0]
    h_star[idx_less_median] = -(synth50[idx_less_median] - h[idx_less_median])
    h_star[h_star == - np.inf] = 0
    h_star[h_star == np.inf] = 0
    return h_star

def color_variant(hex_color, brightness_offset=1):
    """ takes a color like #87c95f and produces a lighter or darker variant """
    if len(hex_color) != 7:
        raise Exception("Passed %s into color_variant(), needs to be in #87c95f format." % hex_color)
    rgb_hex = [hex_color[x:x + 2] for x in [1, 3, 5]]
    new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
    new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int]  # make sure new values are between 0 and 255
    # hex() produces "0x88", we want just "88"
    return rgb2hex(new_rgb_int[0], new_rgb_int[1], new_rgb_int[2])


def detect_outliers(data, threshold=3):
    mean = np.nanmean(data)
    std = np.nanstd(data)
    outliers = []
    for x in data:
        z_score = (x - mean) / std

        if np.abs(z_score) > threshold:
            outliers.append(x)
    return outliers


def open_repo():
    dataset_root_fp = pathlib.Path(
        constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
    return repo


####TIS PART

def reindex_quadrant_mask(quad_mask, mtoc_quad, quad_num=4):
    df = pd.DataFrame(quad_mask)
    df = df.applymap(lambda x: x - mtoc_quad + 1 if x >= mtoc_quad else (x + quad_num - mtoc_quad + 1 if x > 0 else 0))
    quad_mask = np.array(df)
    return quad_mask


def permutations(orig_list):
    if not isinstance(orig_list, list):
        orig_list = list(orig_list)
    yield orig_list
    if len(orig_list) == 1:
        return
    for n in sorted(orig_list):
        new_list = orig_list[:]
        pos = new_list.index(n)
        del (new_list[pos])
        new_list.insert(0, n)
        for resto in permutations(new_list[1:]):
            if new_list[:1] + resto != orig_list:
                yield new_list[:1] + resto


def using_indexed_assignment(x):
    "https://stackoverflow.com/a/5284703/190597 (Sven Marnach)"
    result = np.empty(len(x), dtype=int)
    temp = x.argsort()
    result[temp] = np.arange(len(x))
    return result


def permutations_test(interactions, fwdints, size=4):
    fwdints = fwdints.astype(bool)
    vals = interactions.flatten()
    indx = using_indexed_assignment(vals)
    one_matrix = np.ones((size, size)).astype(int)
    indx_matrix = np.matrix(indx.reshape((size, size)))
    indx_matrix = np.add(indx_matrix, one_matrix)
    ranking = indx_matrix.copy()
    rs0 = np.sum(indx_matrix[fwdints[:]])
    rs1 = np.sum(indx_matrix[fwdints[:] == 0])
    perms = [x for x in itertools.permutations(np.arange(size), size)]
    nps = len(perms)
    rs = []
    for p1 in range(nps):
        for p2 in range(nps):
            test = indx_matrix.copy()
            for i in range(size):
                for j in range(size):
                    np.random.shuffle(test[:, i])
            rs.append(np.sum(test[fwdints[:]]))
    count = 0
    for score in rs:
        if score > rs0:
            count += 1
    p = float(count / float(len(rs)))
    stat = rs1
    return p, stat, ranking


def pearsoncorr(vec1, vec2):
    mu1 = np.mean(vec1)
    mu2 = np.mean(vec2)
    vec1b = vec1 - mu1
    vec2b = vec2 - mu2
    val = np.mean(vec1b * vec2b) / (np.std(vec1) * np.std(vec2))
    return val


def get_forward_interactions(mrna_timepoints, protein_timepoints):
    X = len(mrna_timepoints)
    Y = len(protein_timepoints)
    fwd_interactions = np.zeros((X, Y))
    for x in range(X):
        for y in range(Y):
            if protein_timepoints[y] > mrna_timepoints[x]:
                fwd_interactions[x, y] = 1
    return fwd_interactions



##### muscle data helpers functions

def get_quantized_grid(q, Qx, Qy):
    tmp_x = np.matrix(np.arange(Qx))
    tmp_y = np.matrix(np.arange(Qy))
    qxs = matlib.repmat(tmp_x.transpose(), 1, Qx)
    qys = matlib.repmat(tmp_y, Qy, 1)
    qxs = np.kron(qxs, np.ones((q, q)))
    qys = np.kron(qys, np.ones((q, q)))
    return qxs, qys


def reduce_z_line_mask(z_lines, spots):
    cpt_z = 1
    z_lines_idx = []
    for z_line_mask in z_lines:
        spots_reduced = spots[spots[:, 2] == cpt_z]
        if len(spots_reduced) > 25 and len(spots_reduced) < 2000:
            z_lines_idx.append(cpt_z)
        cpt_z += 1
    return z_lines_idx

def compute_minimal_distance(segment_summed):
    for i in range(15):
        if segment_summed[i] != 0:
            return i

def keep_cell_mask_spots(spots, cell_mask):
    new_spots_list = []
    for spot in spots:
        if cell_mask[spot[1], spot[0]] == 1:
            new_spots_list.append(spot)
    return new_spots_list


def build_density_by_stripe(spots_reduced, z_lines, cell_mask, band_n=100):
    z_lines_idx = reduce_z_line_mask(z_lines, spots_reduced)
    spots_reduced = spots_reduced[z_lines_idx[0] <= spots_reduced[:, 2]]
    spots_reduced = spots_reduced[spots_reduced[:, 2] <= z_lines_idx[len(z_lines_idx) - 1]]
    spot_surfacic_density = len(spots_reduced) / float(np.sum(cell_mask == 1))
    cell_width = cell_mask.shape[1] - 240
    quadrat_edge = cell_width / band_n
    grid_1d = np.zeros((int(band_n)))
    for spot in spots_reduced:
        if spot[0] > 120 and spot[0] < cell_mask.shape[1] - 120:
            x = int(np.floor((spot[0] - 120) / quadrat_edge))
            grid_1d[x] += 1
    grid = [val for val in grid_1d]
    grid_mat = np.matrix(grid).reshape((1, len(grid)))
    grid_mat /= quadrat_edge
    grid_mat /= spot_surfacic_density

    return grid_mat


def calculate_temporal_interaction_score(mrna_data, protein_data, timepoint_num_mrna, timepoint_num_protein):
    S1 = get_forward_interactions(timepoint_num_mrna, timepoint_num_protein)
    interactions = np.zeros((len(timepoint_num_mrna), len(timepoint_num_protein)))
    for i in range(len(timepoint_num_mrna)):
        for j in range(len(timepoint_num_protein)):
            interactions[i, j] = stats.pearsonr(list(mrna_data[i]), list(protein_data[j]))[0]
    (p, stat, ranking) = permutations_test(interactions, S1, size=len(timepoint_num_mrna))
    if len(timepoint_num_mrna)==4:
        #TODO if matrix 4 * 4
        tis = (100 - stat) / 64.0
    else:
        # TODO if matrix 2 * 2
        tis = (12 - stat) / 3.0

    return tis, p, ranking

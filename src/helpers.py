#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import itertools
import math
import pathlib
from typing import List
import numpy as np
import pandas as pd
from colormap import rgb2hex
from numpy import matlib
from scipy import stats
from scipy.spatial import cKDTree
from scipy.special import gamma, digamma
from scipy.special import gammainc
import statsmodels.formula.api as smf

import constants
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint


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


def open_repo():
    dataset_root_fp = pathlib.Path(
        constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
    primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
    secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
    repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
    return repo


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


def median_confidence_interval(a: np.array, cutoff=.95):
    ''' cutoff is the significance level as a decimal between 0 and 1'''
    a = np.sort(a)
    factor = stats.norm.ppf((1 + cutoff) / 2)
    factor *= math.sqrt(len(a))  # avoid doing computation twice

    lix = int(0.5 * (len(a) - factor)) + 1
    uix = int(0.5 * (1 + len(a) + factor)) + 1
    print(lix, uix)
    if lix <= 0:
        lix = 0
    if uix <= len(a)-1:
        uix = len(a)-1
    assert (lix <= len(a)-1), "index " + str(lix) + " is out of bound for array of size " + str(len(a))
    assert (uix <= len(a)-1), "index " + str(uix) + " is out of bound for array of size " + str(len(a))

    return a[lix], a[uix]


def sem(a: np.array, factor=3) -> float:
    """
    sem in presence of extreme outliers (very skewed distribution), factor = 0 gives standard behaviour
    """
    if factor > 0:
        limit = factor * np.std(a)
        a = a[(a < np.mean(a) + limit) & (a > np.mean(a) - limit)]
    return stats.sem(a, ddof=0)


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
    median_upper_idx = np.where(h > synth50)[0]
    h_star[median_upper_idx] = (h[median_upper_idx] - synth50[median_upper_idx]) / delta1[median_upper_idx]
    median_lower_idx = np.where(h < synth50)[0]
    h_star[median_lower_idx] = -(synth50[median_lower_idx] - h[median_lower_idx]) / delta2[median_lower_idx]
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
    idx_equal_median = np.where(h == synth50)[0]
    h_star = np.zeros(max_cell_radius)
    h_star[idx_equal_median] = 0
    median_upper_idx = np.where(h > synth50)[0]
    h_star[median_upper_idx] = (h[median_upper_idx] - synth50[median_upper_idx])
    median_lower_idx = np.where(h < synth50)[0]
    h_star[median_lower_idx] = -(synth50[median_lower_idx] - h[median_lower_idx])
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
        if 120 < spot[0] < cell_mask.shape[1] - 120:
            x = int(np.floor((spot[0] - 120) / quadrat_edge))
            grid_1d[x] += 1
    grid = [val for val in grid_1d]
    grid_mat = np.array(grid).reshape((1, len(grid)))
    grid_mat /= quadrat_edge
    grid_mat /= spot_surfacic_density

    return grid_mat


def get_forward_timepoints(mrna_timepoints: list, protein_timepoints: list) -> np.array:
    '''
    Given two lists of timepoints (numerical values), computes pairs
    where protein measurement happened at a later timepoint than mrna measurement
    returns a 2D numpy array where 1 is for protein_timepoint > mrna_timepoint
    '''
    X = len(mrna_timepoints)
    Y = len(protein_timepoints)
    fwd_interactions = np.zeros((X, Y))
    for x in range(X):
        for y in range(Y):
            if protein_timepoints[y] > mrna_timepoints[x]:
                fwd_interactions[x, y] = 1
    return fwd_interactions


def make_categorical(arr):
    '''
    In a density array replaces values by 0, 1, 2
    0 : within the std, 1 > std, 2 < std
    '''
    categorical_arr = np.zeros(arr.shape[0])
    clustered_indices = np.argwhere(arr > np.mean(arr) + np.std(arr)).flatten()
    underclusstered_indices = np.argwhere(arr < np.mean(arr) - np.std(arr)).flatten()
    categorical_arr[clustered_indices] = 1
    categorical_arr[underclusstered_indices] = 2
    return categorical_arr


def get_neighbors(p, exclude_p=False, shape=None):
    ndim = len(p)
    # generate an (m, ndims) array containing all strings over the alphabet {0, 1, 2}:
    offset_idx = np.indices((3,) * ndim).reshape(ndim, -1).T
    # use these to index into np.array([-1, 0, 1]) to get offsets
    offsets = np.r_[-1, 0, 1].take(offset_idx)
    # optional: exclude offsets of 0, 0, ..., 0 (i.e. p itself)
    if exclude_p:
        offsets = offsets[np.any(offsets, 1)]

    neighbours = p + offsets  # apply offsets to p

    # optional: exclude out-of-bounds indices
    if shape is not None:
        valid = np.all(neighbours < np.array(shape), axis=1) & (neighbours[:, 0] >= 0)
        neighbours = neighbours[valid]

    return neighbours


def neighboring_protein_values(mrna, protein, stripes, quadrants):
    '''
    For each pair of indices i,j in the mrna density vector, collect protein value
    from the protein densities vector in the direct neighborhood of i,j and
    closest in value to the mrna density at this index
    '''
    mrna_stripes = mrna.reshape((stripes, quadrants))
    protein_stripes = protein.reshape((stripes, quadrants))
    padded_stripes = np.pad(protein_stripes, (1, 1), mode='wrap')
    padded_stripes[0, :] = -1;
    padded_stripes[padded_stripes.shape[0] - 1, :] = -1
    all_neighbours = {}
    for i, j in itertools.product(range(1, stripes + 1), range(1, quadrants + 1)):
        neighbors_idx = get_neighbors(np.r_[i, j], shape=padded_stripes.shape)
        neighbors = padded_stripes[tuple(neighbors_idx.T)]
        all_neighbours[(i - 1, j - 1)] = neighbors[neighbors != -1]

    protein_nearest_neighbours = []
    for i, j in itertools.product(range(stripes), range(quadrants)):
        protein_vals = all_neighbours[(i, j)]
        idx = (np.abs(protein_vals - mrna_stripes[i, j])).argmin()
        protein_nearest_neighbours.append(protein_vals[idx])

    return np.array(protein_nearest_neighbours).astype(int)


def neighboring_protein_values_periphery(mrna, protein):
    protein_nearest_neighbours = []
    for i, val in enumerate(mrna):
        indices = range(i - 1, i + 2)
        neighbourhood = protein.take(indices, mode='wrap')
        idx = (np.abs(neighbourhood - val).argmin())
        protein_nearest_neighbours.append(neighbourhood[idx])

    return np.array(protein_nearest_neighbours).astype(int)


def quantile_regression(categorical_mrna, categorical_protein):
    data = pd.DataFrame(columns=['mrna', 'protein'])
    data['mrna'] = categorical_mrna
    data['protein'] = categorical_protein
    mod = smf.quantreg('mrna ~ protein', data)
    res = mod.fit(q=.5)
    return res.prsquared


def calculate_colocalization_score(mrna_data, protein_data, timepoint_num_mrna,
                                   timepoint_num_protein, peripheral_flag, stripes, quadrants):
    num_mrna_tp, num_protein_tp = len(timepoint_num_mrna), len(timepoint_num_protein)
    correlations = np.zeros((num_mrna_tp, num_protein_tp))
    for i, j in itertools.product(range(num_mrna_tp), range(num_protein_tp)):
        categorical_mrna = make_categorical(mrna_data[i])
        categorical_protein = make_categorical(protein_data[j])
        if not peripheral_flag:
            neighbouring_protein = neighboring_protein_values(categorical_mrna, categorical_protein,
                                                              stripes, quadrants)
            correlations[i, j] = stats.pearsonr(categorical_mrna, neighbouring_protein)[0]
        else:
            neighbouring_protein = neighboring_protein_values_periphery(categorical_mrna[0:quadrants],
                                                                        categorical_protein[0:quadrants])
            if (not np.all(categorical_mrna == neighbouring_protein)):
                correlations[i, j] = stats.pearsonr(categorical_mrna[0:quadrants],
                                                    neighbouring_protein[0:quadrants])[0]
            else:  # neighbour-based approach failed, this is a hack
                correlations[i, j] = (categorical_mrna == categorical_protein).sum() / len(categorical_mrna)

    fwd_indices = get_forward_timepoints(timepoint_num_mrna, timepoint_num_protein)
    fwd_correlations = correlations[fwd_indices == 1].flatten()
    other_correlations = correlations[fwd_indices == 0].flatten()
    ranks = stats.rankdata(correlations, method='ordinal').reshape(correlations.shape[0], correlations.shape[1])
    stat, pval = stats.mannwhitneyu(fwd_correlations, other_correlations, use_continuity=False)
    cs = a12(list(fwd_correlations), list(other_correlations))  # effect size

    return cs, pval, ranks


# Stability analysis part

def mean_absolute_deviation(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


def median_absolute_deviation(data, axis=None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)


# color helper

def clamp(val, minimum=0, maximum=255):
    if val < minimum:
        return minimum
    if val > maximum:
        return maximum
    return int(val)


def colorscale(hexstr, scalefactor):
    """
    Scales a hex string by ``scalefactor``. Returns scaled hex string.

    To darken the color, use a float value between 0 and 1.
    To brighten the color, use a float value greater than 1.

    >>> colorscale("#DF3C3C", .5)
    #6F1E1E
    >>> colorscale("#52D24F", 1.6)
    #83FF7E
    >>> colorscale("#4F75D2", 1)
    #4F75D2
    """

    hexstr = hexstr.strip('#')

    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr

    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)

    r = clamp(r * scalefactor)
    g = clamp(g * scalefactor)
    b = clamp(b * scalefactor)

    return "#%02x%02x%02x" % (r, g, b)


def a12(lst1, lst2, rev=True):
    '''
      Non-parametric hypothesis testing using Vargha and Delaney's A12 statistic:
      how often is x in lst1 greater than y in lst2?
    '''
    more = same = 0.0
    for x in lst1:
        for y in lst2:
            if x == y:
                same += 1
            elif rev and x > y:
                more += 1
            elif not rev and x < y:
                more += 1
    return (more + 0.5 * same) / (len(lst1) * len(lst2))


def random_points_in_sphere(center, radius, num_points) -> np.ndarray:
    '''
    Generate num_points random points within a shpere with center in center
    center is a numpy array
    '''
    r = radius
    ndim = center.size
    x = np.random.normal(size=(num_points, ndim))
    ssq = np.sum(x ** 2, axis=1)
    fr = r * gammainc(ndim / 2, ssq / 2) ** (1 / ndim) / np.sqrt(ssq)
    frtiled = np.tile(fr.reshape(num_points, 1), (1, ndim))
    p = center + np.multiply(x, frtiled)
    return p


def roll_densities_mtoc_array(densities, slices=3):
    '''
    Given a 2d numpy array with densities for all slices, where
    first column contains densities and second column encods the mtoc containing quadrants,
    roll all the slices so that the mtoc containing quadrant is the first for all slices
    '''
    slices = np.split(densities, slices)
    for idx, slice in enumerate(slices):
        mtoc_quadrant = np.argwhere(slice[:, 1] == 1)[0][0]
        slice = np.roll(slice, -mtoc_quadrant, 0)
        slices[idx] = slice
    return np.vstack(slices)


def compute_entropy(x, k=1, norm='max', min_dist=0.):
    """
    Estimates the entropy H of a random variable x (in nats) based on
    the kth-nearest neighbour distances between point samples.
    Implementation credits: Paul Brodersen
    @reference:
    Kozachenko, L., & Leonenko, N. (1987). Sample estimate of the entropy of a random vector.
    Problemy Peredachi Informatsii, 23(2), 9–16.
    Arguments:
    ----------
    x: (n, d) ndarray
        n samples from a d-dimensional multivariate distribution
    k: int (default 1)
        kth nearest neighbour to use in density estimate;
        imposes smoothness on the underlying probability distribution
    norm: 'euclidean' or 'max'
        p-norm used when computing k-nearest neighbour distances
    min_dist: float (default 0.)
        minimum distance between data points;
        smaller distances will be capped using this value
    Returns:
    --------
    h: float
        entropy H(X)
    """

    n, d = x.shape

    if norm == 'max':  # max norm:
        p = np.inf
        log_c_d = 0  # volume of the d-dimensional unit ball
    elif norm == 'euclidean':  # euclidean norm
        p = 2
        log_c_d = (d / 2.) * np.log(np.pi) - np.log(gamma(d / 2. + 1))
    else:
        raise NotImplementedError("Variable 'norm' either 'max' or 'euclidean'")

    kdtree = cKDTree(x)

    # query all points -- k+1 as query point also in initial set
    distances, _ = kdtree.query(x, k + 1, eps=0, p=p)
    distances = distances[:, -1]

    # enforce non-zero distances
    distances[distances < min_dist] = min_dist

    sum_log_dist = np.sum(np.log(2 * distances))  # where did the 2 come from? radius -> diameter
    h = -digamma(k) + digamma(n) + log_c_d + (d / float(n)) * sum_log_dist

    return h

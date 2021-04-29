#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import helpers
import numpy as np
from repository import Repository
from image import Image, ImageWithSpots, ImageWithIntensities, ImageWithMTOC, ImageMultiNucleus, \
    ImageMultiNucleusWithSpots
import numexpr
from scipy import signal
from loguru import logger
import constants
import math
import tqdm

from constants import HEIGHT_MAP_PATH_SUFFIX
from constants import CELL_MASK_SLICES_PATH_SUFFIX
from constants import ZERO_LEVEL_PATH_SUFFIX
from constants import NUCLEUS_VOLUME_PATH_SUFFIX
from constants import CELL_VOLUME_PATH_SUFFIX
from constants import SPOTS_PERIPHERAL_DISTANCE_3D_PATH_SUFFIX
from constants import CLUSTERING_INDICES_PATH_SUFFIX
from constants import NUCLEUS_CENTROID_PATH_SUFFIX

class Image3d(Image):
    """ Represents an 3D image, has to have a height map descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        # TODO : check if we need zero level to define a 3D cell.
        #  comment ZERO_LEVEL_PATH_SUFFIX is present() for Cultured data to run volume corrected analysis.
        #  they do not have ZERO LEVEL descriptors and so was rejected as a 3D image
        return repo.is_present(path + HEIGHT_MAP_PATH_SUFFIX)

    def __init__(self, repository: Repository, image_path: str):
        super(Image3d, self).__init__(repository, image_path)
        if not self._repository.is_present(image_path + HEIGHT_MAP_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    def get_height_map(
            self) -> np.ndarray:  # TODO : not restricted to the cell_mask in V0 (same as for the intensities)
        descriptor = self._path + HEIGHT_MAP_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No height map for image %s" % self._path)
        return np.array(self._repository.get(descriptor))

    def get_cytoplasm_height_map(self):
        height_map = self.get_height_map()
        height_map[self.get_cytoplasm_mask() == 0] = 0
        return height_map

    def adjust_height_map(self, cytoplasm=False):
        '''
        adjust the periphery of the cell to be at 0.5 height
        '''
        if cytoplasm == True:
            height_map = self.get_cytoplasm_height_map().astype(float)
            cell_mask = self.get_cytoplasm_mask()
        else:
            height_map = self.get_height_map().astype(float)
            cell_mask = self.get_cell_mask()

        height_map[(cell_mask == 1) & (height_map == 0)] = 0.5
        return height_map

    def get_zero_level(self):
        descriptor = self._path + ZERO_LEVEL_PATH_SUFFIX
        if not self._repository.is_present(descriptor):
            raise LookupError("No zero level for image %s" % self._path)
        return np.array(self._repository.get(descriptor))

    def compute_cell_mask_slices(self, cell_mask=None, height_map=None,
                                 zero_level=None, image_width=None, image_height=None) -> np.ndarray:
        """
        Reconstructs the z-slices given a height_map;
        out of focus slices (defined by zero_level) are not reconstructed
        Was : compute_cell_mask_3d
        :return:
        """
        if cell_mask is None: cell_mask = self.get_cell_mask()
        if height_map is None: height_map = self.get_height_map()
        zero_level = zero_level or self.get_zero_level()
        image_width = image_width or constants.dataset_config["IMAGE_WIDTH"]
        image_height = image_height or constants.dataset_config["IMAGE_HEIGHT"]
        height_map = height_map.astype(float)  # TODO : make this code work without converting to float and using np.nan
        height_map[cell_mask == 0] = np.nan
        reversed_height_map = zero_level - height_map + 1

        if np.nanmin(reversed_height_map) < 0:
            logger.debug("reversed_height_map has negative values for {}", self._path)

        # Create binary cell masks per slice
        cell_masks = np.zeros((image_width, image_height, zero_level))

        # build slice mask
        for slice in range(0, zero_level):
            slice_mask = np.array(reversed_height_map, copy=True)
            with np.errstate(
                    invalid='ignore'):  # TODO : this is to avoid RuntimeWarning due to the comparison with np.nan
                slice_mask[slice_mask > slice] = np.nan
                slice_mask[slice_mask <= float(slice)] = 1
            slice_mask = np.nan_to_num(slice_mask)
            cell_masks[:, :, slice] = slice_mask
        return cell_masks

    @helpers.checkpoint_decorator(CELL_MASK_SLICES_PATH_SUFFIX, dtype=np.int)
    def get_cell_mask_slices(self) -> np.ndarray:
        return self.compute_cell_mask_slices()

    @helpers.checkpoint_decorator(CELL_VOLUME_PATH_SUFFIX, float)
    def get_cell_volume(self) -> float:
        return self.compute_cell_volume()

    def get_peripheral_cell_volume(self) -> float:
        return self.compute_peripheral_cell_volume()

    def compute_cell_volume(self):
        """
        compute cell surface in pixels using the cell mask
        """
        volume_offset = constants.dataset_config['VOLUME_OFFSET']
        cell_mask = self.get_cell_mask()
        height_map = self.get_height_map()
        height_map = height_map + volume_offset
        cell_volume = height_map[np.where(cell_mask[:] == 1)].sum()
        return cell_volume * helpers.volume_coeff()

    def compute_peripheral_cell_volume(self):
        """
        compute cell surface in pixels using the cell mask
        """
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        volume_offset = constants.dataset_config['VOLUME_OFFSET']
        cell_mask = self.get_cell_mask()
        height_map = self.get_height_map()
        height_map = height_map + volume_offset
        height_map_periph = np.multiply(height_map, peripheral_binary_mask)
        peripheral_cell_volume = height_map_periph[np.where(cell_mask[:] == 1)].sum()
        return peripheral_cell_volume * helpers.volume_coeff()

    @helpers.checkpoint_decorator(NUCLEUS_VOLUME_PATH_SUFFIX, float)
    def get_nucleus_volume(self) -> float:
        return self.compute_nucleus_volume()

    def compute_nucleus_volume(self):
        """compute nucleus volume in pixel using nucleus mask"""
        volume_offset = constants.dataset_config['VOLUME_OFFSET']
        nucleus_mask = self.get_nucleus_mask()
        height_map = self.get_height_map()
        height_map = height_map + volume_offset
        nucleus_volume = height_map[np.where(nucleus_mask[:] == 1)].sum()
        return nucleus_volume * helpers.volume_coeff()

    def compute_cytoplasmic_volume(self):
        cell_volume = self.get_cell_volume()
        nucleus_volume = self.get_nucleus_volume()
        cytoplasmic_volume = (cell_volume - nucleus_volume)
        if (cytoplasmic_volume <= 0):
            raise (RuntimeError, "cytoplasm and nucleus have inconsistent volumes for image %s" % self._path)
        return cytoplasmic_volume

    def compute_peripheral_volume(self):
        peripheral_cell_volume = self.get_peripheral_cell_volume()
        if (peripheral_cell_volume <= 0):
            raise (RuntimeError, "peripheral area have inconsistent volumes for image %s" % self._path)
        return peripheral_cell_volume


class Image3dWithSpots(Image3d, ImageWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3d.is_a(repo, path) and ImageWithSpots.is_a(repo, path)

    def compute_spots_peripheral_distance_3d(self) -> np.ndarray:
        """
        Perform the computation in pseudo 3D using the height map
        Return an array of distances
        """
        height_map = self.get_height_map()
        zero_level = self.get_zero_level()
        peripheral_distance_map = self.get_cell_mask_distance_map()
        spots_peripheral_distance = []
        spots = self.get_cytoplasmic_spots()  # was get_spots(), changed to be coherent with the 2D version
        logger.info("Computing 3D peripheral distance for {} spots in image {}",
                    len(spots), self._path)
        peripheral_areas = self.compute_peripheral_areas()
        problematic_spot_num = 0
        for slice_num in range(zero_level, -1, -1):
            height_map_copy = np.array(height_map, copy=True)
            height_map_copy[height_map >= zero_level + 1 - slice_num] = 1
            # TODO : maybe better stay in pixels ?
            slice_area = height_map_copy[height_map_copy == 1].sum() * helpers.surface_coeff()
            slice_index = helpers.find_nearest(peripheral_areas, slice_area)
            spots_in_slice = spots[np.around(spots[:, 2]) == slice_num]
            for j in range(len(spots_in_slice)):
                old_periph_distance = peripheral_distance_map[spots_in_slice[j][1], spots_in_slice[j][0]]
                new_periph_distance = (old_periph_distance - slice_index) * (
                        101 / (101 - slice_index))  # 101 to avoid divide by 0 error

                if new_periph_distance < 0:
                    # this means that the spot falls out of the slice, it is a discrepancy
                    # between the 2D calculation and the slice area estimation
                    spots_peripheral_distance.append(int(1))
                    problematic_spot_num = problematic_spot_num + 1
                else:
                    spots_peripheral_distance.append(int(np.around(new_periph_distance)))

        if problematic_spot_num > 0:
            logger.debug("Peripheral spot distance could not be be determined in 3d for {} spots {}",
                         problematic_spot_num, self._path)

        return np.asarray(spots_peripheral_distance, dtype=np.uint8)

    def compute_average_distance_from_nucleus(self, dsAll, cell_mask) -> float:
        height_map = self.get_height_map()
        height_map += 1  # TODO : why is this ???
        height_map = np.multiply(height_map, cell_mask)
        height_map_dist = np.multiply(height_map, dsAll)

        # S : Average distance of a cytoplasmic voxel from the nucleus centroid
        # TODO : this is false since height_map_dist is computed over all the cell and not the cytoplasm
        return height_map_dist.sum() / height_map.sum()

    @helpers.checkpoint_decorator(SPOTS_PERIPHERAL_DISTANCE_3D_PATH_SUFFIX, dtype=np.float)
    def get_spots_peripheral_distance(self):
        return self.compute_spots_peripheral_distance_3d()

    def ripley_k_point_process(self, nuw: float, my_lambda: float, spots=None, r_max: int = None) -> np.ndarray:
        if spots is None: spots = self.get_spots()
        n_spots = len(spots)
        r_max = r_max or constants.analysis_config["MAX_CELL_RADIUS"]
        pixels_in_slice = numexpr.evaluate(constants.dataset_config["PIXELS_IN_SLICE"]).item()

        K = np.zeros(r_max)
        for i in range(n_spots):
            # TODO : why only z_squared is computed in real size ?
            mask = np.zeros((n_spots, 3));
            mask[i, :] = 1
            other_spots = np.ma.masked_where(mask == 1, np.ma.array(spots, mask=False)).compressed().reshape(
                n_spots - 1, 3)
            x_squared = np.square(spots[i, 0] - other_spots[:, 0])
            y_squared = np.square(spots[i, 1] - other_spots[:, 1])
            z_squared = np.square(pixels_in_slice * (spots[i, 2] - other_spots[:, 2]))
            ds = np.sqrt(x_squared + y_squared + z_squared)
            if n_spots - 1 < r_max:
                for m in range(n_spots - 1):
                    K[math.ceil(ds[m]):r_max] = K[math.ceil(ds[m]):r_max] + 1
            else:
                for m in range(r_max):
                    K[m] = K[m] + ds[ds <= m].sum()
        K = K * (1 / (my_lambda ** 2 * nuw))
        return K

    def compute_random_spots_in_slices(self):
        n_spots = len(self.get_spots())
        slices = self.get_cell_mask_slices()
        x, y, z = np.where(slices == 1)
        idx = np.random.randint(0, len(x), n_spots)  # we chose random indices
        return np.vstack((x[idx], y[idx], z[idx])).T

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripkey-K computation for disks of radius r < MAX_CELL_RADIUS
        :return: clustering indices for all r
        Was : clustering_index_point_process
        """
        logger.info("Running {} simulations of Ripley-K for {}",
                    constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], self._path)
        spots = self.get_spots()
        n_spots = len(spots)
        pixels_in_slice = numexpr.evaluate(constants.dataset_config["PIXELS_IN_SLICE"]).item()

        cell_mask_slices = self.get_cell_mask_slices()
        nuw = (np.sum(cell_mask_slices[:, :, :] == 1)) * pixels_in_slice  # whole volume of the cell
        my_lambda = float(n_spots) / float(nuw)  # spot's volumic density

        k = self.ripley_k_point_process(nuw=nuw, my_lambda=my_lambda)  # TODO : first call for _all_ spots while the subsequent only for those in the height_map
        k_sim = np.zeros((constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], constants.analysis_config["MAX_CELL_RADIUS"]))

        # simulate RIPLEY_K_SIMULATION_NUMBER lists of random spots and run ripley_k
        for t in tqdm.tqdm(range(constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]), desc="Simulations"):
            random_spots = self.compute_random_spots_in_slices()
            tmp_k = self.ripley_k_point_process(spots=random_spots, nuw=nuw, my_lambda=my_lambda).flatten()
            k_sim[t] = tmp_k

        h = np.subtract(np.power(((k * 3) / (4 * math.pi)), 1. / 3),
                        np.arange(1, constants.analysis_config["MAX_CELL_RADIUS"] + 1))
        synth5, synth50, synth95 = helpers.compute_statistics_random_h_star(h_sim=k_sim)
        return helpers.compute_h_star(h, synth5, synth50, synth95)

    @helpers.checkpoint_decorator(CLUSTERING_INDICES_PATH_SUFFIX, dtype=np.float)
    def get_clustering_indices(self):
        return self.compute_clustering_indices()

    def compute_degree_of_clustering(self) -> int:
        h_star = self.get_clustering_indices()
        d_of_c = np.array(h_star[h_star > 1] - 1).sum()
        if int(d_of_c) == 0:
            return 0.0001 # TODO this is a hack so that a downstream log does not fail

        return d_of_c

    def compute_peripheral_density(self):
        # compute mRNA density in the peripheral area
        peripheral_mrna_count = self.get_peripheral_total_spots()
        if peripheral_mrna_count == 0:
            raise RuntimeError("Image contains no spots in periphery %s" % self._path)
        peripheral_volume = self.compute_peripheral_volume()
        return peripheral_mrna_count / peripheral_volume

    def compute_cytoplasmic_density(self):
        # compute mRNA density in the cytoplasm
        cytoplasmic_mrna_count = self.get_cytoplasmic_total_spots()
        if cytoplasmic_mrna_count == 0:
            raise RuntimeError("Image contains no spots %s" % self._path)
        cytoplasmic_volume = self.compute_cytoplasmic_volume()
        return cytoplasmic_mrna_count / cytoplasmic_volume

    # Compare cytoplasmic spread cell with 3D cytoplasmic mrna spread
    # to evaluate degree of spread
    def compute_spots_cytoplasmic_spread(self):
        cytoplasm_mask = self.get_cytoplasm_mask()
        height_map = self.adjust_height_map(cytoplasm=True)
        nucleus_centroid = self.get_nucleus_centroid()
        spots = self.get_cytoplasmic_spots()

        # Compute all possible distance in a matrix [512x512]
        ds1 = np.matlib.repmat(range(0, constants.dataset_config['IMAGE_WIDTH']),
                               constants.dataset_config['IMAGE_WIDTH'], 1) - nucleus_centroid[0]
        ds2 = np.matlib.repmat(np.asmatrix(
            np.arange(0, constants.dataset_config['IMAGE_HEIGHT']).reshape(constants.dataset_config['IMAGE_HEIGHT'],
                                                                           1)), 1,
            constants.dataset_config['IMAGE_HEIGHT']) - nucleus_centroid[1]
        dsAll = np.power(ds1, 2) + np.power(ds2, 2)
        dsAll = np.sqrt(dsAll)

        # Computing spots distance from nucleus centroid
        points_dist_list = []
        counter = 0
        for i in range(len(spots)):
            if cytoplasm_mask[spots[i, 1], spots[i, 0]] == 1:
                dist = 0.0
                for j in range(2):
                    if j == 0:
                        dist += (spots[i, j] - nucleus_centroid[0]) ** 2
                    elif j == 1:
                        dist += (spots[i, j] - nucleus_centroid[1]) ** 2
                points_dist_list.append(math.sqrt(dist))
                counter += 1
        points_dists = np.array(points_dist_list)
        points_dists = points_dists.reshape((counter, 1))
        height_map_dist = np.multiply(height_map, dsAll)
        # S : Average distance of a cytoplasmic voxel from the nucleus centroid
        S = height_map_dist.sum() / height_map.sum()

        # val is average 2D distance from the nucleus centroid of cytoplasmic mRNAs
        # normalized by the cytoplasmic cell spread (taking a value 1 when mRNAs are evenly
        # distributed across the cytoplasm).
        normalized_average_2d_distance = np.mean(points_dists) / S
        return normalized_average_2d_distance


class Image3dWithIntensities(Image3d, ImageWithIntensities):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3d.is_a(repo, path) and ImageWithIntensities.is_a(repo, path)

    def compute_average_distance_proportional_intensity(self, dsAll, cell_mask) -> float:
        height_map = self.get_height_map()
        height_map += 1
        height_map = np.multiply(height_map, cell_mask)
        dsCellular = np.multiply(height_map, dsAll)
        return dsCellular.sum() / height_map.sum()  # TODO : what for ?

    def compute_cytoplasmic_density(self):
        # compute signal density of the cytoplasm
        cytoplasmic_intensity_count = self.get_cytoplasmic_total_intensity()
        cytoplasmic_volume = self.compute_cytoplasmic_volume()
        return cytoplasmic_intensity_count / cytoplasmic_volume

    def compute_peripheral_density(self):
        # compute mRNA density in the peripheral area
        peripheral_intensity_count = self.get_peripheral_total_intensity()
        peripheral_volume = self.compute_peripheral_volume()
        return peripheral_intensity_count / peripheral_volume

    def compute_cell_density(self):
        # compute density of the cell
        intensity_count = self.get_total_intensity()
        volume = self.get_cell_volume()
        return intensity_count / volume

    def compute_intensities_cytoplasmic_spread(self):
        height_map = self.adjust_height_map(cytoplasm=True)  #get_height_map()
        nucleus_centroid = self.get_nucleus_centroid()
        IF = self.get_cytoplasmic_intensities()
        ds1 = np.matlib.repmat(range(0, constants.dataset_config['IMAGE_WIDTH']),
                               constants.dataset_config['IMAGE_WIDTH'], 1) - nucleus_centroid[0]
        ds2 = np.matlib.repmat(np.asmatrix(
            np.arange(0, constants.dataset_config['IMAGE_HEIGHT']).reshape(constants.dataset_config['IMAGE_HEIGHT'],
                                                                           1)), 1,
            constants.dataset_config['IMAGE_HEIGHT']) - nucleus_centroid[1]

        dsAll = np.power(ds1, 2) + np.power(ds2, 2)
        dsAll = np.sqrt(dsAll)
        height_map_dist = np.multiply(height_map, dsAll)
        S = height_map_dist.sum() / height_map.sum()
        dist_IF = np.multiply(IF, dsAll)
        val = dist_IF.sum() / (IF.sum() * S)
        return val

    def compute_clustering_indices(self) -> np.ndarray:
        """
        Point process Ripkey-K computation for disks of radius r < MAX_CELL_RADIUS
        :return: clustering indices for all r
        Was : clustering_index_point_process
        """
        logger.info("Running {} simulations of Ripley-K for {}",
                    constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], self._path)

        pixels_in_slice = numexpr.evaluate(constants.dataset_config["PIXELS_IN_SLICE"]).item()
        IF = self.get_intensities()
        cell_mask = self.get_cell_mask()
        IF = IF.astype(float) * cell_mask
        # TODO in VO we do not multiply by pixels_in_slice ???
        nuw = (np.sum(cell_mask[:, :] == 1))  # * pixels_in_slice  # whole surface of the cell
        my_lambda = float(np.sum(IF)) / float(nuw)  # volumic density
        k = self.ripley_k_random_measure_2D(IF, my_lambda, nuw)
        k_sim = np.zeros(
            (constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"], constants.analysis_config["MAX_CELL_RADIUS"]))
        # simulate RIPLEY_K_SIMULATION_NUMBER list of random intensities and run ripley_k
        indsAll = np.where(cell_mask[:, :] == 1)
        for t in tqdm.tqdm(range(constants.analysis_config["RIPLEY_K_SIMULATION_NUMBER"]), desc="Simulations"):
            inds_permuted = np.random.permutation(range(len(indsAll[0])))
            I_samp = np.zeros(IF.shape)
            for u in range(len(inds_permuted)):
                I_samp[indsAll[0][inds_permuted[u]], indsAll[1][inds_permuted[u]]] = IF[indsAll[0][u], indsAll[1][u]]
            k_sim[t, :] = self.ripley_k_random_measure_2D(I_samp, my_lambda, nuw).flatten()

        h = np.subtract(np.sqrt(k / math.pi), np.arange(1, constants.analysis_config["MAX_CELL_RADIUS"] + 1).reshape(
            (constants.analysis_config["MAX_CELL_RADIUS"], 1))).flatten()
        synth5, synth50, synth95 = helpers.compute_statistics_random_h_star_2d(k_sim)
        return helpers.compute_h_star_2d(h, synth5, synth50, synth95)

    def ripley_k_random_measure_2D(self, IF, my_lambda, nuw):
        IF_rev = IF[::-1, ::-1]
        P = signal.convolve(IF, IF_rev)
        dMap = np.zeros((P.shape[0], P.shape[1]))
        p, q = np.meshgrid(range(P.shape[0]), range(P.shape[0]))
        dMap = np.sqrt((p - IF.shape[0]) ** 2 + (q - IF.shape[1]) ** 2)
        # sum convolution using dMap
        K = np.zeros((constants.analysis_config["MAX_CELL_RADIUS"], 1))
        for dist in range(constants.analysis_config["MAX_CELL_RADIUS"]):
            K[dist] = P[dMap[:, :] <= dist].sum()
        K = K * (1 / (my_lambda * nuw)) - (1 / my_lambda)

        return K


class Image3dWithMTOC(Image3d, ImageWithMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithMTOC.is_a(repo, path) and Image3d.is_a(repo, path)

    def compute_peripheral_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4):
        """
        compute volumic density per quadrant;
        return values of density paired with the MTOC presence flag (0/1)
        Note that this is not relative density
        """
        peripheral_fraction_threshold = constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"]
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        quadrant_mask = quadrant_mask * peripheral_binary_mask
        return self.compute_density_per_quadrant(mtoc_quad, quadrant_mask, quadrants_num)

    def compute_density_per_quadrant_and_slices(self, mtoc_quad, quadrant_mask, stripes, quadrants_num=4,
                                                peripheral_flag=False):
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        slices_per_stripe = np.floor(100.0 / stripes)  # number of isolines per stripe
        if peripheral_flag:
            slices_per_stripe = np.floor(constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"] / stripes)
        arr = np.empty((0, 2), float)
        for stripe_num in range(1, stripes + 1):
            stripe_mask = (cell_mask_dist_map > (stripe_num - 1) * slices_per_stripe) & \
                          (cell_mask_dist_map <= stripe_num * slices_per_stripe + 1).astype(int)
            stripe_quadrant_mask = quadrant_mask * stripe_mask
            res = self.compute_density_per_quadrant(mtoc_quad, stripe_quadrant_mask, quadrants_num)
            arr = np.append(arr, res[res[:, 1].argsort()[::-1]], axis=0)  # MTOC quadrant slice always first
        assert arr.shape[0] == quadrants_num * stripes, "Incorrect shape"
        return arr

    def compute_peripheral_density_per_quadrant_and_slices(self, mtoc_quad_num, quadrant_mask, stripes, quadrants_num=4):
        peripheral_fraction_threshold = constants.analysis_config["PERIPHERAL_FRACTION_THRESHOLD"]
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        quadrant_mask = quadrant_mask * peripheral_binary_mask
        return self.compute_density_per_quadrant_and_slices(mtoc_quad_num, quadrant_mask, stripes, quadrants_num,
                                                            peripheral_flag=True)


class Image3dWithSpotsAndMTOC(Image3dWithMTOC, Image3dWithSpots):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithMTOC.is_a(repo, path) and Image3dWithSpots.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4) -> np.ndarray:
        """
        compute volumic density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        volume_coeff = helpers.volume_coeff()
        height_map = self.adjust_height_map(cytoplasm=True)
        spots = self.get_cytoplasmic_spots()
        density_per_quadrant = np.zeros((quadrants_num, 2))
        for spot in spots:
            spot_quad = quadrant_mask[spot[1], spot[0]]
            if spot_quad == 0: continue
            density_per_quadrant[spot_quad - 1, 0] += 1

        if (density_per_quadrant[:,0].sum() == 0):
            logger.warning("No spots in image within quadrants {}", self._path)
            return density_per_quadrant

        # mark the mtoc quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            quadrant_volume = np.sum(height_map[quadrant_mask == quad_num + 1]) * volume_coeff
            density_per_quadrant[quad_num, 0] = density_per_quadrant[quad_num, 0] / quadrant_volume

        if density_per_quadrant[:, 1].sum() != 1.0:
            raise (RuntimeError, "error in the MTOC quadrant detection for image %s" % self._path)

        return density_per_quadrant


class Image3dWithIntensitiesAndMTOC(Image3dWithMTOC, Image3dWithIntensities):

    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithMTOC.is_a(repo, path) and Image3dWithIntensities.is_a(repo, path)

    def compute_density_per_quadrant(self, mtoc_quad, quadrant_mask, quadrants_num=4) -> np.ndarray:
        """
        compute volumic density per quadrant;
        return an array of values of density paired with the MTOC presence flag (0/1)
        """
        IF = self.get_cytoplasmic_intensities()
        height_map = self.adjust_height_map(cytoplasm=True)
        density_per_quadrant = np.zeros((quadrants_num, 2))
        # mark the MTOC quadrant
        density_per_quadrant[mtoc_quad - 1, 1] = 1
        for quad_num in range(quadrants_num):
            height = np.sum(height_map[quadrant_mask == quad_num + 1])
            quadrant_volume = height * helpers.volume_coeff()
            density_per_quadrant[quad_num, 0] = np.sum(IF[quadrant_mask == (quad_num + 1)]) / quadrant_volume

        if density_per_quadrant[:, 1].sum() != 1.0:
            raise (RuntimeError, "error in the MTOC quadrant detection for image %s" % self._path)

        return density_per_quadrant


class Image3dWithSpotsAndIntensitiesAndMTOC(Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithSpotsAndMTOC.is_a(repo, path) and Image3dWithIntensitiesAndMTOC.is_a(repo, path)


class Image3dMultiNucleus(Image3d):
    """ Represents a generic image with one cell mask and multiple nucleus. It has at least a cell_mask, a nucleus_mask and one or more nucleus_centroid """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return repo.is_multiple(path, path + NUCLEUS_CENTROID_PATH_SUFFIX) and Image3d.is_a(repo, path)

    def __init__(self, repository: Repository, image_path: str):
        super(Image3dMultiNucleus, self).__init__(repository, image_path)
        if not self._repository.is_multiple(image_path, image_path + NUCLEUS_CENTROID_PATH_SUFFIX):
            raise AttributeError("Incorrect format for image %s" % image_path)

    @helpers.checkpoint_decorator(CELL_VOLUME_PATH_SUFFIX, float)
    def get_cell_volume(self) -> float:
        return self.compute_cell_volume()

    def get_peripheral_cell_volume(self) -> float:
        return self.compute_peripheral_cell_volume()

    def compute_cell_volume(self):
        """
        compute cell surface in pixels using the cell mask
        """
        volume_offset = constants.dataset_config['VOLUME_OFFSET']
        cell_mask = self.get_cell_mask()
        height_map = self.get_height_map()
        height_map = height_map + volume_offset
        cell_volume = height_map[np.where(cell_mask[:] == 1)].sum()
        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return cell_volume / len(self.get_multiple_nucleus_centroid()) * helpers.volume_coeff()
        else:
            return cell_volume * helpers.volume_coeff()

    def compute_peripheral_cell_volume(self):
        """
        compute cell surface in pixels using the cell mask
        """
        peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
        volume_offset = constants.dataset_config['VOLUME_OFFSET']
        cell_mask_dist_map = self.get_cell_mask_distance_map()
        peripheral_binary_mask = (cell_mask_dist_map > 0) & \
                                 (cell_mask_dist_map <= peripheral_fraction_threshold).astype(int)
        cell_mask = self.get_cell_mask()
        height_map = self.get_height_map()
        height_map = height_map + volume_offset
        height_map_periph = np.multiply(height_map, peripheral_binary_mask)
        peripheral_cell_volume = height_map_periph[np.where(cell_mask[:] == 1)].sum()

        if (len(self.get_multiple_nucleus_centroid())) > 1:
            return peripheral_cell_volume / len(self.get_multiple_nucleus_centroid()) * helpers.volume_coeff()
        else:
            return peripheral_cell_volume * helpers.volume_coeff()


class Image3dMultiNucleusWithSpots(Image3dMultiNucleus, Image3dWithSpots):
    """ Represents an image with identified spots (e.g. from FISH), has to have spots descriptor """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dMultiNucleus.is_a(repo, path) and Image3dWithSpots.is_a(repo, path)

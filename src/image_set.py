#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

from typing import List, Union

import numpy as np
import tqdm
from loguru import logger

import constants
from image import Image, ImageWithMTOC
from image3d import Image3d, Image3dWithSpots, Image3dWithIntensities, \
    Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC, Image3dWithSpotsAndIntensitiesAndMTOC, \
    Image3dMultiNucleus, Image3dMultiNucleusWithSpots
from imageWithIntensities import ImageWithIntensities, ImageWithIntensitiesAndMTOC
from imageWithSpots import ImageWithSpots, ImageWithSpotsAndMTOC
from imageWithZlines import imageMultiNucleusWithSpotsAndZlines, imageWithSpotsAndZlines, ImageMultiNucleus, \
    ImageMultiNucleusWithSpots, ImageWithSpotsAndIntensitiesAndMTOC, ImageWithSpotsAndIntensities
from repository import Repository


class ImageSet(object):
    def __init__(self, repository: Repository, path_list: List[str], force2D=False):
        """
        Sets instance variable path to the image location (image_path) in an HDF5 file object.
        """
        logger.info("Initializing image set from {}", path_list)
        self._repository = repository
        for p in path_list:
            if not self._repository.is_present(p):
                raise LookupError("Cannot init ImageSet, path %s is absent from repository" % p)
        self.path_list = path_list
        self.images = []
        all_paths = self._repository.get_image_path_list(self.path_list)
        for p in all_paths:
            has_multi_nucleus_spots_and_zlines = imageMultiNucleusWithSpotsAndZlines.is_a(self._repository, p)
            has_spots_and_zlines = imageWithSpotsAndZlines.is_a(self._repository, p)
            has_multi_nucleus = ImageMultiNucleus.is_a(self._repository, p)
            has_spots_multi_nucleus = ImageMultiNucleusWithSpots.is_a(self._repository, p)
            has_spots = ImageWithSpots.is_a(self._repository, p)
            has_intensities = ImageWithIntensities.is_a(self._repository, p)
            is_3d = Image3d.is_a(self._repository, p)
            is_3d_with_spots = Image3dWithSpots.is_a(self._repository, p)
            is_3d_with_intensities = Image3dWithIntensities.is_a(self._repository, p)
            is_3d_with_spots_and_MTOC = Image3dWithSpotsAndMTOC.is_a(self._repository, p)
            is_with_intensities_and_MTOC = ImageWithIntensitiesAndMTOC.is_a(self._repository, p)
            is_with_spots_and_intensities_and_MTOC = ImageWithSpotsAndIntensitiesAndMTOC.is_a(self._repository, p)
            is_3d_with_spots_and_intensities_and_MTOC = Image3dWithSpotsAndIntensitiesAndMTOC.is_a(self._repository, p)
            is_3d_with_intensities_and_MTOC = Image3dWithIntensitiesAndMTOC.is_a(self._repository, p)
            is_3d_multi_nucleus = Image3dMultiNucleus.is_a(self._repository, p)
            is_3d_spots_multi_nucleus = Image3dMultiNucleusWithSpots.is_a(self._repository, p)
            is_with_spots_and_MTOC = ImageWithSpotsAndMTOC.is_a(self._repository, p)
            # TODO : Refactor, test and debug the logical branching
            #logger.debug("creating Image for {} ",p)
            if force2D:
                if has_spots:
                    img = ImageWithSpots(self._repository, p)
                elif has_intensities:
                    img = ImageWithIntensities(self._repository, p)
                elif is_with_spots_and_MTOC:
                    img = ImageWithSpotsAndMTOC(self._repository, p)
                elif is_with_spots_and_MTOC:
                    img = ImageWithSpotsAndMTOC(self._repository, p)
                    if has_intensities:
                        img = ImageWithSpotsAndIntensitiesAndMTOC(self._repository, p)
                elif is_with_spots_and_intensities_and_MTOC:
                    img = ImageWithSpotsAndIntensitiesAndMTOC(self._repository, p)
                elif has_spots and has_intensities:
                    img = ImageWithSpotsAndIntensities(self._repository, p)
                elif has_spots_and_zlines:
                    img = imageWithSpotsAndZlines(self._repository, p)
                elif has_spots_multi_nucleus:
                    img = ImageMultiNucleusWithSpots(self._repository, p)
                elif has_multi_nucleus:
                    img = ImageMultiNucleus(self._repository, p)
                else:
                    raise NotImplemented("Couldn't deduce type for image %s" % p)
            else:
                if is_3d_spots_multi_nucleus:
                    img = Image3dMultiNucleusWithSpots(self._repository, p)
                elif is_3d_multi_nucleus:
                    img = Image3dMultiNucleus(self._repository, p)
                elif is_3d_with_spots_and_intensities_and_MTOC:
                    img = Image3dWithSpotsAndIntensitiesAndMTOC(self._repository, p)
                elif is_3d_with_spots_and_MTOC:
                    img = Image3dWithSpotsAndMTOC(self._repository, p)
                elif is_3d_with_intensities_and_MTOC:
                    img = Image3dWithIntensitiesAndMTOC(self._repository, p)
                elif is_3d_with_spots:
                    img = Image3dWithSpots(self._repository, p)
                elif is_3d_with_intensities:
                    img = Image3dWithIntensities(self._repository, p)
                elif is_3d:
                    img = Image3d(self._repository, p)
                elif is_with_spots_and_intensities_and_MTOC:
                    img = ImageWithSpotsAndIntensitiesAndMTOC(self._repository, p)
                elif has_multi_nucleus_spots_and_zlines:
                    img = imageMultiNucleusWithSpotsAndZlines(self._repository, p)
                elif has_spots_and_zlines:
                    img = imageWithSpotsAndZlines(self._repository, p)
                elif is_with_spots_and_MTOC:
                    img = ImageWithSpotsAndMTOC(self._repository, p)
                elif is_with_intensities_and_MTOC:
                    img = ImageWithIntensitiesAndMTOC(self._repository, p)
                elif has_spots and has_intensities:
                    img = ImageWithSpotsAndIntensities(self._repository, p)
                elif has_spots:
                    img = ImageWithSpots(self._repository, p)
                elif has_intensities:
                    img = ImageWithIntensities(self._repository, p)
                elif has_spots_multi_nucleus:
                    img = ImageMultiNucleusWithSpots(self._repository, p)
                elif has_multi_nucleus:
                    img = ImageMultiNucleus(self._repository, p)

                else:
                    raise NotImplemented("Couldn't deduce type for image %s" % p)

                # TODO: should we fix it? this is a hack
                if is_3d and np.sum(img.get_height_map()) < constants.analysis_config['MIN_HEIGHT_MAP_AREA']:
                    logger.debug("Image has bad height map {}", img._path)
                    continue
            # TODO: should we fix it? this is a hack
            if has_spots and len(img.get_cytoplasmic_spots()) < constants.analysis_config['MIN_SPOT_NUM']:
                logger.debug("Image contains too few spots {} - threshold = {}", img._path,
                             constants.analysis_config['MIN_SPOT_NUM'])
                continue

            if has_spots_multi_nucleus and len(img.get_cytoplasmic_spots()) > constants.analysis_config['MAX_SPOT_NUM']:
                logger.debug("Image contains too much spots ({}) - {} - threshold = {}", len(img.get_cytoplasmic_spots()), img._path,
                             constants.analysis_config['MAX_SPOT_NUM'])
                continue

            if has_intensities and img.signal_to_noise() < constants.analysis_config['MIN_SNR']:
                logger.debug("Insufficient signal to noise ratio for image {}", img._path)
                continue

            self.images.append(img)

        logger.info("Initialized image set from {} with {} images", path_list, len(self.images))

    def __sizeof__(self):
        return len(self.images)

    def get_images(self):
        return self.images

    def compute_spots_fractions_per_periphery(self):
        all_signals = self.compute_signal_from_periphery()
        image: ImageWithIntensities
        spot_counts = [len(image.get_spots()) for image in self.images]
        return np.array(all_signals) / np.array(spot_counts)[:, None]

    def compute_cytoplsamic_spots_fractions_per_periphery(self):
        all_signals = self.compute_signal_from_periphery()
        image: ImageWithSpots
        #cytoplasmic_densities = [image.compute_cytoplasmic_density() for image in self.images]
        #all_areas = self.compute_areas_from_periphery()
        spot_counts = [len(image.get_cytoplasmic_spots()) for image in self.images]
        return np.array(all_signals) / np.array(spot_counts)[:, None]
        #return np.divide(np.array(all_signals), np.array(all_areas)) / np.array(cytoplasmic_densities)[:, None]

    def compute_intensities_fractions_per_periphery(self):
        all_signals = self.compute_signal_from_periphery()
        image: ImageWithIntensities
        intensities_counts = [np.multiply(image.get_intensities(), image.get_cell_mask()).sum()
                              for image in self.images]
        return np.array(all_signals) / np.array(intensities_counts)[:, None]

    def compute_cytoplsamic_intensities_fractions_per_periphery(self):
        all_signals = self.compute_signal_from_periphery()
        image: ImageWithIntensities
        intensities_counts = [image.compute_cytoplasmic_intensities().sum() for image in self.images]
        return np.array(all_signals) / np.array(intensities_counts)[:, None]

    def compute_signal_from_periphery(self):
        arr = np.zeros((self.__sizeof__(), constants.analysis_config['NUM_CONTOURS']))
        image: Union[ImageWithSpots, ImageWithIntensities]
        for image_num, image in tqdm.tqdm(enumerate(self.images), desc="Images", total=self.__sizeof__()):
            arr[image_num] = image.compute_signal_from_periphery()
        return arr

    def compute_areas_from_periphery(self):
        arr = np.zeros((self.__sizeof__(), 100))
        image: Image
        for image_num, image in tqdm.tqdm(enumerate(self.images), desc="Images", total=self.__sizeof__()):
            arr[image_num] = image.compute_areas_from_periphery()
        return arr

    def compute_volumes_from_periphery(self):
        arr = np.zeros((self.__sizeof__(), 100))
        image: Image3d
        for image_num, image in tqdm.tqdm(enumerate(self.images), desc="Images", total=self.__sizeof__()):
            arr[image_num] = image.compute_volumes_from_periphery()
        return arr

    def compute_cytoplasmic_spots_counts(self) -> List[int]:
        cytoplasmic_spots_counts = []
        image: ImageWithSpots
        for image in self.images:
            cytoplasmic_spots_counts.append(image.compute_cytoplasmic_total_spots())
        return cytoplasmic_spots_counts

    def compute_cytoplasmic_intensities(self) -> List[float]:
        total_cytoplasmic_intensities = []
        image: ImageWithIntensities
        for image in self.images:
            total_cytoplasmic_intensities.append(image.compute_cytoplasmic_total_intensity())
        return total_cytoplasmic_intensities

    def compute_cytoplasmic_spots_centrality(self) -> List[float]:
        centralities = np.array([])
        image: Union[ImageWithSpots, Image3dWithSpots]
        for image in self.images:
            centralities= np.append(centralities, image.compute_spots_normalized_distance_to_centroid())
        valid_centralities = centralities[~np.isnan(centralities)]
        if len(valid_centralities) < len(centralities):
            logger.warning("spots out of hull for {} images out of {}",
                           len(centralities)-len(valid_centralities), self.__sizeof__())
        l = len(valid_centralities[valid_centralities > 1])
        if l > 0:
            logger.debug("normalized distance to centroid is > 1 for {} images out of {}",
                           l, self.__sizeof__())
        return valid_centralities

    def compute_cytoplasmic_spots_spread(self) -> List[float]:
        spots_spread = np.array([])
        image: Union[ImageWithSpots, Image3dWithSpots]
        for image in self.images:
            spots_spread = np.append(spots_spread, image.compute_spots_normalized_cytoplasmic_spread())
        l = len(spots_spread[spots_spread > 1])
        if l > 0:
            logger.debug("normalized distance to centroid is > 1 for {} images out of {}",
                           l, self.__sizeof__())
        return spots_spread

    def compute_intensities_cytoplasmic_centrality(self) -> List[float]:
        centralities = np.array([])
        image: Union[ImageWithIntensities, Image3dWithIntensities]
        for image in self.images:
            centralities = np.append(centralities, image.compute_intensities_normalized_spread_to_centroid())
        valid_centralities = centralities[~np.isnan(centralities)]
        if len(valid_centralities) < len(centralities):
            logger.warning("problematic intensity spread for {} images out of {}",
                           len(centralities) - len(valid_centralities), self.__sizeof__())
        l = len(centralities[centralities>1])
        if l > 0:
            logger.debug("normalized distance to centroid is > 1 for {} images out of {}",
                            l, self.__sizeof__())
        return valid_centralities

    def compute_intensities_cytoplasmic_spread(self) -> List[float]:
        spreads = np.array([])
        image: Union[ImageWithIntensities, Image3dWithIntensities]
        for image in self.images:
            spreads = np.append(spreads, image.compute_intensities_normalized_cytoplasmic_spread())
        l = len(spreads[spreads > 1])
        if l > 0:
            logger.debug("normalized distance to centroid is > 1 for {} images out of {}",
                           l, self.__sizeof__())
        return spreads

    def compute_degree_of_clustering(self):
        image: Union[ImageWithSpots, Image3dWithSpots]
        return [image.compute_degree_of_clustering() for image in self.images]

    def compute_normalised_quadrant_densities(self, quadrants_num=4, peripheral_flag=False,
                                              stripes=3, stripes_flag = False) -> np.array:
        """
        computes normalized densities per slice (quadrant by default) for all images
        """
        all_densities = np.empty((0, 2), float)
        for image in tqdm.tqdm(self.images, desc="Images"):
            cytoplasmic_density = image.compute_cytoplasmic_density()
            mdmq = image.get_or_compute_quadrant_densities(quadrants_num, peripheral_flag, stripes, stripes_flag)
            mdmq[:, 0] = mdmq[:, 0] / cytoplasmic_density
            all_densities = np.append(all_densities, mdmq[mdmq[:, 1].argsort()[::-1]], axis=0)

        mtoc_num = all_densities[all_densities[:,1]==1][:,1].sum()
        non_mtoc_num = len(all_densities[all_densities[:,1]==0][:,1])
        logger.debug("\nMTOC density {} and non MTOC density {} per element (quadrant or slice)",
                     all_densities[all_densities[:,1]==1][:,0].sum() / mtoc_num,
                     all_densities[all_densities[:,1]==0][:,0].sum() / non_mtoc_num)
        return all_densities

    def mtoc_is_in_leading_edge(self):
        image: ImageWithMTOC
        return [image.mtoc_is_in_leading_edge() for image in self.images]

    def compute_cytoplasmic_density(self):
        return [image.compute_cytoplasmic_density() for image in self.images]

    def compute_spots_peripheral_distance(self):
        return np.array([image.compute_spots_peripheral_distance() for image in self.images])

    def compute_zline_distance(self, z_line_spacing):
        image: imageMultiNucleusWithSpotsAndZlines
        total_profile = []
        image_counter = 0
        for i, image in enumerate(tqdm.tqdm(self.images, desc="Images")):
            z_line_distance_profile=image.get_minimal_z_line_distance(z_line_spacing)
            total_profile.append(z_line_distance_profile)
            image_counter += 1
        total_profile = np.array(total_profile).reshape((image_counter, z_line_spacing))

        return total_profile

    # Implements equation 17 (supplemental)of padovan-merhar et al. 2015
    def compute_volume_corrected_nm(self):
        #image: ImageMultiNucleusWithSpots
        cell_volume= [image.compute_cell_volume() for image in self.images]
        transcount = [len(image.get_spots()) for image in self.images]
        coeffs = np.polyfit(cell_volume, transcount, 1)
        a = coeffs[1]
        b = coeffs[0]
        expected_mrnas = [(a + (b * volume)) for volume in cell_volume]
        variance_expected_mrnas = np.std(expected_mrnas) ** 2
        variance_mrnas = np.std(transcount) ** 2
        exp_mrnas = (np.mean(transcount) ** 2)
        nm = (variance_mrnas - variance_expected_mrnas) / exp_mrnas
        return nm

    # Implements equation 17 (supplemental)of padovan-merhar et al. 2015
    def compute_surface_corrected_nm(self):
        image: ImageMultiNucleusWithSpots
        cell_surface = [image.compute_cell_area() for image in self.images]
        transcount = [len(image.get_spots()) for image in self.images]
        coeffs = np.polyfit(cell_surface, transcount, 1)
        a = coeffs[1]
        b = coeffs[0]
        expected_mrnas = [(a + (b * volume)) for volume in cell_surface]
        variance_expected_mrnas = np.std(expected_mrnas) ** 2
        variance_mrnas = np.std(transcount) ** 2
        exp_mrnas = (np.mean(transcount) ** 2)
        nm = (variance_mrnas - variance_expected_mrnas) / exp_mrnas
        return nm

    def compute_cell_mask_between_nucleus_centroid(self):
        image: imageMultiNucleusWithSpotsAndZlines
        nuc_dist = []
        nucs_pos = []
        cell_masks = []
        nucs_dist = []
        for im in self.images:
            cell_mask = im.get_cell_mask()
            nucleus_centroids = im.get_multiple_nucleus_centroid()
            nucleus_centroid = np.sort(nucleus_centroids, axis=0)
            im_mask = []
            nucs = []
            nuc_pos = []
            for nuc_n in range(len(nucleus_centroid) - 1):
                cell_mask_copy = cell_mask.copy()
                nuc_pos.append([nucleus_centroid[nuc_n][0], nucleus_centroid[nuc_n + 1][0]])
                nuc_dist.append(nucleus_centroid[nuc_n + 1][0] - nucleus_centroid[nuc_n][0])
                nucs.append(nucleus_centroid[nuc_n + 1][0] - nucleus_centroid[nuc_n][0])
                cell_mask_copy[:, 0:nucleus_centroid[nuc_n][0]] = 0
                cell_mask_copy[:, nucleus_centroid[nuc_n + 1][0]::] = 0
                im_mask.append(cell_mask_copy)
            nucs_dist.append(nucs)
            cell_masks.append(im_mask)
            nucs_pos.append(nuc_pos)

        return nuc_dist, nucs_dist, cell_masks, nucs_pos
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import sys
from typing import List
import tqdm
from loguru import logger
import numpy as np
import pandas as pd
import constants
import scipy.stats as stats
from repository import Repository
from image import ImageWithSpots, ImageWithIntensities, ImageWithMTOC, \
    ImageWithSpotsAndIntensities, ImageWithSpotsAndMTOC, ImageWithSpotsAndIntensitiesAndMTOC, \
    ImageWithIntensitiesAndMTOC, ImageMultiNucleus, ImageMultiNucleusWithSpots, imageWithSpotsAndZlines, imageMultiNucleusWithSpotsAndZlines
from image3d import Image3d, Image3dWithSpots, Image3dWithIntensities, \
    Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC, Image3dWithSpotsAndIntensitiesAndMTOC, Image3dMultiNucleus, Image3dMultiNucleusWithSpots


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
                if  is_3d_spots_multi_nucleus:
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

            # # TODO: should we fix it? this is a hack
            # if is_3d and np.sum(img.get_height_map()) < constants.analysis_config['MIN_HEIGHT_MAP_AREA']:
            #     logger.debug("Image has bad height map {}", img._path)
            #     continue

            if has_intensities and img.signal_to_noise() < constants.analysis_config['MIN_SNR']:
                logger.debug("Insufficient signal to noise ratio for image {}", img._path)
                continue

            self.images.append(img)
            #logger.debug("Image type : {} ", img.__class__.__name__)

        logger.info("Initialized image set from {} with {} images", path_list, len(self.images))

    def __sizeof__(self):
        return len(self.images)

    def get_images(self):
        return self.images

    # TODO : change this function's name !!!
    def compute_histogram_spots_peripheral_fraction(self):
        arr = []
        for image in tqdm.tqdm(self.images, desc="Images"):
            spots_peripheral_distance = image.get_spots_peripheral_distance()
            peripheral_profile = np.zeros(100)
            for i in range(0, 100):
                peripheral_profile[i] = float(len(np.where(spots_peripheral_distance <= (i + 1))[0])) / float(
                    len(spots_peripheral_distance))
            arr.append(peripheral_profile)
        return arr

    # TODO why here are counts and for intensities fractions ?
    # TODO why here it is not normalized and for intensities it is ?
    def compute_histogram_spots_peripheral_counts(self):
        arr = []
        for image in tqdm.tqdm(self.images, desc="Images"):
            spots_peripheral_distance = image.get_spots_peripheral_distance()
            arr.append(len(spots_peripheral_distance[spots_peripheral_distance <=
                                                     constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']]) / len(
                spots_peripheral_distance))
        return arr

    # TODO this only compute in 2D  ??????
    def compute_histogram_intensities_peripheral_fractions(self):
        arr = []
        image: Image3dWithIntensitiesAndMTOC
        for image in tqdm.tqdm(self.images, desc="Images"):
            cell_mask = image.get_cell_mask()
            intensities = image.get_intensities()
            # TODO this line below not present in V1, we have to add it to reproduce nocodazole protein peripheral fraction results
            # intensities = np.multiply(intensities, cell_mask)
            cytoplasmic_intensity = intensities[
                image.get_nucleus_mask() == 0].sum()  # TODO : should be = image.get_cytoplasmic_total_intensity()
            peripheral_intensities = image.compute_peripheral_intensities()
            peripheral_intensities = peripheral_intensities[constants.analysis_config[
                'PERIPHERAL_FRACTION_THRESHOLD']] / cytoplasmic_intensity
            arr.append(peripheral_intensities)
        return arr

    def compute_cytoplasmic_spots_counts(self) -> List[int]:
        """
        was previously compute_cytoplasmic_total
        :return:
        """
        cytoplasmic_spots_counts = []
        image: ImageWithSpots
        for image in self.images:
            cytoplasmic_spots_counts.append(image.get_cytoplasmic_total_spots())
        return cytoplasmic_spots_counts

    def compute_cytoplasmic_intensities(self) -> List[float]:
        """
        was previously compute_cytoplasmic_total
        :return:
        """
        total_cytoplasmic_intensities = []
        image: ImageWithIntensities
        for image in self.images:
            total_cytoplasmic_intensities.append(image.get_cytoplasmic_total_intensity())
        return total_cytoplasmic_intensities

    def compute_spots_cytoplasmic_spread(self) -> List[float]:
        cytoplasmic_spread_list = []
        image: Image3dWithSpots
        for image in self.images:
            cytoplasmic_spread_list.append(image.compute_spots_cytoplasmic_spread())
        return cytoplasmic_spread_list

    def compute_intensities_cytoplasmic_spread(self) -> List[float]:
        cytoplasmic_spread_list = []
        image: Image3dWithIntensities
        for image in self.images:
            cytoplasmic_spread_list.append(image.compute_intensities_cytoplasmic_spread())
        return cytoplasmic_spread_list

    def compute_degree_of_clustering(self):
        image: Image3dWithSpots
        return [image.compute_degree_of_clustering() for image in self.images]

    def compute_mtoc_dependent_degree_of_clustering(self):
        image: Image3dWithSpotsAndMTOC
        return [image.compute_degree_of_clustering() for image in self.images if image.mtoc_is_in_leading_edge()]

    def compute_normalised_quadrant_densities(self, quadrant_labels: list,
                                              mtoc_quadrant_label='MTOC', quadrants_num=4) -> dict:
        """
        builds a dictionary of densities per slice (quadrant by default) for all images
        The dictionary's keys are the quadrant_labels + mtoc_quadrant_label
        """
        if quadrants_num != len(quadrant_labels) + 1:
            raise RuntimeError("Quandrants number quadrants_num and total labels number have to be the same")

        gene_dict = {key: list([]) for key in quadrant_labels + [mtoc_quadrant_label]}
        mtoc_count = 0
        for image in tqdm.tqdm(self.images, desc="Images"):
            cytoplasmic_density = image.compute_cytoplasmic_density()
            mdmq = image.get_quadrants_densities(quadrants_num)
            mdmq[:, 0] = mdmq[:, 0] / cytoplasmic_density
            # add MTOC density (flag 1 in mdmq)
            if (mdmq[:, 0][np.where(mdmq[:, 1] == 1)][0] > np.median(mdmq[:, 0][np.where(mdmq[:, 1] == 0)])):
                mtoc_count = mtoc_count + 1
            gene_dict[mtoc_quadrant_label].append(mdmq[:, 0][np.where(mdmq[:, 1] == 1)][0])
            # add non MTOC densities
            for label, val in zip(quadrant_labels, mdmq[:, 0][np.where(mdmq[:, 1] == 0)]):
                gene_dict[label].append(val)

        logger.debug("\nMTOC enrichment {} out of {} ; relative {} \n", mtoc_count, len(self.images),
                     mtoc_count / len(self.images))
        return gene_dict

    def compute_peripheral_normalised_quadrant_densities(self, quadrant_labels: list,
                                                         mtoc_quadrant_label='MTOC', quadrants_num=4) -> dict:
        """
        builds a dictionary of densities per slice (quadrant by default) for all images
        The dictionary's keys are the quadrant_labels + mtoc_quadrant_label
        """
        if quadrants_num != len(quadrant_labels) + 1:
            raise RuntimeError("Quandrants number quadrants_num and total labels number have to be the same")

        gene_dict = {key: list([]) for key in quadrant_labels + [mtoc_quadrant_label]}
        mtoc_count = 0
        for image in tqdm.tqdm(self.images, desc="Images"):
            try:
                peripheral_fraction_threshold = constants.analysis_config['PERIPHERAL_FRACTION_THRESHOLD']
                if "/mrna/pard3/" in image._path:
                    peripheral_fraction_threshold = 50;
                cytoplasmic_density = image.compute_cytoplasmic_density()
                mdmq = image.get_peripheral_quadrants_densities(quadrants_num,
                                                                peripheral_fraction_threshold=peripheral_fraction_threshold)
                mdmq[:, 0] = mdmq[:, 0] / cytoplasmic_density
                # add MTOC density (flag 1 in mdmq)
                if (mdmq[:, 0][np.where(mdmq[:, 1] == 1)][0] > np.median(mdmq[:, 0][np.where(mdmq[:, 1] == 0)])):
                    mtoc_count = mtoc_count + 1
                gene_dict[mtoc_quadrant_label].append(mdmq[:, 0][np.where(mdmq[:, 1] == 1)][0])
                # add non MTOC densities
                for label, val in zip(quadrant_labels, mdmq[:, 0][np.where(mdmq[:, 1] == 0)]):
                    gene_dict[label].append(val)
            except RuntimeError as rte:
                print(rte)

        logger.debug("\nMTOC enrichment {} out of {} ; relative {} \n", mtoc_count, len(self.images),
                     mtoc_count / len(self.images))
        return gene_dict

    def compute_normalized_quadrant_and_slice_densities(self, quadrants_num=4, stripes=3):
        """
        build an array of densities per slices (quadrant by default) for all images
        the array size is an n * m array where n is the number of images and m the number of slices.
        For each image, the first slice is the one that contains the MTOC.
        """
        arr = np.zeros((len(self.images), stripes * quadrants_num))
        for i, image in enumerate(tqdm.tqdm(self.images, desc="Images")):
            try:
                cytoplasmic_density = image.compute_cytoplasmic_density()
                arr[i, :] = image.get_quadrants_and_slices_densities(quadrants_num, stripes)
                arr[i, :] = arr[i, :] / cytoplasmic_density
            except RuntimeError as rte:
                print(rte)
        return arr

    def compute_peripheral_normalized_quadrant_and_slice_densities(self, quadrants_num=4, stripes=3):
        """
        build an array of densities per slices (quadrant by default) for all images
        the array size is an n * m array where n is the number of images and m the number of slices.
        For each image, the first slice is the one that contains the MTOC.
        """
        arr = np.zeros((len(self.images), stripes * quadrants_num))
        for i, image in enumerate(tqdm.tqdm(self.images, desc="Images")):
            try:
                peripheral_density=image.compute_peripheral_density()
                arr[i, :] = image.get_peripheral_quadrants_and_slices_densities(quadrants_num, stripes)
                arr[i, :] = arr[i, :] / peripheral_density
            except RuntimeError as rte:
                print(rte)
        return arr

    def mtoc_is_in_leading_edge(self):
        image: ImageWithMTOC
        return [image.mtoc_is_in_leading_edge() for image in self.images]

    def compute_cytoplasmic_density(self):
        return [image.compute_cytoplasmic_density() for image in self.images]

    def compute_spots_peripheral_distance(self):
        return np.array([image.get_spots_peripheral_distance() for image in self.images])

    # def compute_spots_peripheral_distance_2D(self):
    #     image: ImageWithSpotsAndMTOC
    #     return np.array([image.get_spots_peripheral_distance() for image in self.images])

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
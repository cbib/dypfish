#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

'''Build a dataframe containing mrna and protein normalised distributions between MTOC quadrant,
   MTOC leading edge quadrant and non MTOC quadrants'''

import h5py
import numpy as np
import pandas as pd
import tqdm
import argparse
import src.helpers as helps
import src.image_descriptors as idsc
import src.path as path
from src.utils import enable_logger, check_dir, loadconfig

pd.set_option('display.max_rows', 500)

parser = argparse.ArgumentParser()
parser.add_argument("--peripheral", "-p", help='boolean flag: perform peripheral computation or not',
                    action="store_true", default=False)
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)

args = parser.parse_args()
is_periph = args.peripheral
input_dir_name = args.input_dir_name


def main():
    check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/')
    configData = loadconfig(input_dir_name)

    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    mtoc_file_name = configData["MTOC_FILE_NAME"]

    # Required descriptors: spots, IF, cell mask and height_map
    # annotations for the number of the quadrant containing the mtoc comes from the mtoc5.h5
    with h5py.File(path.data_dir + input_dir_name + '/' + basic_file_name, "r") as file_handler, \
            h5py.File(path.data_dir + input_dir_name + '/' + secondary_file_name, "r") as second_file_handler, \
            h5py.File(path.data_dir + input_dir_name + '/' + mtoc_file_name, "r") as mtoc_file_handler:
        compute_mrna_counts_per_quadrant(file_handler, is_periph, mtoc_file_handler, second_file_handler, configData)
        compute_protein_counts_per_quadrant(file_handler, is_periph, mtoc_file_handler, second_file_handler, configData)


def compute_mrna_counts_per_quadrant(file_handler, is_periph, mtoc_file_handler, second_file_handler, configData):
    # mrna part
    logger = enable_logger()
    molecule_type = ['/mrna']
    mrnas = configData["GENES"]
    timepoints = configData["TIMEPOINTS_MRNA"]
    peripheral_fraction_threshold = configData["PERIPHERAL_FRACTION_THRESHOLD"]
    image_width = configData["IMAGE_WIDTH"]
    image_height = configData["IMAGE_HEIGHT"]
    volume_coeff = configData["VOLUME_COEFF"]
    global_mtoc = []
    global_non_mtoc1 = []
    global_non_mtoc2 = []
    global_non_mtoc3 = []
    global_mtoc_leading = []
    global_mrna = []
    global_image = []
    global_index = []
    global_timepoint = []
    for mrna in tqdm.tqdm(mrnas, desc="Processing mRNAs"):
        logger.info("Computing for mRNA %s", mrna)
        for timepoint in tqdm.tqdm(timepoints, desc="Processing timepoints"):
            logger.info("Computing for mRNA %s : timepoint %s", mrna, timepoint)
            image_list = helps.preprocess_image_list3(file_handler, molecule_type, mrna, [timepoint])
            for image in tqdm.tqdm(image_list, desc="Processing images"):
                if is_periph:
                    spot_by_quad = idsc.search_periph_mrna_quadrants(file_handler, second_file_handler, peripheral_fraction_threshold, image_width, image_height, volume_coeff, image)
                else:
                    spot_by_quad = idsc.search_mrna_quadrants(file_handler, second_file_handler,image_width, image_height, volume_coeff, image)
                num_mtoc_quadrant = idsc.get_mtoc_quad(mtoc_file_handler, image)
                mtoc_spot = spot_by_quad[:, :, 1] == 1
                non_mtoc_spot = spot_by_quad[:, :, 1] == 0
                for i in range(90):
                    global_index.append(image.split("/")[4] + "_" + str(i + 1))
                    global_image.append(image.split("/")[4])
                    global_mrna.append(mrna)
                    global_timepoint.append(timepoint)
                global_mtoc.extend(spot_by_quad[mtoc_spot][:, 0].flatten())
                for i in range(0, 270, 3):
                    # TODO : has to be extend
                    global_non_mtoc1.append(spot_by_quad[non_mtoc_spot][:, 0].flatten()[i:i + 3][0])
                    global_non_mtoc2.append(spot_by_quad[non_mtoc_spot][:, 0].flatten()[i:i + 3][1])
                    global_non_mtoc3.append(spot_by_quad[non_mtoc_spot][:, 0].flatten()[i:i + 3][2])
                # quadrant number is 1 if it is in the leading edge
                if num_mtoc_quadrant == 1:
                    global_mtoc_leading.extend(spot_by_quad[mtoc_spot][:, 0].flatten())
                else:
                    for i in range(90):
                        global_mtoc_leading.append(np.nan)
    # TODO: add non_MTOC 1, 2 and 3
    df = pd.DataFrame(
        {'Image': global_image, 'Gene': global_mrna, 'timepoint': global_timepoint, 'MTOC': global_mtoc,
         'MTOC leading edge': global_mtoc_leading, 'Non MTOC1': global_non_mtoc1, 'Non MTOC2': global_non_mtoc2,
         'Non MTOC3': global_non_mtoc3}, index=global_index)
    if is_periph:
        df.to_csv(
            check_dir(
                path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_mrna_all_NMTOC.csv')
    else:
        df.to_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_mrna_all_NMTOC.csv')


def compute_protein_counts_per_quadrant(file_handler, is_periph, mtoc_file_handler, sec_file_handler, configData):
    # protein part
    logger = enable_logger()
    molecule_type = ['/protein']
    proteins = configData["PROTEINS"]
    timepoints = configData["TIMEPOINTS_PROTEIN"]
    peripheral_fraction_threshold = configData["PERIPHERAL_FRACTION_THRESHOLD"]
    image_width = configData["IMAGE_WIDTH"]
    image_height = configData["IMAGE_HEIGHT"]
    volume_coeff = configData["VOLUME_COEFF"]
    global_protein = []
    global_mtoc = []
    global_non_mtoc1 = []
    global_non_mtoc2 = []
    global_non_mtoc3 = []
    global_mtoc_leading = []
    global_image = []
    global_index = []
    global_timepoint = []
    for protein in tqdm.tqdm(proteins, desc="Processing proteins"):
        for timepoint in tqdm.tqdm(timepoints, desc="Processing timepoints"):
            image_list = helps.preprocess_image_list3(file_handler, molecule_type, protein, [timepoint])
            for image in tqdm.tqdm(image_list, desc="Processing images"):
                intensity_by_quad = \
                    idsc.search_periph_protein_quadrants(file_handler, sec_file_handler,peripheral_fraction_threshold, image_width, image_height, volume_coeff, image) if is_periph \
                        else \
                    idsc.search_protein_quadrants(file_handler, sec_file_handler,image_width, image_height, volume_coeff, image)
                mtoc_intensity = intensity_by_quad[:, :, 1] == 1
                non_mtoc_intensity = intensity_by_quad[:, :, 1] == 0
                num_mtoc_quadrant = idsc.get_mtoc_quad(mtoc_file_handler, image)

                for i in range(90):
                    # needed for final boxplot
                    global_index.append(image.split("/")[4] + "_" + str(i + 1))
                    global_image.append(image.split("/")[4])
                    global_protein.append(protein)
                    global_timepoint.append(timepoint)
                global_mtoc.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                for i in range(0, 270, 3):
                    # TODO: has to be extend
                    global_non_mtoc1.append(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 3][0])
                    global_non_mtoc2.append(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 3][1])
                    global_non_mtoc3.append(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 3][2])
                    # global_non_mtoc.append(np.mean(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 3]))
                if num_mtoc_quadrant == 1:
                    global_mtoc_leading.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                else:
                    for i in range(90):
                        global_mtoc_leading.append(np.nan)
    # TODO: add non_MTOC 1, 2 and 3
    df = pd.DataFrame({'Image': global_image, 'Gene': global_protein, 'timepoint': global_timepoint,
                       'MTOC': global_mtoc, 'MTOC leading edge': global_mtoc_leading, 'Non MTOC1': global_non_mtoc1,
                       'Non MTOC2': global_non_mtoc2, 'Non MTOC3': global_non_mtoc3},
                      index=global_index)
    if is_periph:
        df.to_csv(
            check_dir(
                path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_protein_all_NMTOC.csv')
    else:
        df.to_csv(
            check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_protein_all_NMTOC.csv')


if __name__ == "__main__":
    main()


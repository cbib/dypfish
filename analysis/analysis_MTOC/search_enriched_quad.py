#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import numpy as np
import h5py
import pandas as pd
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, check_dir

pd.set_option('display.max_rows', 500)

'''build dataframe containing mrna and protein normalised distributions between MTOC quadrant, 
   MTOC leading edge quadrant and non MTOC quadrant'''


def main(is_periph=False):
    enable_logger()
    check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/')

    # Required descriptors: spots, IF, cell mask an height_map
    with h5py.File(path.basic_file_path, "r") as file_handler, \
            h5py.File(path.secondary_file_path, "r") as sec_file_handler, \
            h5py.File(path.mtoc_file_path, "r") as mtoc_file_handler:
        molecule_type = ['/mrna']
        mrnas = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]

        # mrna part
        global_mtoc = []
        global_nmtoc = []
        global_mtoc_leading = []
        global_mrna = []
        global_image = []
        global_index = []
        global_timepoint = []
        for mrna in mrnas:
            print(mrna)

            timepoints = ["2h", "3h", "4h", "5h"]

            for timepoint in timepoints:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, mrna, [timepoint])

                for image in image_list:

                    spot_by_quad = idsc.search_periph_mrna_quadrants(file_handler, sec_file_handler, image) \
                        if is_periph else idsc.search_mrna_quadrants(file_handler, image)
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
                    mtoc_spot = spot_by_quad[:, :, 1] == 1
                    non_mtoc_spot = spot_by_quad[:, :, 1] == 0
                    for i in range(90):
                        global_index.append(image.split("/")[4] + "_" + str(i + 1))
                        global_image.append(image.split("/")[4])
                        global_mrna.append(mrna)
                        global_timepoint.append(timepoint)

                    global_mtoc.extend(spot_by_quad[mtoc_spot][:, 0].flatten())
                    for i in range(0, 270, 3):
                        global_nmtoc.append(np.mean(spot_by_quad[non_mtoc_spot][:, 0].flatten()[i:i + 2]))

                    if mtoc_quad_j == 1:
                        global_mtoc_leading.extend(spot_by_quad[mtoc_spot][:, 0].flatten())
                    else:
                        for i in range(90):
                            global_mtoc_leading.append(np.nan)

        df = pd.DataFrame(
            {'Image': global_image, 'Gene': global_mrna, 'timepoint': global_timepoint, 'Non MTOC': global_nmtoc,
             'MTOC': global_mtoc, 'MTOC leading edge': global_mtoc_leading}, index=global_index)

        if is_periph:

            df.to_csv(
                check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_mrna.csv')
        else:
            df.to_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_mrna.csv')
        # protein part
        molecule_type = ['/protein']
        proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

        global_protein = []
        global_mtoc = []
        global_nmtoc = []
        global_mtoc_leading = []
        global_image = []
        global_index = []
        global_timepoint = []
        for protein in proteins:

            ##Normal conditions
            timepoints = ["2h", "3h", "5h", "7h"]

            for timepoint in timepoints:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, protein, [timepoint])

                for image in image_list:

                    intensity_by_quad = idsc.search_periph_protein_quadrants(
                        file_handler, sec_file_handler, protein, image, path.path_data) if is_periph else \
                        idsc.search_protein_quadrants(file_handler, mtoc_file_handler, protein, image)

                    mtoc_intensity = intensity_by_quad[:, :, 1] == 1
                    non_mtoc_intensity = intensity_by_quad[:, :, 1] == 0
                    mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)

                    for i in range(90):
                        # needed for final boxplot
                        global_index.append(image.split("/")[4] + "_" + str(i + 1))
                        global_image.append(image.split("/")[4])
                        global_protein.append(protein)
                        global_timepoint.append(timepoint)

                    global_mtoc.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                    for i in range(0, 270, 3):
                        global_nmtoc.append(np.mean(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 2]))

                    if mtoc_quad_j == 1:
                        global_mtoc_leading.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
                    else:
                        for i in range(90):
                            global_mtoc_leading.append(np.nan)

        df = pd.DataFrame({'Image': global_image, 'Gene': global_protein, 'timepoint': global_timepoint,
                           'Non MTOC': global_nmtoc, 'MTOC': global_mtoc, 'MTOC leading edge': global_mtoc_leading},
                          index=global_index)

        if is_periph:

            df.to_csv(
                check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_protein.csv')
        else:
            df.to_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_protein.csv')


if __name__ == "__main__":
    # main()
    main(is_periph=True)

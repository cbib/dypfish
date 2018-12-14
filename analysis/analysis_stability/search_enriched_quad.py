#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pandas as pd

pd.set_option('display.max_rows', 500)

import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, check_dir


def plot_bar_profile(data, genes, y_limit, ylabel, figname, colors):
    ## third technic
    fig = plt.figure()
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ## the data
    N = len(genes)
    y_lim_max = np.max(data) + 0.2
    y_lim_min = np.min(data) - 0.2

    ## necessary variables
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars

    ## the bars
    rects1 = ax.bar(ind, data, width, color=colors)
    # axes and labels
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(y_lim_min, y_lim_max)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    ax.set_xticks(ind)
    plt.legend([gene for gene in genes], loc='upper right')
    ax.legend(rects1, genes, prop={'size': 8})
    plt.savefig(figname, format='svg')
    plt.show()


def main():
    enable_logger()

    # Required descriptors: spots, IF, cell mask an height_map
    with h5py.File(path.basic_file_path, "r") as file_handler, h5py.File(
            path.mtoc_file_path, "r") as mtoc_file_handler:
        molecule_type = ['/mrna']
        mrnas = ["arhgdia", "arhgdia_cultured"]
        # mrnas = [ "arhgdia_nocodazole","pard3_nocodazole"]

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
            spots_mtoc_all = []
            spots_non_mtoc_all = []

            timepoints = ["3h"]
            # timepoints = ["3h", "5h"]

            for timepoint in timepoints:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, mrna, [timepoint])

                for image in image_list:
                    # if image != '/mrna/arhgdia/4h/46':
                    #     continue
                    spot_by_quad = idsc.search_mrna_quadrants(file_handler, image)
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

        df.to_csv(check_dir(path.analysis_dir + 'analysis_stability/dataframe/') + 'global_mtoc_file_all_mrna.csv')

        # # protein part
        # molecule_type = ['/protein']
        # proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        # #proteins = [ "pard3"]
        #
        # global_protein = []
        # global_mtoc = []
        # global_nmtoc = []
        # global_mtoc_leading = []
        # global_image = []
        # global_index = []
        # global_timepoint = []
        # for protein in proteins:
        #
        #     ##Normal conditions
        #     timepoints = ["2h", "3h", "5h", "7h"]
        #     ##nocodazole condition
        #     #timepoints = ["3h"]
        #
        #     for timepoint in timepoints:
        #         image_list = helps.preprocess_image_list3(file_handler, molecule_type, protein, [timepoint])
        #
        #         for image in image_list:
        #             # if image != '/protein/beta_actin/5h/15':
        #             #     continue
        #
        #
        #             intensity_by_quad = idsc.search_protein_quadrants(file_handler, mtoc_file_handler, protein, image,
        #                                                               path.path_data)
        #             mtoc_intensity = intensity_by_quad[:, :, 1] == 1
        #             non_mtoc_intensity = intensity_by_quad[:, :, 1] == 0
        #             mtoc_quad_j = idsc.get_mtoc_quad(mtoc_file_handler, image)
        #
        #             for i in range(90):
        #                 # needed for final boxplot
        #                 global_index.append(image.split("/")[4] + "_" + str(i + 1))
        #                 global_image.append(image.split("/")[4])
        #                 global_protein.append(protein)
        #                 global_timepoint.append(timepoint)
        #
        #             global_mtoc.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
        #             for i in range(0, 270, 3):
        #                 global_nmtoc.append(np.mean(intensity_by_quad[non_mtoc_intensity][:, 0].flatten()[i:i + 2]))
        #
        #             if mtoc_quad_j == 1:
        #                 global_mtoc_leading.extend(intensity_by_quad[mtoc_intensity][:, 0].flatten())
        #             else:
        #                 for i in range(90):
        #                     global_mtoc_leading.append(np.nan)
        #
        # df = pd.DataFrame({'Image': global_image, 'Gene': global_protein, 'timepoint': global_timepoint,
        #                    'Non MTOC': global_nmtoc, 'MTOC': global_mtoc, 'MTOC leading edge': global_mtoc_leading},
        #                   index=global_index)
        #
        # df.to_csv(check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_protein.csv')


if __name__ == "__main__":
    main()

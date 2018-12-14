#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
import os
import os.path
import src.plot as plot

import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import check_dir

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
# if "log" not in globals():
# logger = Logger.init_logger('REFORMAT_%s'%(cfg.language_code), load_config())
logger.info("Running %s", sys.argv[0])


def cytoplasmic_spread(file_handler, molecule_type,genes,path_data,colors,y_lim):
    cyt_spreads = []
    #dfname = path.analysis_dir + "analysis_cytoplasmic_spread/df/"+ molecule_type[0] +"_cyt_spread.csv"
    figname = check_dir(path.analysis_dir + '/analysis_nocodazole/figures/cyt_spread/') + molecule_type[0].split('/')[1] + '_cytoplasmic_spread.png'

    for gene in genes:
        image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread(image_list, file_handler, path_data))
    stan.plot_bar_profile(cyt_spreads, genes, y_lim, 'Cytoplasmic spread', figname, colors)

    # if os.path.isfile(path.analysis_dir+"analysis_cytoplasmic_spread/df/mrna_cyt_spread.csv"):
    #     df=pd.read_csv(dfname)
    #     stan.plot_bar_profile(cyt_spreads, genes, 1.60, 'Cytoplasmic spread', figname, colors)
    # else:
    #     dict={}
    #     for gene in genes:
    #
    #         image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
    #         cyt_spreads.append(adsc.compute_cytoplasmic_spread(file_handler, image_list, path_data))
    #         dict[gene]=adsc.compute_cytoplasmic_spread(file_handler, image_list, path_data)
    #     df = pd.DataFrame(dict)
    #     df.to_csv(dfname)


if __name__ == "__main__":

    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map

    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    ''' 
    1-You need to create a password.txt file before running to connect via ssh
    '''



    basic_file_basename = 'basic'
    sec_file_basename = 'secondary'
    path_data = path.raw_data_dir

    # basic_file_basename = 'basics_scratch_data'
    # sec_file_basename = 'secondary_scratch_data'
    # path_data=path.scratch_data_dir


    basic_file_path = path.analysis_data_dir + basic_file_basename + '.h5'
    secondary_file_path = path.analysis_data_dir + sec_file_basename + '.h5'

    # Compute bar plot cytoplasmic spread

    # micropatterned data

    # colors = ['#1E95bb', '#1ec5d4']# colors for arhgdia
    colors = ['#F16c1b', '#f1bc1b']# colors for pard3


    # genes = ["arhgdia", "arhgdia_nocodazole"]
    genes = ["pard3", "pard3_nocodazole"]

    # proteins = ["arhgdia", "arhgdia_nocodazole", "arhgdia_cytod", "pard3", "pard3_nocodazole","pard3_cytod"]
    # proteins = ["arhgdia", "arhgdia_nocodazole"]
    # proteins = ["pard3", "pard3_nocodazole"]





    timepoints = [ "3h", "5h"]
    timepoints_protein = ["3h",  "5h"]
    cell_type = 'micropatterned/'

    # scratch data
    # genes = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    # proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    # timepoints = ["1h", "3h", "5h"]
    # timepoints_protein = ["1h", "3h",'7h']
    # colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B']
    # cell_type='cultured/'



    with h5py.File(basic_file_path, "r") as file_handler:


        molecule_type = ['/mrna']
        cytoplasmic_spread(file_handler,molecule_type,genes,path_data,colors,0.4)

        molecule_type = ['/protein']
        cytoplasmic_spread(file_handler,molecule_type,genes,path_data,colors,0.9)

    # # This part produces plot interpolation of cytoplasmic spread by timepoint
    # with h5py.File(basic_file_path, "r") as file_handler:
    #
    #     molecule_type = ['/mrna']
    #     genes = ["arhgdia_nocodazole", "arhgdia", "pard3_nocodazole", "pard3"]
    #     proteins = ["arhgdia_nocodazole", "arhgdia", "pard3_nocodazole", "pard3"]
    #     timepoints = ["3h", "5h"]
    #     timepoints_protein = ["3h", "5h"]
    #     #colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
    #     colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']
    #
    #     for i in range(len(genes)):
    #         fig = plt.figures()
    #         ax = plt.axes()
    #         ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    #         ax.set_xlim(1, 8)
    #         ax.set_ylim(0.75, 1.40)
    #         x_mrna = np.arange(3, 5, 0.01)
    #         mrna_data = np.zeros((3, 2))
    #         counter = 0
    #         timepoints = ["3h", "5h"]
    #         for timepoint in timepoints:
    #             print(genes[i], '_', timepoint)
    #             image_list = helps.preprocess_image_list3(file_handler, molecule_type, genes[i], [timepoint])
    #             cyt_spread = adsc.compute_cytoplasmic_spread(file_handler, image_list,path_data)
    #             cyt_spread_median = np.median(cyt_spread)
    #             print(cyt_spread_median)
    #             err = np.median(np.abs(np.tile(np.median(cyt_spread), (1, len(cyt_spread))) - cyt_spread))
    #             upp_env = cyt_spread_median + err
    #             low_env = cyt_spread_median - err
    #             mrna_data[0, counter] = cyt_spread_median
    #             mrna_data[1, counter] = upp_env
    #             mrna_data[2, counter] = low_env
    #             counter += 1
    #         mrna_data = mrna_data / np.mean(mrna_data[0, :])
    #
    #         print(mrna_data[0, :],mrna_data[1, :],mrna_data[2, :])
    #         spl = interpolate.UnivariateSpline([3, 5], mrna_data[0, :],k=1)
    #         spl_upp = interpolate.UnivariateSpline([3, 5], mrna_data[1, :],k=1)
    #         spl_low = interpolate.UnivariateSpline([3, 5], mrna_data[2, :],k=1)
    #         m_y_new = spl(x_mrna)
    #         m_y_new_upp = spl_upp(x_mrna)
    #         m_y_new_down = spl_low(x_mrna)
    #
    #         solid_mrna, = plt.plot(x_mrna, m_y_new, linestyle="-", color=colors[i], linewidth=2)
    #         plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=colors[i])
    #         plt.plot(x_mrna, m_y_new_down, linestyle="-", color=colors[i])
    #         ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=colors[i], alpha=0.5, interpolate=False)
    #
    #         if genes[i] in proteins:
    #
    #             counter = 0
    #             x_protein = np.arange(3, 5, 0.01)
    #             protein_data = np.zeros((3, 2))
    #             timepoints = ["3h", "5h"]
    #             for timepoint in timepoints:
    #                 print(genes[i], '_', timepoint)
    #                 image_list = helps.preprocess_image_list3(file_handler, ['/protein'], genes[i], [timepoint])
    #                 cyt_spread = adsc.compute_cytoplasmic_spread(file_handler, image_list,path_data)
    #                 cyt_spread_median = np.median(cyt_spread)
    #                 err = np.median(np.abs(np.tile(np.median(cyt_spread), (1, len(cyt_spread))) - cyt_spread))
    #                 upp_env = cyt_spread_median + err
    #                 low_env = cyt_spread_median - err
    #                 protein_data[0, counter] = cyt_spread_median
    #                 protein_data[1, counter] = upp_env
    #                 protein_data[2, counter] = low_env
    #                 counter += 1
    #             protein_data = protein_data / np.mean(protein_data[0, :])
    #
    #             spl = interpolate.UnivariateSpline([3, 5], protein_data[0, :],k=1)
    #             spl_upp = interpolate.UnivariateSpline([3, 5], protein_data[1, :],k=1)
    #             spl_low = interpolate.UnivariateSpline([3, 5], protein_data[2, :],k=1)
    #             p_y_new = spl(x_protein)
    #             p_y_new_upp = spl_upp(x_protein)
    #             p_y_new_down = spl_low(x_protein)
    #
    #             dashed_protein, = plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=colors[i],linewidth=2)
    #             plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=colors[i])
    #             plt.plot(x_protein, p_y_new_down, linestyle="--", color=colors[i])
    #             ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=colors[i], alpha=0.25,interpolate=False)
    #             plt.legend([dashed_protein, solid_mrna], ['Protein', 'Mrna'])
    #         else:
    #             plt.legend([solid_mrna], ['Mrna'])
    #         ax.set_ylabel('Cytoplasmic spread')
    #         ax.set_xlabel('Time(hrs)')
    #         ax.set_title(genes[i])
    #         plt.savefig(path.analysis_dir + '/analysis_nocodazole/figures/cyt_spread_' + genes[i]+'.png',format='png')
    #         plt.show()
    #
    #

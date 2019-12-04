#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import h5py
import math
import argparse
import json
import numpy as np
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
import src.plot as plot
from src.utils import enable_logger, plot_colors, check_dir, loadconfig

''' 
2-This script is supposed to be ran after compute h_star in dypfish. It produces:
 -bar plot for degree of clustering 
'''
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def main():
    # Required h5 file:  hstar.h5
    # Required descriptors:  h_star
    enable_logger()
    configData = loadconfig(input_dir_name)
    molecule_types = configData["MOLECULE_TYPES"]
    genes = configData["GENES"]
    timepoints = configData["TIMEPOINTS_MRNA"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    h_star_file_name = configData["HSTAR_FILE_NAME"]

    # produce bar plot for degree of clustering
    with h5py.File(path.data_dir+input_dir_name+'/'+h_star_file_name, "a") as input_file_handler:

        for molecule_type in molecule_types:
            print(molecule_type)

            base = math.log(0.5)
            global_median = []
            global_err = []
            for gene in genes:
                image_list = helps.preprocess_image_list2(input_file_handler, molecule_type, gene)
                dof = adsc.compute_degree_of_clustering_wo_MTOC(image_list,input_file_handler)
                global_median.append(math.log(np.median(dof)) - base)
                err = np.median(np.abs(np.tile(np.median(dof), (1, len(dof))) - dof))
                global_err.append(math.log(np.median(dof) + err) - math.log(np.median(dof)) - base)

            # plot figures
            figname = check_dir(path.analysis_dir + 'analysis_chx/figures/') + molecule_type+ '_degree_of_clustering.png'
            plot.bar_profile_median(global_median, genes, molecule_type, figname, global_err)

        # produce plot interpolation of degree of clustering by timepoint
        #degree_of_clustering_dynamic_profile(input_file_handler, genes, molecule_type,timepoints,timepoints_num_mrna,timepoints_num_protein)

        data_generator = plot.data_extractor_generic(genes, genes, timepoints, timepoints, input_file_handler, adsc.compute_degree_of_clustering_wo_MTOC, input_file_handler)
        for mrna_data, protein_data, i in data_generator:
            figpath = check_dir(path.analysis_dir + 'analysis_chx/figures/') + '/degree_of_clustering_' + genes[i] + '.png'
            plot.dynamic_profiles(mrna_data, protein_data, timepoints_num_mrna, timepoints_num_protein, genes[i], plot_colors[i], 'Time(hrs)', 'Degree of clustering(*)',figpath)


if __name__ == "__main__":
    main()

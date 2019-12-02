#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import h5py
import math
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


def degree_of_clustering_dynamic_profile(input_file_handler,genes, molecule_type, timepoints, timepoints_num_mrna, timepoints_num_protein):

        data_generator = plot.data_extractor_generic(genes, genes, timepoints, input_file_handler, adsc.compute_degree_of_clustering_wo_MTOC, input_file_handler)
        #data_generator = plot.data_extractor(genes, genes, input_file_handler, adsc.compute_degree_of_clustering_wo_MTOC, input_file_handler )
        print(data_generator)
        for mrna_data, protein_data, i in data_generator:
            figpath = check_dir(path.analysis_dir + 'analysis_chx/figures/') + '/dof_' + genes[i] + '.png'
            plot.dynamic_profiles(mrna_data, protein_data, timepoints_num_mrna, timepoints_num_protein, genes[i], plot_colors[i], 'Time(hrs)', 'Degree of clustering(*)',figpath)


def main():
    # Required descriptors: cell_mask, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()
    configData = loadconfig("chx")
    molecule_types = configData["MOLECULE_TYPES"]
    genes = configData["GENES"]
    timepoints = configData["TIMEPOINTS"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]


    # produce bar plot for degree of clustering
    with h5py.File(path.h_star_chx_file_path, "a") as input_file_handler:

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
        degree_of_clustering_dynamic_profile(input_file_handler, genes, molecule_type,timepoints,timepoints_num_mrna,timepoints_num_protein)



if __name__ == "__main__":
    main()

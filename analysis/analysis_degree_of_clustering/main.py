#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import matplotlib.pyplot as plt
import h5py
import math
import numpy as np
import sys
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
import src.plot as plot
from scipy import interpolate
from src.utils import enable_logger, plot_colors, check_dir

''' 
2-This script is supposed to be ran after compute h_star in dypfish. It produces:
 -bar plot for degree of clustering 
'''


def degree_of_clustering_dynamic_profile(input_file_handler,mtoc_file_handler,genes):

        data_generator = plot.data_extractor(genes, genes, input_file_handler,
                                             adsc.compute_degree_of_clustering, input_file_handler, mtoc_file_handler)
        print(data_generator)
    # TODO bad bad bad arguments to data_extractor function
        for mrna_data, protein_data, i in data_generator:
            figpath = check_dir(path.analysis_dir + 'analysis_degree_of_clustering/figures/') + '/dof_' + genes[i] + '.png'
            plot.dynamic_profiles(mrna_data, protein_data, genes[i], plot_colors[i], 'Time(hrs)', 'Degree of clustering(*)',figpath)


def main():
    # Required descriptors: cell_mask, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    enable_logger()

    # produce bar plot for degree of clustering
    with h5py.File(path.h_star_file_path, "a") as input_file_handler, h5py.File(path.mtoc_file_path, "a") as mtoc_file_handler:
        # mrna part
        molecule_type = ['/mrna']
        genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
        base = math.log(0.5)
        mrna_median = []
        mrna_err = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(input_file_handler, molecule_type[0], gene)
            dof = adsc.compute_degree_of_clustering(image_list,input_file_handler, mtoc_file_handler)
            mrna_median.append(math.log(np.median(dof)) - base)
            err = np.median(np.abs(np.tile(np.median(dof), (1, len(dof))) - dof))
            mrna_err.append(math.log(np.median(dof) + err) - math.log(np.median(dof)) - base)

        # protein part
        molecule_type = ['/protein']
        proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
        base = math.log(0.01)
        protein_median = []
        protein_err = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(input_file_handler, molecule_type[0], protein)
            dof = adsc.compute_degree_of_clustering(image_list,input_file_handler, mtoc_file_handler)
            protein_median.append(math.log(np.median(dof)) - base)
            err = np.median(np.abs(np.tile(np.median(dof), (1, len(dof))) - dof))
            protein_err.append(math.log(np.median(dof) + err) - math.log(np.median(dof)) - base)

        #plot figures
        figname = check_dir(path.analysis_dir + 'analysis_degree_of_clustering/figures/') + 'mrna_degree_of_clustering.png'

        plot.bar_profile_median(mrna_median, genes, 'mrna', figname, mrna_err)

        figname = check_dir(path.analysis_dir + 'analysis_degree_of_clustering/figures/') + 'protein_degree_of_clustering.png'
        plot.bar_profile_median(protein_median, proteins, 'protein', figname, protein_err)

        # produce plot interpolation of degree of clustering by timepoint
        degree_of_clustering_dynamic_profile(input_file_handler,mtoc_file_handler,proteins)


if __name__ == "__main__":
    main()

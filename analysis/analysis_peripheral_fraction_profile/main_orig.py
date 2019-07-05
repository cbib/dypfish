#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import h5py
import src.plot as plot
import src.constants as constants
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, plot_colors, check_dir


def peripheral_profile(mean_profiles, timepoint, genes, colors):
    timepoint = timepoint if timepoint else 'All'
    plot_filename = 'peripheral_fraction_' + timepoint + '_timepoint' + str(constants.NUM_CONTOURS) + 'contours.png'
    mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
    figure_title = ' peripheral profile ('+ timepoint+')'
    plot.profile(mean_profiles, genes, constants.NUM_CONTOURS,check_dir(path.analysis_dir + 'analysis_peripheral_fraction_profile/figures/')+ plot_filename, figure_title, colors, True)


def peripheral_profiles(file_handler, molecule_type, genes, colors, compute_peripheral_fraction_profiles, timepoints=[False]):
    for timepoint in timepoints:
        print(timepoint)
        mean_profiles = []
        for gene in genes:
            if timepoint:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, gene, [timepoint])
            else:
                image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            profiles=compute_peripheral_fraction_profiles(file_handler, image_list)
            mean_profiles.append(np.average(profiles, axis=0))
        peripheral_profile(mean_profiles, timepoint, genes, colors)

def peripheral_fraction_profile(sec_file_handler,molecule_type,genes,fraction,colors,file_handler):
    fractions = []
    for gene in genes:
        print(gene)
        image_list = helps.preprocess_image_list2(sec_file_handler, molecule_type[0], gene)
        fractions.append(adsc.build_histogram_periph_fraction(sec_file_handler, image_list, fraction,path.path_data,file_handler))
    print(fractions)
    figname = path.analysis_dir + 'analysis_peripheral_fraction_profile/figures/peripheral_fraction_'+str(fraction)+'.png'
    #plot.fraction_profile(fractions, fraction, genes, figname,colors)
    plot.bar_profile(fractions, genes, figname)

def histogram_peripheral_profile(basic_file_handler,secondary_file_handler,genes, proteins, colors, path_data):
    molecule_type = ['/mrna']
    periph_fraction = []
    for gene in genes:
        print(gene)
        image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], gene)
        periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
    figname = path.analysis_dir + 'analysis_peripheral_fraction_profile/figures/'+molecule_type[
        0] +'_peripheral_fraction.png'
    plot.bar_profile(periph_fraction, genes, figname)
    molecule_type = ['/protein']
    periph_fraction = []
    for protein in proteins:
        image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], protein)
        periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD,path_data))
    figname = path.analysis_dir + 'analysis_peripheral_fraction_profile/figures/'+molecule_type[0] +'_peripheral_fraction.png'
    plot.bar_profile(periph_fraction, proteins, figname)

def peripheral_fraction_dynamic_profile(basic_file_handler,secondary_file_handler, genes,proteins, colors, mrna_tp, protein_tp, path_data):
    data_generator = plot.data_extractor(genes, proteins, secondary_file_handler,
                                         adsc.compute_periph_fraction, basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD,path_data)
    for mrna_data, protein_data, i in data_generator:
        figpath = check_dir(
            path.analysis_dir + 'analysis_peripheral_fraction_profile/figures/') + '/peripheral_fraction_' + \
                  genes[i] + '.png'
        plot.dynamic_profiles(mrna_data, protein_data, genes[i], plot_colors[i], 'Time(hrs)', 'Peripheral fraction',
                              figpath)

def main():
    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots
    enable_logger()
    ## Build peripheral profile plot either for each or for all timepoint
    molecule_type = ['/mrna']
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    timepoints = ["2h", "3h", "4h", "5h"]
    timepoints_protein = ["2h", "3h", "5h", "7h"]
    with h5py.File(path.basic_file_path, "r") as basic_file_handler,h5py.File(path.secondary_file_path, "r") as secondary_file_handler:
        peripheral_profiles(secondary_file_handler, molecule_type, genes, plot_colors,adsc.compute_peripheral_fraction_profiles_3D)

        ## Section to build peripheral profile fraction 10 and 30
        # peripheral_fraction_profile(secondary_file_handler, molecule_type, genes, 10, plot_colors, basic_file_handler)

        ## Section to compute bar plot peripheral fraction
        # histogram_peripheral_profile(basic_file_handler, secondary_file_handler, genes, proteins, plot_colors, path.path_data)

        ## Section to produce plot interpolation (dynamic profile) of peripheral fraction by timepoint
        # peripheral_fraction_dynamic_profile(basic_file_handler, secondary_file_handler, genes, proteins, plot_colors, timepoints, timepoints_protein, path.path_data)

if __name__ == "__main__":
    main()

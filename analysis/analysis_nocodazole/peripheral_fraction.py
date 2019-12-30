#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import src.plot as plot
import src.constants as constants
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.helpers as helps
from src.utils import loadconfig,check_dir
import hashlib

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


def peripheral_fraction_profile(file_handler,molecule_type,genes,fraction,colors,image_type,gene_root_name,basic_file_handler):
    if len(image_type) == 0:
        fractions = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type, gene)
            if molecule_type=="mrna":
                fractions.append(adsc.build_histogram_mrna_periph_fraction(file_handler, image_list, fraction,path_data,basic_file_handler))
            else:
                fractions.append(adsc.build_histogram_protein_periph_fraction(file_handler, image_list, fraction,path_data,basic_file_handler))
        print(fractions)
        figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/')+molecule_type+'_peripheral_fraction_'+str(fraction)+'.png'
        #plot.fraction_profile(fractions, fraction, genes, figname,colors)
        plot.bar_profile(fractions,genes,figname)
    else:
        for image_t in image_type:
            fractions = []
            for gene in genes:
                image_list = helps.preprocess_image_list5(file_handler, molecule_type, gene,image_t)
                if molecule_type == "mrna":
                    fractions.append(adsc.build_histogram_mrna_periph_fraction(file_handler, image_list, fraction, path_data,basic_file_handler))
                else:
                    fractions.append(adsc.build_histogram_protein_periph_fraction(file_handler, image_list, fraction, path_data,basic_file_handler))
            figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/') +gene_root_name+'_peripheral_fraction_' +image_t+'_'+ str(fraction) + '.png'
            #plot.fraction_profile(fractions, fraction, genes, figname, colors)
            plot.bar_profile(fractions, genes, figname)


def histogram_peripheral_profile(basic_file_handler,secondary_file_handler, genes, proteins, colors, path_data,cell_type,image_type,periph_fraction_cst):
    if len(image_type) == 0:
        molecule_type = ['/mrna']
        periph_fraction = []
        for gene in genes:
            print(gene)
            image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], gene)
            periph_fraction.append(adsc.compute_mrna_periph_fraction(image_list,basic_file_handler, secondary_file_handler, periph_fraction_cst))
        print(periph_fraction)
        figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/')+molecule_type[
            0] +'_peripheral_fraction.png'
        plot.bar_profile(periph_fraction, genes, figname)
        molecule_type = ['/protein']
        periph_fraction = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], protein)
            periph_fraction.append(adsc.compute_protein_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  periph_fraction_cst))
        figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/') +molecule_type[0] +'_peripheral_fraction.png'
        plot.bar_profile(periph_fraction, proteins,  figname)
    else:
        for image_t in image_type:
            print(image_t)
            periph_fraction = []
            molecule_type = ['/mrna']
            for gene in genes:
                print(gene)
                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], gene,image_t)
                periph_fraction.append(
                    adsc.compute_mrna_periph_fraction(image_list,basic_file_handler, secondary_file_handler,periph_fraction_cst))
            figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/') + \
                      molecule_type[0] +'_'+image_t+ 'peripheral_fraction.svg'
            plot.bar_profile(periph_fraction, genes, figname)
            molecule_type = ['/protein']
            periph_fraction = []
            for protein in proteins:
                print(protein)
                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], protein,image_t)
                periph_fraction.append(adsc.compute_protein_periph_fraction(image_list,basic_file_handler, secondary_file_handler, periph_fraction_cst))
            figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/') + \
                      molecule_type[0] + '_' + image_t + '_peripheral_fraction.svg'
            plot.bar_profile(periph_fraction, proteins, figname)

if __name__ == "__main__":
    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots


    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"]
    proteins = configData["PROTEINS"]
    mrna_timepoints = configData["TIMEPOINTS_MRNA"]
    prot_timepoints = configData["TIMEPOINTS_PROTEIN"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    mtoc_file_name = configData["MTOC_FILE_NAME"]
    colors = configData["COLORS"]
    molecule_types=configData["MOLECULE_TYPES"]
    periph_fraction_cst = configData["PERIPHERAL_FRACTION_THRESHOLD"]
    cell_type = 'micropatterned'
    image_type = []
    gene_root_name = ""
    path_data = path.raw_data_dir

    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name, "r") as basic_file_handler, \
            h5py.File(path.data_dir+input_dir_name+'/'+secondary_file_name, "r") as secondary_file_handler:

        for molecule_type in molecule_types:
            peripheral_fraction_profile(secondary_file_handler, molecule_type, mrnas, 10,colors, image_type,gene_root_name,basic_file_handler)
            peripheral_fraction_profile(secondary_file_handler, molecule_type, mrnas, 30,colors, image_type,gene_root_name,basic_file_handler)
        histogram_peripheral_profile(basic_file_handler, secondary_file_handler, mrnas, proteins, colors, path_data,cell_type,image_type,periph_fraction_cst)

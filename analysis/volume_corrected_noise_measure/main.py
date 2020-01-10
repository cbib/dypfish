#!/usr/bin/python
# encoding: UTF-8

import sys
import h5py
import argparse
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.plot as plot
from src.utils import enable_logger, loadconfig, check_dir


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


def main():
    enable_logger()
    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"][0:4]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    volume_offset = configData["VOLUME_OFFSET"]
    size_coefficient = configData["SIZE_COEFFICIENT"]
    colors=configData["COLORS"]
    volume_coeff = ((1 / size_coefficient) ** 2) * 0.3


    ## Build spots density relative to cell area for arhgdia and arhgdia cultured
    ## Compare cell and nucleus area for arhgdia and arhgdia cultured
    with h5py.File(path.data_dir + input_dir_name + '/' + basic_file_name, "a") as file_handler, \
            h5py.File(path.data_dir + input_dir_name + '/' + secondary_file_name, "a") as sec_file_handler:

        ###First part
        # stan.compare_spots_density(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        # stan.compare_cell_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        # stan.compare_cell_volume(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        # stan.compare_nucleus_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        check_dir(path.analysis_dir + "/volume_corrected_noise_measure/figures/")


        ##Second part
        arhgdia = helps.build_image_list_2(file_handler, 'mrna', "arhgdia",["3h"])
        arhgdia_cultured = helps.build_image_list_2(file_handler, 'mrna', "arhgdia_cultured",["1h","3h"])
        #nm_arhgdia=stan.compute_volume_corrected_nm(file_handler, arhgdia,volume_offset, volume_coeff)
        nm_arhgdia=stan.compute_surface_corrected_nm(file_handler, arhgdia, size_coefficient)

        nm_arhgdia_cultured=stan.compute_surface_corrected_nm(file_handler, arhgdia_cultured, size_coefficient)
        plot.histogram_noise_measured(nm_arhgdia, nm_arhgdia_cultured)

        ##Build dynamic profiles
        molecule_type = ['/mrna']
        genes = ["beta_actin", "arhgdia", "gapdh", "pard3","pkp4","rab13"]
        for i in range(len(genes)):
            timepoints = ["2h", "3h", "4h", "5h"]
            nms=[]
            for timepoint in timepoints:
                #print(genes[i], '_', timepoint)
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, genes[i], [timepoint])
                print(image_list)
                nm = stan.compute_volume_corrected_nm(file_handler, image_list, volume_offset, volume_coeff)
                nms.append(nm)
            plot.noise_measured_dynamic_profile(nms,genes[i],colors[i])
            figname=path.analysis_dir + "volume_corrected_noise_measure/figures/nm_barplot_profile_" + genes[i]+".png"
            plot.bar_profile_simple(nms, figname, colors[i])

        nms = []
        for i in range(len(genes)):
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], genes[i])
            nm = stan.compute_volume_corrected_nm(file_handler, image_list, volume_offset, volume_coeff)
            nms.append(nm)

if __name__ == "__main__":
    main()

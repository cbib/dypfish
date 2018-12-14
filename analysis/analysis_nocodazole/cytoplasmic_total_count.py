#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.acquisition_descriptors as adsc
from src.path import basic_file_path,path_data, analysis_dir
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger

if __name__ == "__main__":

    # Required descriptors: spots, IF, cell mask an height_map
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    enable_logger()

    # Compute bar plot cytoplasmic total transcripts
    with h5py.File(basic_file_path, "r") as file_handler:

        molecule_type = ['/mrna']
        #colors = ['#1E95bb', '#1ec5d4']
        colors = ['#F16c1b', '#f1bc1b']
        #genes = ["arhgdia", "arhgdia_nocodazole"]
        genes = ["pard3", "pard3_nocodazole"]
        cytoplasmic_transcripts=[]
        for gene in genes:
            print(gene)
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path_data))
        figname= analysis_dir + 'analysis_nocodazole/figures/' + molecule_type[0] +'_total_cytoplasmic_transcript.png'
        stan.plot_bar_profile(cytoplasmic_transcripts, genes, 200, "Cytoplasmic total transcript", figname,colors)

        molecule_type = ['/protein']
        #proteins = ["arhgdia", "arhgdia_nocodazole", "arhgdia_cytod", "pard3", "pard3_nocodazole", "pard3_cytod"]
        #proteins = ["arhgdia", "arhgdia_nocodazole"]
        proteins = ["pard3", "pard3_nocodazole"]
        cytoplasmic_intensities = []
        for protein in proteins:
            print(protein)
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], protein)
            cytoplasmic_intensities.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path_data))
        print(cytoplasmic_intensities)
        figname = analysis_dir + 'analysis_nocodazole/figures/' + molecule_type[0] + '_total_cytoplasmic_intensity.png'
        stan.plot_bar_profile(cytoplasmic_intensities, genes, 100000000, "Cytoplasmic total intensity (arbitrary units)", figname,colors)
        figname = path.analysis_dir + 'analysis_nocodazole/figures/' + molecule_type[0] + '_total_cytoplasmic_intensity.png'
        stan.plot_bar_profile(cytoplasmic_intensities, proteins, 100000000, "Cytoplasmic total intensity (arbitrary units)", figname,colors)

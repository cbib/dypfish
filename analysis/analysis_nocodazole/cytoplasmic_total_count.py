#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.acquisition_descriptors as adsc
from src.path import basic_file_path,path_data, analysis_dir
import src.path as path
import src.helpers as helps
from src.utils import enable_logger
from src.utils import check_dir
import src.plot as plot

if __name__ == "__main__":
    # Required descriptors: spots, IF, cell mask an height_map
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    enable_logger()

    # Compute bar plot cytoplasmic total transcripts
    with h5py.File(basic_file_path, "r") as file_handler:
        molecule_type = ['/mrna']
        colors = ['#F16c1b', '#f1bc1b']
        genes = ["pard3", "pard3_nocodazole"]
        cytoplasmic_transcripts=[]
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path_data))
        figname= check_dir(analysis_dir + 'analysis_nocodazole/figures/cyt_total/') + molecule_type[0] +'_total_cytoplasmic_transcript.png'
        plot.bar_profile(cytoplasmic_transcripts, genes, figname)
        molecule_type = ['/protein']
        proteins = ["pard3", "pard3_nocodazole"]
        cytoplasmic_intensities = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], protein)
            cytoplasmic_intensities.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path_data))
        figname = check_dir(path.analysis_dir + 'analysis_nocodazole/figures/cyt_total') + molecule_type[0] + '_total_cytoplasmic_intensity.png'
        plot.bar_profile(cytoplasmic_intensities, proteins,figname)

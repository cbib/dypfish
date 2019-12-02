#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, check_dir, loadconfig, plot_colors_cytoD
import src.plot as plot


def main():
    enable_logger()

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts

    configData = loadconfig("cytoD")
    genes = configData["GENES"]


    with h5py.File(path.basic_file_path, "r") as file_handler:
        molecule_type = ['/mrna']
        cytoplasmic_transcripts = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            cytoplasmic_transcripts.append(adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data))
        plot_filename = molecule_type[0] + '_total_cytoplasmic_transcript_arhgdia_cytod.png'
        figname = check_dir(path.analysis_dir + 'analysis_cytoD/figures/') + plot_filename

        plot.bar_profile(cytoplasmic_transcripts, genes, figname)
if __name__ == "__main__":
    main()

#!/usr/bin/python
# encoding: UTF-8

import h5py
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.plot as plot
from src.utils import enable_logger

def main():
    enable_logger()

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8


    ## Build spots density relative to cell area for arhgdia and arhgdia cultured
    ## Compare cell and nucleus area for arhgdia and arhgdia cultured
    molecule_type = ['/mrna']
    with h5py.File(path.basic_file_path, "a") as file_handler, h5py.File(path.secondary_file_path, "a") as sec_file_handler:

        stan.compare_spots_density(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        #stan.compare_spots_volume_density(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        #stan.compare_spots_density_by_gene_and_timepoint(file_handler, "arhgdia", "arhgdia_cultured",['2h', '3h', '4h', '5h'], ['1h', '3h'])
        stan.compare_cell_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        stan.compare_cell_volume(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")
        stan.compare_nucleus_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_cultured")

        # stan.compare_spots_density(file_handler, sec_file_handler, "arhgdia", "arhgdia_scratch")
        # stan.compare_spots_volume_density(file_handler, sec_file_handler, "arhgdia", "arhgdia_scratch")
        # stan.compare_spots_density_by_gene_and_timepoint(file_handler, "arhgdia", "arhgdia_scratch",
        #                                                  ['2h', '3h', '4h', '5h'], ['1h', '3h', '5h'])
        # stan.compare_cell_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_scratch")
        # stan.compare_cell_volume(file_handler, sec_file_handler, "arhgdia", "arhgdia_scratch")
        # stan.compare_nucleus_area(file_handler, sec_file_handler, "arhgdia", "arhgdia_scratch")



        #compute and plot volume corrected noise measure
        arhgdia = helps.build_image_list_2(file_handler, 'mrna', "arhgdia",["3h"])
        arhgdia_cultured = helps.build_image_list_2(file_handler, 'mrna', "arhgdia_scratch",["3h"])
        nm_arhgdia=stan.compute_volume_corrected_nm(file_handler, arhgdia)
        nm_arhgdia_cultured=stan.compute_surface_corrected_nm(file_handler, arhgdia_cultured)

        plot.histogram_noise_measured(nm_arhgdia, nm_arhgdia_cultured)

        molecule_type = ['/mrna']
        genes = ["beta_actin", "arhgdia", "gapdh", "pard3","pkp4","rab13"]
        colors=['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']


        for i in range(len(genes)):
            timepoints = ["2h", "3h", "4h", "5h"]
            nms=[]
            for timepoint in timepoints:
                print(genes[i], '_', timepoint)
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, genes[i], [timepoint])
                nm = stan.compute_volume_corrected_nm(file_handler, image_list)
                nms.append(nm)

            plot.noise_measured_dynamic_profile(nms,genes[i],colors[i])


        nms = []
        for i in range(len(genes)):
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], genes[i])
            nm = stan.compute_volume_corrected_nm(file_handler, image_list)
            nms.append(nm)
        #print(nms)




if __name__ == "__main__":
    main()

#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues


import numpy as np
import h5py
from src import acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, check_dir, loadconfig
import src.plot as plot

def pearsoncorr(vec1, vec2):
    mu1 = np.mean(vec1)
    mu2 = np.mean(vec2)
    vec1b = vec1 - mu1
    vec2b = vec2 - mu2
    val = np.mean(vec1b * vec2b) / (np.std(vec1) * np.std(vec2))
    return val
'''this analysis aims to compare normalized profile of different descriptors peripheral fraction, cytoplasmic spread, cytoplasmic total count
and degree of clustering. We investigate correlations between pairs of such normaized distributions using the Pearson Correlation Coefficient'''

def main():
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()
    configData = loadconfig("original")
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    # micropatterned data

    with h5py.File(path.basic_file_path, "r") as file_handler,h5py.File(path.h_star_file_path, "r") as hstar_file_handler,h5py.File(path.mtoc_file_path, "r") as mtoc_file_handler,h5py.File(path.secondary_file_path, "r") as secondary_file_handler:
        molecule_type = ['/mrna']
        prof_m = np.zeros((4, 4, 4))
        cpt_g = 0
        for i in range(len(proteins)):
            cpt_tp=0
            for timepoint in timepoints_mrna:
                print(proteins[i], '_', timepoint)
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, proteins[i], [timepoint])

                cyt_spread = adsc.compute_cytoplasmic_spread(image_list,file_handler , path.path_data)
                prof_m[cpt_g,0,cpt_tp]=np.median(cyt_spread)

                periph_frac = adsc.compute_periph_fraction(image_list, file_handler, secondary_file_handler,  10,path.path_data)
                prof_m[cpt_g, 1, cpt_tp]=np.median(periph_frac)

                cyt_total=adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data)
                prof_m[cpt_g, 2, cpt_tp]=np.median(cyt_total)

                dof = adsc.compute_degree_of_clustering(image_list,hstar_file_handler, mtoc_file_handler)
                prof_m[cpt_g, 3, cpt_tp]=np.median(dof)

                cpt_tp+=1
            prof_m[cpt_g, 0, :] /= np.mean(prof_m[cpt_g, 0, :])
            prof_m[cpt_g, 1, :] /= np.mean(prof_m[cpt_g, 1, :])
            prof_m[cpt_g, 2, :] /= np.mean(prof_m[cpt_g, 2, :])
            prof_m[cpt_g, 3, :] /= np.mean(prof_m[cpt_g, 3, :])
            cpt_g+=1

        molecule_type = ['/protein']
        prof_p = np.zeros((4, 4, 4))
        cpt_g = 0
        for i in range(len(proteins)):
            cpt_tp = 0
            for timepoint in timepoints_protein:
                print(proteins[i], '_', timepoint)
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, proteins[i], [timepoint])

                cyt_spread = adsc.compute_cytoplasmic_spread(image_list, file_handler, path.path_data)
                prof_p[cpt_g, 0, cpt_tp] = np.median(cyt_spread)

                periph_frac = adsc.compute_periph_fraction(image_list,file_handler, secondary_file_handler,  10,path.path_data)
                prof_p[cpt_g, 1, cpt_tp] = np.median(periph_frac)

                cyt_total = adsc.compute_cytoplasmic_total(image_list, file_handler, path.path_data)
                prof_p[cpt_g, 2, cpt_tp] = np.median(cyt_total)

                dof = adsc.compute_degree_of_clustering(image_list, mtoc_file_handler, hstar_file_handler)
                prof_p[cpt_g, 3, cpt_tp] = np.median(dof)

                cpt_tp += 1
            prof_p[cpt_g, 0, :] /= np.mean(prof_p[cpt_g, 0, :])
            prof_p[cpt_g, 1, :] /= np.mean(prof_p[cpt_g, 1, :])
            prof_p[cpt_g, 2, :] /= np.mean(prof_p[cpt_g, 2, :])
            prof_p[cpt_g, 3, :] /= np.mean(prof_p[cpt_g, 3, :])
            cpt_g += 1

        vec_match = []
        vec_non_match = []
        for i in range(4):
            for j in range(4):
                vec1 = prof_m[i, :, :].flatten()
                vec2 = prof_p[j, :, :].flatten()
                v = pearsoncorr(vec1,vec2)
                if i == j:
                    vec_match.append(v)
                else:
                    vec_non_match.append(v)
        figpath = check_dir(path.analysis_dir + 'analysis_correlation_profile/figures/') + 'profile_correlation.png'
        plot.boxplot(vec_match,vec_non_match,figpath)

if __name__ == "__main__":
    main()

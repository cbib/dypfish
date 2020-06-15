#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import pprint as pp
import pandas as pd
from loguru import logger

import constants
import plot
import helpers
import numpy as np
from repository import H5RepositoryWithCheckpoint
from image_set import ImageSet
# this should be called as soon as possible
from path import global_root_dir
import itertools
import matplotlib.pyplot as plt
from scipy import stats



def pearsoncorr(vec1, vec2):
    mu1 = np.mean(vec1)
    mu2 = np.mean(vec2)
    vec1b = vec1 - mu1
    vec2b = vec2 - mu2
    val = np.mean(vec1b * vec2b) / (np.std(vec1) * np.std(vec2))
    print(val)
    return val

def calculate_temporal_interaction_score(mrna_data, protein_data, timepoint_num_mrna, timepoint_num_protein):
    S1 = helpers.get_forward_interactions(timepoint_num_mrna, timepoint_num_protein)
    interactions = np.zeros((len(timepoint_num_mrna), len(timepoint_num_protein)))
    for i in range(len(timepoint_num_mrna)):
        for j in range(len(timepoint_num_protein)):
            #interactions[i, j] = pearsoncorr(list(mrna_data[i]), list(protein_data[j]))
            interactions[i, j] = stats.pearsonr(list(mrna_data[i]), list(protein_data[j]))[0]
    (p, stat, ranking) = helpers.permutations_test(interactions, S1, size=len(timepoint_num_mrna))
    if len(timepoint_num_mrna)==4:
        #TODO if matrix 4 * 4
        tis = (100 - stat) / 64.0
    else:
        # TODO if matrix 2 * 2
        tis = (12 - stat) / 3.0

    return tis, p, ranking



def compute_protein_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num = 4):
    prot_tis_dict={}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['PROTEINS']:
        prot_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_PROTEIN']:
            image_set = ImageSet(analysis_repo, ["protein/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num, stripes = stripes)
            mrna_tp_df = pd.DataFrame(arr)
            prot_median.append(mrna_tp_df.mean(axis=0).values)
        prot_tis_dict[gene]=prot_median

    return prot_tis_dict


def compute_mrna_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num = 4):
    mrna_tis_dict={}
    stripes = constants.analysis_config['STRIPE_NUM']
    for gene in constants.analysis_config['MRNA_GENES']:
        mrna_median = []
        for timepoint in constants.dataset_config['TIMEPOINTS_MRNA']:
            image_set = ImageSet(analysis_repo, ["mrna/{0}/{1}/".format(gene, timepoint)])
            arr = image_set.compute_normalized_quadrant_and_slice_densities(quadrants_num=quadrants_num, stripes = stripes)
            print(len(arr[0]))
            mrna_tp_df = pd.DataFrame(arr)
            mrna_median.append(mrna_tp_df.mean(axis=0).values)
        mrna_tis_dict[gene]=mrna_median

    return mrna_tis_dict


def compute_heatmap(ranking, gene, size=4, xtickslabel=['2h', '3h', '5h', '7h'], ytickslabel = ['2h', '3h', '4h', '5h']):
    im = np.flipud(np.kron(ranking, np.ones((10, 10))))
    plt.imshow(im, extent=[0, size, 0, size], cmap='GnBu', interpolation='nearest')
    ax = plt.axes()
    ax.set_ylabel("mRNA  - Time (hrs)")
    ax.set_xlabel("Protein  - Time (hrs)")
    myxticklabels = xtickslabel
    ax.xaxis.set(ticks=np.arange(0.5, size + 0.5, 1), ticklabels=myxticklabels)
    myyticklabels = ytickslabel
    ax.yaxis.set(ticks=np.arange(0.5, size + 0.5, 1), ticklabels=myyticklabels)
    ax.set_title(gene)
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS'].format(gene=gene)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    plt.savefig(tgt_fp)
    plt.close()


# # # ################################## ################################## ################################## #################################
# # # ################################## ################################## ################################## #################################
# # # # #                               Figure 5C et 5D Analysis TIS for original data
# # # ################################## ################################## ################################## #################################
# # # ################################## ################################## ################################## #################################


logger.info("Temporal interaction score for the mRNA original data")
constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/temporal_interactions/config_original.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)


target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), "original_mrna_dataframe.csv")

#mrna_df=pd.DataFrame(columns=['Gene'])
mrna_tis_dict = compute_mrna_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
#mrna_df = pd.concat([mrna_df, pd.DataFrame(mrna_tis_dict)])
#mrna_df.to_csv(target_df_fp)
#mrna_df = pd.read_csv(target_df_fp, index_col=None)

target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), "original_protein_dataframe.csv")
#prot_df=pd.DataFrame(columns=['Gene'])
prot_tis_dict = compute_protein_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
#prot_df = pd.concat([prot_df, pd.DataFrame(prot_tis_dict)])
#prot_df.to_csv(target_df_fp)
#prot_df = pd.read_csv(target_df_fp, index_col=None)


tiss=[]
p_vals=[]
for gene in constants.analysis_config['PROTEINS']:
    mrna_list= mrna_tis_dict[gene]
    prot_list = prot_tis_dict[gene]
    print("gene:", gene)
    #mrna_list = mrna_df['Gene'].values
    #prot_list = prot_df['Gene'].values
    (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list, constants.dataset_config['TIMEPOINTS_NUM_MRNA'], constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
    tiss.append(tis)
    p_vals.append(p)
    compute_heatmap(ranking, gene)

tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS_HISTOGRAM']
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),tgt_image_name)
plot.bar_profile(tiss, constants.analysis_config['PROTEINS'], tgt_fp)



# # # ################################## ################################## ################################## #################################
# # # ################################## ################################## ################################## #################################
# # # # #                               Figure 6E Analysis TIS for nocodazole arhgdia et pard3 data
# # # ################################## ################################## ################################## #################################
# # # ################################## ################################## ################################## #################################


# logger.info("Temporal interaction score for the mRNA nocodazole data")
# constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/temporal_interactions/config_nocodazole_arhgdia.json"))
# dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
# primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
# secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
# analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
#
# logger.info("Analysis Temporal interaction score for nocodazole arhgdia FISH data")
# target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
#                             "arhdgia_nocodazole_mrna_dataframe.csv")
# #mrna_df=pd.DataFrame(columns=['Gene'])
# mrna_tis_dict = compute_mrna_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
# #mrna_df = pd.concat([mrna_df, pd.DataFrame(mrna_tis_dict)])
# #mrna_df.to_csv(target_df_fp)
#
# #mrna_df = pd.read_csv(target_df_fp, index_col=None)
#
# logger.info("Analysis MTOC enrichment for the nocodazole arhgdia IF data")
# target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
#                             "arhdgia_nocodazole_protein_dataframe.csv")
#
# #prot_df=pd.DataFrame(columns=['Gene'])
# prot_tis_dict = compute_protein_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
# #prot_df = pd.concat([prot_df, pd.DataFrame(prot_tis_dict)])
# #prot_df.to_csv(target_df_fp)
#
# #prot_df = pd.read_csv(target_df_fp, index_col=None)
#
# tiss=[]
# p_vals=[]
# for gene in constants.analysis_config['PROTEINS']:
#     mrna_list= mrna_tis_dict[gene]
#     prot_list = prot_tis_dict[gene]
#     #mrna_list = mrna_df['Gene'].values
#     #prot_list = prot_df['Gene'].values
#     (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list, constants.dataset_config['TIMEPOINTS_NUM_MRNA'], constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
#     tiss.append(tis)
#     p_vals.append(p)
#     compute_heatmap(ranking, gene)
#
# tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS_HISTOGRAM']
# tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),tgt_image_name)
# plot.bar_profile(tiss, constants.analysis_config['PROTEINS'], tgt_fp)





# logger.info("Temporal interaction score for the mRNA pard3 nocodazole data")
# constants.init_config(analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/temporal_interactions/config_nocodazole_pard3.json"))
# dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
# primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
# secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
# analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)
#
# logger.info("Analysis Temporal interaction score for nocodazole pard3 FISH data")
# target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
#                             "pard3_nocodazole_mrna_dataframe.csv")
# #mrna_df=pd.DataFrame(columns=['Gene'])
# mrna_tis_dict = compute_mrna_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
# #mrna_df = pd.concat([mrna_df, pd.DataFrame(mrna_tis_dict)])
# #mrna_df.to_csv(target_df_fp)
#
# #mrna_df = pd.read_csv(target_df_fp, index_col=None)
#
# logger.info("Analysis MTOC enrichment for the nocodazole pard3 IF data")
# target_df_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
#                             "pard3_nocodazole_protein_dataframe.csv")
#
# #prot_df=pd.DataFrame(columns=['Gene'])
# prot_tis_dict = compute_protein_relative_density_per_quadrants_and_slices(analysis_repo, quadrants_num=8)
# #prot_df = pd.concat([prot_df, pd.DataFrame(prot_tis_dict)])
# #prot_df.to_csv(target_df_fp)
#
# #prot_df = pd.read_csv(target_df_fp, index_col=None)
#
# tiss=[]
# p_vals=[]
# for gene in constants.analysis_config['PROTEINS']:
#     mrna_list= mrna_tis_dict[gene]
#     prot_list = prot_tis_dict[gene]
#     #mrna_list = mrna_df['Gene'].values
#     #prot_list = prot_df['Gene'].values
#     (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list, constants.dataset_config['TIMEPOINTS_NUM_MRNA'], constants.dataset_config['TIMEPOINTS_NUM_PROTEIN'])
#     tiss.append(tis)
#     p_vals.append(p)
#     compute_heatmap(ranking, gene, size=2, xtickslabel=['3h','5h'], ytickslabel=['3h','5h'])
#
# tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_TIS_HISTOGRAM']
# tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),tgt_image_name)
# plot.bar_profile(tiss, constants.analysis_config['PROTEINS'], tgt_fp)


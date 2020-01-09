#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import src.path as path

from src.utils import enable_logger, plot_colors,check_dir, loadconfig


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

def permutations(orig_list):
    if not isinstance(orig_list, list):
        orig_list = list(orig_list)
    yield orig_list
    if len(orig_list) == 1:
        return
    for n in sorted(orig_list):
        new_list = orig_list[:]
        pos = new_list.index(n)
        del (new_list[pos])
        new_list.insert(0, n)
        for resto in permutations(new_list[1:]):
            if new_list[:1] + resto <> orig_list:
                yield new_list[:1] + resto

def using_indexed_assignment(x):
    "https://stackoverflow.com/a/5284703/190597 (Sven Marnach)"
    result = np.empty(len(x), dtype=int)
    temp = x.argsort()
    result[temp] = np.arange(len(x))
    return result

def permutations_test(interactions, fwdints):
    import itertools
    fwdints = fwdints.astype(bool)
    vals = interactions.flatten()
    indx = using_indexed_assignment(vals)
    one_matrix = np.ones((2, 2)).astype(int)
    indx_matrix = np.matrix(indx.reshape((2, 2)))
    indx_matrix = np.add(indx_matrix, one_matrix)
    ranking = indx_matrix.copy()
    rs0 = np.sum(indx_matrix[fwdints[:]])
    rs1 = np.sum(indx_matrix[fwdints[:] == 0])
    perms = [x for x in itertools.permutations([0, 1], 4)]
    nps = len(perms)
    rs = []
    for p1 in range(nps):
        for p2 in range(nps):
            test = indx_matrix.copy()
            for i in range(4):
                np.random.shuffle(test[:, i])
            rs.append(np.sum(test[fwdints[:]]))
    count = 0
    for score in rs:
        if score > rs0:
            count += 1
    if count==0:
        p=0
    stat = rs1
    return p, stat, ranking

def pearsoncorr(vec1, vec2):
    mu1 = np.mean(vec1)
    mu2 = np.mean(vec2)
    vec1b = vec1 - mu1
    vec2b = vec2 - mu2
    val = np.mean(vec1b * vec2b) / (np.std(vec1) * np.std(vec2))
    return val

def get_forward_interactions(mrna_timepoints, protein_timepoints):
    X = len(mrna_timepoints)
    Y = len(protein_timepoints)
    fwd_interactions = np.zeros((X, Y))
    for x in range(X):
        for y in range(Y):
            if protein_timepoints[y] > mrna_timepoints[x]:
                fwd_interactions[x, y] = 1
    return fwd_interactions

def calculate_temporal_interaction_score(mrna_data, protein_data):
    S1 = get_forward_interactions([3, 5],[3, 5])
    interactions = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            interactions[i, j] = stats.pearsonr(list(mrna_data[i]), list(protein_data[j]))[0]
    (p, stat,ranking) = permutations_test(interactions, S1)
    tis = (100 - stat) / 64.0;
    return (tis, p, ranking)

def plot_bar_profile(data, genes, figname, colors):
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    y_lim = np.max(data) + 0.3
    ind = np.arange(N)
    width = 0.35
    rects1 = ax.bar(ind, data, width,
                    color=colors)
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim + 0.1)
    ax.set_title('')
    ax.set_xticks(ind)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    plt.yticks(fontsize=25)
    plt.savefig(figname, format='svg')
    plt.close()
def main():


    # WARNING you need to run before MTOC analysis script called search enriched quad

    enable_logger()
    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"]
    proteins = configData["PROTEINS"]
    mrna_timepoints = configData["TIMEPOINTS_MRNA"]
    prot_timepoints = configData["TIMEPOINTS_PROTEIN"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    mtoc_file_name = configData["MTOC_FILE_NAME"]
    colors = configData["COLORS"]

    #mrnas = ["arhgdia_nocodazole", "pard3_nocodazole"]
    #mrna_timepoints = ["3h", "5h"]
    #proteins = ["arhgdia_nocodazole", "pard3_nocodazole"]
    #prot_timepoints = ["3h", "5h"]

    tiss = []
    p_vals = []
    count_gene = 0
    for mrna in mrnas:
        mrna_list = []
        prot_list = []
        for timepoint in mrna_timepoints:
            mrna_df = pd.read_csv(path.analysis_dir+"nocodazole/dataframe/periph_"+mrna + '_' + timepoint  + "_mrna.csv", index_col=0)
            mrna_list.append(mrna_df.median(axis=0).values)
        for timepoint in prot_timepoints:
            prot_df = pd.read_csv(path.analysis_dir+"nocodazole/dataframe/periph_"+mrna + '_' + timepoint  +"_protein.csv", index_col=0)
            prot_list.append(prot_df.median(axis=0).values)
        (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list)
        tiss.append(tis)
        p_vals.append(p)
        im = np.flipud(np.kron(ranking, np.ones((10, 10))))
        plt.imshow(im, extent=[0, 2, 0, 2], cmap='GnBu', interpolation='nearest')
        ax = plt.axes()
        ax.set_ylabel("mRNA  - Time (hrs)")
        ax.set_xlabel("Protein  - Time (hrs)")
        myxticklabels = ['3h', '5h']
        ax.xaxis.set(ticks=np.arange(0.5, 2.5, 1), ticklabels=myxticklabels)
        myyticklabels = ['3h', '5h']
        ax.yaxis.set(ticks=np.arange(0.5, 2.5, 1), ticklabels=myyticklabels)
        ax.set_title(mrna)
        figname = check_dir(path.analysis_dir + 'nocodazole/figures/TIS/') + 'periph_' + mrna + '_TIS_correlation_ranking.svg'
        plt.savefig(figname, format='svg')
        plt.close()
        count_gene += 1

    ylabel = 'Global temporal interaction score'
    figname = check_dir(path.analysis_dir + 'nocodazole/figures/TIS/')+'periph_' + 'TIS.svg'
    plot_bar_profile(tiss, mrnas,figname, plot_colors)

if __name__ == "__main__":
    main()

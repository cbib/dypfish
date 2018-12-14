#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import math

import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
from skimage.draw import circle
from skimage import measure
import pandas as pd

import os
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
import src.constants

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
np.set_printoptions(precision=4)
# if "log" not in globals():
# logger = Logger.init_logger('REFORMAT_%s'%(cfg.language_code), load_config())
logger.info("Running %s", sys.argv[0])


def permutations (orig_list):
    if not isinstance(orig_list, list):
        orig_list = list(orig_list)

    yield orig_list

    if len(orig_list) == 1:
        return

    for n in sorted(orig_list):
        new_list = orig_list[:]
        pos = new_list.index(n)
        del(new_list[pos])
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

def permutations_test(interactions,fwdints):
    fwdints = fwdints.astype(bool)
    ranks = np.zeros(interactions.shape)
    vals=interactions.flatten()
    print(vals)
    indx=using_indexed_assignment(vals)
    print(indx)
    import itertools
    one_matrix=np.ones((2,2)).astype(int)

    indx_matrix=np.matrix(indx.reshape((2,2)))
    indx_matrix=np.add(indx_matrix,one_matrix)
    print(indx_matrix)
    ranking=indx_matrix.copy()

    rs0=np.sum(indx_matrix[fwdints[:]])
    rs1 = np.sum(indx_matrix[fwdints[:]==0])

    perms = [x for x in itertools.permutations([0, 1], 4)]
    print(perms)
    nps = len(perms)
    print(nps)
    rs = []

    for p1 in range(nps):
        for p2 in range(nps):
            test=indx_matrix.copy()
            for i in range(4):
                np.random.shuffle(test[:, i])
            rs.append(np.sum(test[fwdints[:]]))
    print(rs)
    count=0
    for score in rs:
        if score> rs0:
            count+=1
    p=0
    stat=rs1
    print(stat)
    print(ranking)

    return (p, stat,ranking)


def pearsoncorr(vec1,vec2):
    mu1 = np.mean(vec1);
    mu2 = np.mean(vec2);

    vec1b = vec1 - mu1
    vec2b = vec2 - mu2

    val = np.mean(vec1b * vec2b) / (np.std(vec1) * np.std(vec2));
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
    #print(S1)
    interactions = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            #interactions[i, j] = pearsoncorr(list(mrna_data[i]), list(protein_data[j]))
            interactions[i, j] = stats.pearsonr(list(mrna_data[i]), list(protein_data[j]))[0]
    (p, stat,ranking) = permutations_test(interactions, S1)
    tis = (100 - stat) / 64.0;
    #print(tis)



    return (tis, p, ranking)


def plot_bar_profile(data,genes,ylabel,figname,colors):
    ## third technic
    fig = plt.figure()
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ## the data
    N = len(genes)

    y_lim=np.max(data)+0.3

    # for l in data:
    #     max_v=np.median(l)
    #     if max_v > y_lim:
    #         y_lim=max_v
    #     dataMedians.append(np.median(l))
    #     dataStd.append(np.std(l))
    ## necessary variables
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars
    #colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']


    ## the bars
    rects1 = ax.bar(ind, data, width,
                    color=colors)

    # axes and labels
    ax.set_xlim(-width, len(ind) + width)

    ax.set_ylim(0, y_lim+0.1)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    xTickMarks = ["" for i in range(0, N)]
    ax.set_xticks(ind)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.legend([gene for gene in genes], loc='upper right')
    ax.legend(rects1, genes,prop={'size':8})
    plt.savefig(figname,format='svg')
    #plt.show()
    plt.close()


if __name__ == "__main__":

    # Required descriptors: spots, IF, cell mask an height_map

    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    ''' 
    1-You need to create a password.txt file before running to connect via ssh
    '''
    basic_file_path = path.analysis_data_dir + 'basic.h5'
    secondary_file_path = path.analysis_data_dir + 'secondary.h5'
    mtoc_file_path = path.analysis_data_dir + 'mtoc.h5'
    path_data = path.raw_data_dir


    mrnas = [ "arhgdia_nocodazole","pard3_nocodazole"]
    mrna_timepoints = [ "3h", "5h"]

    proteins = ["arhgdia_nocodazole","pard3_nocodazole"]
    prot_timepoints = [ "3h", "5h"]
    colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B']

    tiss=[]
    p_vals=[]
    #ranks=np.zeros((4,4,4))
    count_gene=0
    for mrna in mrnas:
        mrna_list = []
        prot_list = []
        for timepoint in mrna_timepoints:
            print(mrna, timepoint)
            mrna_df = pd.read_csv(path.analysis_dir+"analysis_nocodazole/df/"+mrna + '_' + timepoint + "_mrna.csv", index_col=0)
            #mrna_df=pd.DataFrame(np.matrix(mrna_df)[:,::4])
            mrna_list.append(mrna_df.median(axis=0).values)

        for timepoint in prot_timepoints:
            prot_df = pd.read_csv(path.analysis_dir+"analysis_nocodazole/df/"+mrna + '_' + timepoint + "_protein.csv", index_col=0)
            #prot_df = pd.DataFrame(np.matrix(prot_df)[:, ::4])
            prot_list.append(prot_df.median(axis=0).values)
        (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list)
        tiss.append(tis)
        p_vals.append(p)
        #ranks[:,:,count_gene]=ranking
        print(mrna)
        print(ranking)

        # print(np.kron(ranking, np.ones((10, 10))))
        im = np.flipud(np.kron(ranking, np.ones((10, 10))))

        plt.imshow(im, extent=[0, 2, 0, 2],cmap='GnBu',interpolation='nearest')
        #plt.imshow(ranking,cmap='GnBu')

        ax = plt.axes()
        ax.set_ylabel("mRNA  - Time (hrs)")
        ax.set_xlabel("Protein  - Time (hrs)")

        myxticklabels=['3h', '5h']
        #ax.set_xticklabels(['2h', '3h', '4h', '5h'], minor=True)
        ax.xaxis.set(ticks=np.arange(0.5,2.5,1), ticklabels=myxticklabels)
        #ax.set_yticklabels(['2h', '3h', '5h', '7h'], minor=True)
        myyticklabels=['3h', '5h']

        ax.yaxis.set(ticks=np.arange(0.5,2.5,1), ticklabels=myyticklabels)
        ax.set_title(mrna)
        #plt.show()
        figname = path.analysis_dir + 'analysis_nocodazole/figures/'+mrna+'_TIS_correlation_ranking.svg'
        plt.savefig(figname,format='svg')
        plt.close()
        count_gene+=1


    # print(p_vals)
    # print(tiss)

    ylabel='Global temporal interaction score'
    figname = path.analysis_dir + 'analysis_nocodazole/figures/TIS.svg'
    plot_bar_profile(tiss,mrnas,ylabel,figname,colors)


    #print(ranks)



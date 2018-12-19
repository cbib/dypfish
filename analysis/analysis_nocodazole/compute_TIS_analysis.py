#!/usr/bin/python
# encoding: UTF-8
# author: benjamin Dartigues


import logging
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import src.path as path

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
np.set_printoptions(precision=4)
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
    perms = [x for x in itertools.permutations([0, 1], 2)]
    print(perms)
    nps = len(perms)
    print(nps)
    rs = []

    for p1 in range(nps):
        for p2 in range(nps):
            test=indx_matrix.copy()
            for i in range(2):
                np.random.shuffle(test[:, i])
            rs.append(np.sum(test[fwdints[:]]))
    count=0
    for score in rs:
        if score> rs0:
            count+=1
    p = float(count / float(len(rs)))
    print(p)
    p=0
    stat=rs1
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
    interactions = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            interactions[i, j] = stats.pearsonr(list(mrna_data[i]), list(protein_data[j]))[0]
    (p, stat,ranking) = permutations_test(interactions, S1)
    tis = (100 - stat) / 64.0;
    print(p)
    return (tis, p, ranking)

def plot_bar_profile(data,genes,ylabel,figname,colors):
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    y_lim=np.max(data)+0.3
    ind = np.arange(N)
    width = 0.35
    rects1 = ax.bar(ind, data, width,
                    color=colors)
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim+0.1)
    ax.set_ylabel(ylabel)
    ax.set_title('')
    ax.set_xticks(ind)
    plt.legend([gene for gene in genes], loc='upper right')
    ax.legend(rects1, genes,prop={'size':8})
    plt.savefig(figname,format='svg')
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
    mrnas = ["arhgdia", "pard3"]
    mrna_timepoints = [ "3h", "5h"]
    proteins = ["arhgdia","pard3"]
    prot_timepoints = [ "3h", "5h"]
    colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B']
    tiss=[]
    p_vals=[]
    count_gene=0
    for mrna in mrnas:
        mrna_list = []
        prot_list = []
        for timepoint in mrna_timepoints:
            print(mrna, timepoint)
            mrna_df = pd.read_csv(path.analysis_dir+"analysis_nocodazole/df/"+mrna + '_' + timepoint + "_mrna.csv", index_col=0)
            mrna_list.append(mrna_df.median(axis=0).values)
        for timepoint in prot_timepoints:
            prot_df = pd.read_csv(path.analysis_dir+"analysis_nocodazole/df/"+mrna + '_' + timepoint + "_protein.csv", index_col=0)
            prot_list.append(prot_df.median(axis=0).values)
        (tis, p, ranking) = calculate_temporal_interaction_score(mrna_list, prot_list)
        tiss.append(tis)
        p_vals.append(p)
        print(mrna)
        print(ranking)
        print (p_vals)
        im = np.flipud(np.kron(ranking, np.ones((10, 10))))
        plt.imshow(im, extent=[0, 2, 0, 2],cmap='GnBu',interpolation='nearest')
        ax = plt.axes()
        ax.set_ylabel("mRNA  - Time (hrs)")
        ax.set_xlabel("Protein  - Time (hrs)")
        myxticklabels=['3h', '5h']
        ax.xaxis.set(ticks=np.arange(0.5,2.5,1), ticklabels=myxticklabels)
        myyticklabels=['3h', '5h']
        ax.yaxis.set(ticks=np.arange(0.5,2.5,1), ticklabels=myyticklabels)
        ax.set_title(mrna)
        figname = path.analysis_dir + 'analysis_nocodazole/figures/'+mrna+'_TIS_correlation_ranking.svg'
        plt.savefig(figname,format='svg')
        plt.close()
        count_gene+=1

    ylabel='Global temporal interaction score'
    figname = path.analysis_dir + 'analysis_nocodazole/figures/TIS.svg'
    plot_bar_profile(tiss,mrnas,ylabel,figname,colors)
#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import math
import pandas as pd
pd.set_option('display.max_rows', 500)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import interpolate
from scipy import stats
from skimage.draw import circle
from skimage import measure

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


def box_plot_concentrations(k1, k2, title, figname):
    #k1 = sum(k1[:, 0], [])
    #print(k1)
    #k2 = sum(k2[:, 0], [])
    #k3 = sum(k3[:, 0], [])
    #k1_k2med = k1 - np.median(k2)
    fig, ax = plt.subplots(1, 1)
    plt.boxplot([np.matrix(k1), np.matrix(k2)], 'gD')
    plt.ylabel('Mean mrna concentration (total mrna spots) by quadrant')
    plt.title(title)
    y_lim = np.max(k1)
    plt.axis([0, 2.5, 0, y_lim + y_lim / 10])
    plt.axis([0, 2.5, -1, y_lim])

    MTOC_patch = mpatches.Patch(label='1 - Quadrants with MTOC')
    noMTOC_patch = mpatches.Patch(label='2 - Other quadrants')
    plt.legend(handles=[MTOC_patch, noMTOC_patch])
    ax.yaxis.grid(color='gray', linestyle='dashed')
    plt.savefig(figname)
    #plt.show()
    plt.close()

# function for setting the colors of the box plots pairs
def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')


def plot_bar_profile(data,random_data,data_noc, random_data_noc,genes,noc_genes, y_limit, ylabel, figname, colors):
    ## third technic
    fig = plt.figure()
    # ax = fig.add_subplot(111)
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    ## the data
    N = len(genes)

    dataMedians = []
    dataStd = []
    dataMediansnoc = []
    dataStdnoc = []
    y_lim_max = np.max(data) + 0.4
    y_lim_min = np.min(data) - 0.4
    print(y_lim_max)

    for l in random_data:
        dataMedians.append(np.median(l))
        dataStd.append(np.std(l))
    for l in random_data_noc:
        dataMediansnoc.append(np.median(l))
        dataStdnoc.append(np.std(l))


    ## necessary variables
    ind = np.arange(N)  # the x locations for the groups
    width = 0.20  # the width of the bars
    # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
    # colors =['#0A3950','#1E95BB','#A1BA6D','#F16C1B','#C02A18','#E9CB45']
    #colors = ['#1E95bb', '#1ec5d4']
    colors = ['#F16c1b', '#f1bc1b']
    ## the bars
    rects1 = ax.bar(ind, data, width,
                    color=['#F16c1b','#1ec5d4'],yerr=dataStd,
                    error_kw=dict(elinewidth=1, ecolor='black'))
    rects2 = ax.bar(ind+width, data_noc, width,
                    color=['#f1bc1b','#1E95BB'],yerr=dataStdnoc,
                    error_kw=dict(elinewidth=1, ecolor='black'))

    # axes and labels\
    print(len(ind))
    ax.set_xlim(-width, len(ind)- 0.4)

    ax.set_ylim(y_lim_min, y_lim_max)
    #ax.set_ylabel(ylabel)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ax.set_title('')

    ax.set_xticks(ind + width)
    #ax.set_xticklabels(('arhgdia', 'pard3'))
    ax.set_xticklabels(('pard3','arhgdia'))

    #ax.set_xticklabels('arhgdia')
    #ax.legend((rects1[0], rects2[0]), ('control', 'Nocodazole'))
    #plt.legend([gene for gene in genes], loc='upper right')
    #ax.legend(rects1, genes, prop={'size': 8})
    #ax.legend(rects2, noc_genes, prop={'size': 8})
    plt.yticks(fontsize=25)
    plt.savefig(figname, format='png')
    plt.show()

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
    basic_md5_path=path.analysis_data_dir + 'basic.md5'
    secondary_md5_path=path.analysis_data_dir + 'secondary.md5'
    mtoc_md5_path = path.analysis_data_dir + 'mtoc.md5'
    # pwd_path = path.analysis_data_dir + 'password.txt'
    path_data = path.raw_data_dir
    # f = open(pwd_path, 'r')
    # loginpwd = f.readline()
    # login = loginpwd.split(":")[0]
    # pwd = loginpwd.split(":")[1]
    # helps.check_data(basic_file_path, 'basic', login, pwd)
    # helps.check_data(mtoc_file_path, 'mtoc', login, pwd)
    with h5py.File(basic_file_path, "r") as file_handler,h5py.File(secondary_file_path, "r") as second_file_handler,h5py.File(mtoc_file_path, "r") as mtoc_file_handler:
        from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes
        molecule_type=['/mrna']
        #mrnas = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
        colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']

        # Paired MTOC quadrants with non MTOC quadrants
        #mrnas_noc = [ "arhgdia_nocodazole"]
        mrnas_noc = ["pard3_nocodazole"]

        #df = pd.read_csv('df/old_df/global_mtoc_file_all_mrna_nocodazole')
        df = pd.read_csv('df/global_mtoc_file_mrna_all.csv')

        df_sorted_noc=df.sort('MTOC', ascending=False).groupby(['Gene','timepoint','Image'], as_index=False).first()
        mpis_noc=[]
        mpis_random_noc=[]
        for gene, line in df_sorted_noc.groupby(['Gene']):
            if gene in mrnas_noc:
                gene_random_mpis = []
                for i in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
                    gene_random_mpis.append(mpi)
                mpis_random_noc.append(gene_random_mpis)
                mpi,p=stan.calculate_mpi(line['MTOC'].values,line['Non MTOC'].values)
                mpis_noc.append(mpi)

        #mrnas = ["arhgdia"]
        mrnas = ["pard3"]

        df = pd.read_csv('df/global_mtoc_file_mrna_all.csv')
        df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        mpis = []
        random_mpis=[]
        for gene, line in df_sorted.groupby(['Gene']):
            if gene in mrnas:
                gene_random_mpis = []
                for i in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values / 100000, line['Non MTOC'].values / 100000)
                    gene_random_mpis.append(mpi)
                random_mpis.append(gene_random_mpis)
                mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'].values)
                mpis.append(mpi)
        figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[0] + '_paired_mpis_nocodazole.png'
        plot_bar_profile(mpis,random_mpis, mpis_noc,mpis_random_noc, mrnas_noc, mrnas, 0.7, 'Cytoplasmic MPI', figname,colors)

###########################################################################3


        # #MTOC Quadrant with random selcted non-MTOC quadrant
        # mrnas = ["arhgdia", "beta_actin", "gapdh", "pard3", "pkp4", "rab13"]
        # mrnas = ["arhgdia_nocodazole", "pard3_nocodazole"]
        # df = pd.read_csv('df/global_mtoc_file_all_mrna_nocodazole')
        # df_sorted=df.sort('MTOC', ascending=False).groupby(['Gene','timepoint','Image'], as_index=False).first()
        # random_mpis_noc=[]
        # for gene, line in df_sorted.groupby(['Gene']):
        #     gene_random_mpis = []
        #     for i in range(100):
        #         mpi,p=stan.calculate_random_mpi(line['MTOC'].values,line['Non MTOC'].values)
        #         gene_random_mpis.append(mpi)
        #     random_mpis_noc.append(np.median(np.array(gene_random_mpis)))
        #
        # mrnas = ["arhgdia", "pard3"]
        # df = pd.read_csv('df/global_mtoc_file_mrna_all')
        # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # random_mpis = []
        # for gene, line in df_sorted.groupby(['Gene']):
        #     if gene == "arhgdia" or gene == "pard3":
        #         gene_random_mpis = []
        #         for i in range(100):
        #             mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
        #             gene_random_mpis.append(mpi)
        #         random_mpis.append(np.median(np.array(gene_random_mpis)))
        # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[0] + '__bootstrap_mpis_mpis_nocodazole.svg'
        # plot_bar_profile(random_mpis, random_mpis_noc, mrnas_noc, mrnas, 0.7, 'Cytoplasmic MPI', figname, colors)

#######################################################################################

        df = pd.read_csv('df/global_mtoc_file_mrna_all.csv')
        array = ['arhgdia', 'arhgdia_nocodazole']
        df=df.loc[df['Gene'].isin(array)]

        df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        import seaborn as sns
        df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'], var_name='Quadrants')
        my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}

        box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants',palette=my_pal)
        box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
        box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
        #box.set(ylabel='log2 volumic density')
        box.set(xlabel='')
        box.set(ylabel='')
        box.set_xlabel("", fontsize=15)
        box.set_ylabel("", fontsize=15)
        plt.yticks(fontsize=25)

        box.legend_.remove()

        figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'
        plt.savefig(figname, format='svg')
        plt.show()

        # df = pd.read_csv('df/global_mtoc_file_all_mrna_nocodazole')
        # df_sorted=df.sort('MTOC', ascending=False).groupby(['Gene','timepoint','Image'], as_index=False).first()
        # df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        # df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        # df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],var_name='Quadrants')
        #
        # box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants')
        # #box.set(ylim=(0, 5))
        # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[0] + '_boxplot_MTOC_enrichment_nocodazole.svg'
        # plt.savefig(figname, format='svg')
        # plt.show()
        #
        #
        #
        #
        #
        #
        #
        #
        #
        #
        # molecule_type = ['/protein']
        # # mrnas = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
        # colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']
        # # Paired MTOC quadrants with non MTOC quadrants
        #
        #
        # proteins_noc = ["arhgdia_nocodazole", "pard3_nocodazole"]
        # df = pd.read_csv('df/global_mtoc_file_all_protein_nocodazole')
        # df_sorted_noc = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # mpis_noc = []
        # for gene, line in df_sorted_noc.groupby(['Gene']):
        #     mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'].values)
        #     mpis_noc.append(mpi)
        #
        # proteins = ["arhgdia", "pard3"]
        # df = pd.read_csv('df/global_mtoc_file_protein_all')
        # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # mpis = []
        # for gene, line in df_sorted.groupby(['Gene']):
        #     if gene == "arhgdia" or gene == "pard3":
        #         mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'].values)
        #         mpis.append(mpi)
        #
        # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[0] + '_paired_mpis_nocodazole.svg'
        # plot_bar_profile(mpis, mpis_noc, proteins_noc, proteins, 0.7, 'Cytoplasmic MPI', figname, colors)



####################################################################################33















        # # MTOC Quadrant with random selcted non-MTOC quadrant
        # mrnas = ["arhgdia_nocodazole", "pard3_nocodazole"]
        # df = pd.read_csv('df/global_mtoc_file_all_protein_nocodazole')
        # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # random_mpis_noc = []
        # for gene, line in df_sorted.groupby(['Gene']):
        #     gene_random_mpis = []
        #     for i in range(100):
        #         mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
        #         gene_random_mpis.append(mpi)
        #     random_mpis_noc.append(np.median(np.array(gene_random_mpis)))
        #
        # mrnas = ["arhgdia", "pard3"]
        # df = pd.read_csv('df/global_mtoc_file_protein_all')
        # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # random_mpis = []
        # for gene, line in df_sorted.groupby(['Gene']):
        #     if gene == "arhgdia" or gene == "pard3":
        #         gene_random_mpis = []
        #         for i in range(100):
        #             mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
        #             gene_random_mpis.append(mpi)
        #         np.median(np.array(gene_random_mpis))
        #         random_mpis.append(np.median(np.array(gene_random_mpis)))
        #
        # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[
        #     0] + '__bootstrap_mpis_nocodazole.svg'
        # plot_bar_profile(random_mpis, random_mpis_noc, mrnas_noc, mrnas, 0.7, 'Cytoplasmic MPI', figname, colors)
        #
        # # df = pd.read_csv('df/global_mtoc_file_protein_all')
        # # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # # import seaborn as sns
        # #
        # # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],
        # #              var_name='Quadrants')
        # # box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants')
        # # box.set(ylim=(0, 3.5))
        # # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[
        # #     0] + '_boxplot_MTOC_enrichment.svg'
        # # plt.savefig(figname, format='svg')
        # # plt.show()
        # #
        # # df = pd.read_csv('df/global_mtoc_file_all_protein_nocodazole')
        # # df_sorted = df.sort('MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
        # # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],
        # #              var_name='Quadrants')
        # # box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants')
        # # box.set(ylim=(0, 5))
        # # figname = path.analysis_dir + 'analysis_nocodazole/figures/MPI/' + molecule_type[
        # #     0] + '_boxplot_MTOC_enrichment_nocodazole.svg'
        # # plt.savefig(figname, format='svg')
        # # plt.show()








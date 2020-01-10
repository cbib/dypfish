#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import argparse
import numpy as np
import h5py
import pandas as pd
pd.set_option('display.max_rows', 500)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import src.path as path
import src.statistical_analysis as stan
import src.plot as plot
import seaborn as sns
from src.utils import loadconfig,check_dir
from pylab import setp

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
np.set_printoptions(precision=4)
logger.info("Running %s", sys.argv[0])

parser = argparse.ArgumentParser()
parser.add_argument("--peripheral", "-p", help='boolean flag: perform peripheral computation or not',
                    action="store_true", default=False)
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name
is_periph = args.peripheral

def box_plot_concentrations(k1, k2, title, figname):
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
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
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
    ind = np.arange(N)
    width = 0.20
    ax.bar(ind, [dataMedians, dataMediansnoc], width, color=colors)
    print(len(ind))
    ax.set_xlim(-width, len(ind)- 0.4)
    ax.set_ylim(y_lim_min, y_lim_max)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ax.set_title('')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('pard3','arhgdia'))
    plt.yticks(fontsize=25)
    plt.savefig(figname, format='png')

if __name__ == "__main__":

    configData = loadconfig(input_dir_name)
    mrnas = configData["GENES"]
    proteins = configData["PROTEINS"]
    mrna_timepoints = configData["TIMEPOINTS_MRNA"]
    prot_timepoints = configData["TIMEPOINTS_PROTEIN"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    mtoc_file_name = configData["MTOC_FILE_NAME"]
    colors = configData["COLORS"]
    bootstrap_mpi = configData["BOOTSTRAP_MPI"]

    path_data = path.raw_data_dir
    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name, "r") as file_handler,\
            h5py.File(path.data_dir+input_dir_name+'/'+secondary_file_name, "r") as second_file_handler,\
            h5py.File(path.data_dir+input_dir_name+'/'+mtoc_file_name, "r") as mtoc_file_handler:
        molecule_type=['mrna']
        colors = ['#F16c1b', '#f1bc1b']
        mrnas = ['pard3']
        mrnas_noc = ['pard3_nocodazole']

        check_dir(path.analysis_dir + 'nocodazole/figures/')

        if is_periph:
            mrna_df = pd.read_csv(path.analysis_dir + 'nocodazole/dataframe/periph_global_mtoc_file_mrna.csv')
            protein_df = pd.read_csv(path.analysis_dir + 'nocodazole/dataframe/periph_global_mtoc_file_protein.csv')
        else:
            mrna_df = pd.read_csv(path.analysis_dir + 'nocodazole/dataframe/global_mtoc_file_mrna.csv')
            protein_df = pd.read_csv(path.analysis_dir + 'nocodazole/dataframe/global_mtoc_file_protein.csv')



        ####Compute MPI
        #df = pd.read_csv('df/global_mtoc_file_mrna.csv')
        df_sorted_noc = mrna_df.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()
        mpis_noc=[]
        mpis_random_noc=[]
        for gene, line in df_sorted_noc.groupby(['Gene']):
            if gene in mrnas_noc:
                gene_random_mpis = []
                nonMTOC = line['Non MTOC1'].values + line['Non MTOC2'].values + line['Non MTOC3'].values

                for i in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values, nonMTOC, bootstrap_mpi)
                    gene_random_mpis.append(mpi)
                mpis_random_noc.append(gene_random_mpis)
                mpi,p=stan.calculate_mpi(line['MTOC'].values,nonMTOC)
                mpis_noc.append(mpi)


        df_sorted = mrna_df.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()

        mpis = []
        random_mpis=[]
        for gene, line in df_sorted.groupby(['Gene']):
            if gene in mrnas:
                gene_random_mpis = []
                nonMTOC = line['Non MTOC1'].values + line['Non MTOC2'].values + line['Non MTOC3'].values

                for i in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values , nonMTOC, bootstrap_mpi)
                    gene_random_mpis.append(mpi)
                random_mpis.append(gene_random_mpis)
                mpi, p = stan.calculate_mpi(line['MTOC'].values, nonMTOC)
                mpis.append(mpi)

        if is_periph:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'periph_' + molecule_type[0] + '_mpis_nocodazole.png'
        else:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + molecule_type[0] + '_mpis_nocodazole.png'

        #plot_bar_profile(mpis,random_mpis, mpis_noc,mpis_random_noc, mrnas_noc, mrnas, 0.7, 'Cytoplasmic MPI', figname,colors)

        ####Compute MTOC quadrant enrichment fro mRNA

        #df = pd.read_csv('df/global_mtoc_file_mrna.csv')
        array = ['arhgdia', 'arhgdia_nocodazole']
        mrna_df_c=mrna_df.loc[mrna_df['Gene'].isin(array)]
        df_sorted = mrna_df_c.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()

        # df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        # df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        # df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        #dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'], var_name='Quadrants')

        dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC1', 'Non MTOC2', 'Non MTOC3', 'MTOC leading edge'],var_name='Quadrants')
        dd = dd.replace('Non MTOC1', 'Non MTOC')
        dd = dd.replace('Non MTOC2', 'Non MTOC')
        dd = dd.replace('Non MTOC3', 'Non MTOC')
        dd = dd.replace(0.000000, np.nan)
        dd['value'] = dd['value'].apply(np.log2)

        my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}
        box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants',palette=my_pal)
        box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
        box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
        box.set(xlabel='')
        box.set(ylabel='')
        box.set_xlabel("", fontsize=15)
        box.set_ylabel("", fontsize=15)
        plt.yticks(fontsize=25)
        box.legend_.remove()
        figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') +'arhgdia_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'

        if is_periph:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'periph_arhgdia_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'
        else:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'arhgdia_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'

        plt.savefig(figname, format='svg')
        plt.close()




        array = ['pard3', 'pard3_nocodazole']
        mrna_df = mrna_df.loc[mrna_df['Gene'].isin(array)]
        df_sorted = mrna_df.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()

        # df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        # df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        # df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],var_name='Quadrants')

        dd = pd.melt(df_sorted, id_vars=['Gene'],value_vars=['MTOC', 'Non MTOC1', 'Non MTOC2', 'Non MTOC3', 'MTOC leading edge'],var_name='Quadrants')
        dd = dd.replace('Non MTOC1', 'Non MTOC')
        dd = dd.replace('Non MTOC2', 'Non MTOC')
        dd = dd.replace('Non MTOC3', 'Non MTOC')
        dd = dd.replace(0.000000, np.nan)
        dd['value'] = dd['value'].apply(np.log2)


        my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}
        box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants', palette=my_pal)
        box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
        box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
        box.set(xlabel='')
        box.set(ylabel='')
        box.set_xlabel("", fontsize=15)
        box.set_ylabel("", fontsize=15)
        plt.yticks(fontsize=25)
        box.legend_.remove()

        if is_periph:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'periph_pard3_' + molecule_type[
                0] + '_boxplot_MTOC_enrichment.svg'
        else:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'pard3_' + molecule_type[
                0] + '_boxplot_MTOC_enrichment.svg'

        plt.savefig(figname, format='svg')
        plt.close()



        ####Compute MTOC quadrant enrichment for proteins
        molecule_type = ['protein']
        array = ['arhgdia', 'arhgdia_nocodazole']
        protein_df_c = protein_df.loc[protein_df['Gene'].isin(array)]
        df_sorted = protein_df_c.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()

        # df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        # df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        # df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'], var_name='Quadrants')

        dd = pd.melt(df_sorted, id_vars=['Gene'],
                     value_vars=['MTOC', 'Non MTOC1', 'Non MTOC2', 'Non MTOC3', 'MTOC leading edge'],
                     var_name='Quadrants')
        dd = dd.replace('Non MTOC1', 'Non MTOC')
        dd = dd.replace('Non MTOC2', 'Non MTOC')
        dd = dd.replace('Non MTOC3', 'Non MTOC')
        dd = dd.replace(0.000000, np.nan)
        dd['value'] = dd['value'].apply(np.log2)

        my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}
        box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants', palette=my_pal)
        box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
        box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
        box.set(xlabel='')
        box.set(ylabel='')
        box.set_xlabel("", fontsize=15)
        box.set_ylabel("", fontsize=15)
        plt.yticks(fontsize=25)
        box.legend_.remove()
        if is_periph:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'periph_arhgdia_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'
        else:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'arhgdia_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'


        plt.savefig(figname, format='svg')
        plt.close()
        array = ['pard3', 'pard3_nocodazole']
        protein_df = protein_df.loc[protein_df['Gene'].isin(array)]
        df_sorted = protein_df.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()

        # df_sorted['MTOC'] = df_sorted['MTOC'].apply(np.log2)
        # df_sorted['Non MTOC'] = df_sorted['Non MTOC'].apply(np.log2)
        # df_sorted['MTOC leading edge'] = df_sorted['MTOC leading edge'].apply(np.log2)
        # dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],var_name='Quadrants')

        dd = pd.melt(df_sorted, id_vars=['Gene'],
                     value_vars=['MTOC', 'Non MTOC1', 'Non MTOC2', 'Non MTOC3', 'MTOC leading edge'],
                     var_name='Quadrants')
        dd = dd.replace('Non MTOC1', 'Non MTOC')
        dd = dd.replace('Non MTOC2', 'Non MTOC')
        dd = dd.replace('Non MTOC3', 'Non MTOC')
        dd = dd.replace(0.000000, np.nan)
        dd['value'] = dd['value'].apply(np.log2)

        my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}
        box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants', palette=my_pal)
        box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
        box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
        box.set(xlabel='')
        box.set(ylabel='')
        box.set_xlabel("", fontsize=15)
        box.set_ylabel("", fontsize=15)
        plt.yticks(fontsize=25)
        box.legend_.remove()

        if is_periph:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'periph_pard3_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'
        else:
            figname = check_dir(path.analysis_dir + 'nocodazole/figures/MTOC/') + 'pard3_' + molecule_type[0] + '_boxplot_MTOC_enrichment.svg'

        plt.savefig(figname, format='svg')
        plt.close()


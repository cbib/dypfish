#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

'''Plot Figures XXX and YYY for the manuscript bla bla'''

import numpy as np
import pandas as pd
import argparse

from src.plot import *
import src.path as path
import src.statistical_analysis as stan
from src.utils import check_dir, enable_logger, loadconfig, my_pal

pd.set_option('display.max_rows', 500)

parser = argparse.ArgumentParser()
parser.add_argument("--peripheral", "-p", help = 'boolean flag: perform peripheral computation or not',
                        action = "store_true", default = False)
args = parser.parse_args()
is_periph = args.peripheral

def compute_mpis(df_sorted):
    mpis = []
    random_mpis = []
    for gene, line in df_sorted.groupby(['Gene']):
        gene_random_mpis = []
        for i in range(100):
            mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'])
            gene_random_mpis.append(mpi)
        random_mpis.append(gene_random_mpis)
        mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'])
        mpis.append(mpi)
    err = []
    for l in random_mpis:
        err.append(np.std(l))
    return mpis,err

def main(is_periph=False):

    configData = loadconfig("original")
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    mime_type = configData["PNG_IMAGES_MIME_TYPE"]


    check_dir(path.analysis_dir + 'analysis_MTOC/figures/')
    enable_logger()
    if is_periph:
        print ('Periph mode')
    if is_periph:
        df_filename_m=check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_mrna.csv'
        df_filename_p=check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'periph_global_mtoc_file_all_protein.csv'
    else:
        print ('Global mode')
        df_filename_m=check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_mrna.csv'
        df_filename_p=check_dir(path.analysis_dir + 'analysis_MTOC/dataframe/') + 'global_mtoc_file_all_protein.csv'

    #mrna part
    df_m = pd.read_csv(df_filename_m)
    df_sorted=df_m.sort_values(by='MTOC',ascending=False).groupby(['Gene','timepoint','Image'], as_index=False).first()

    # plot bar plot cytoplasmic mpi
    mpis,err=compute_mpis(df_sorted)

    if is_periph:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_mrna_paired_mpis.png'
    else:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'mrna_paired_mpis.png'
    bar_profile_median(mpis, genes, 'mrna', figname, err)

    # plot boxplot cytoplasmic mpi
    dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'], var_name='Quadrants')

    if is_periph:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_mrna_boxplot_MTOC_enrichment.png'
    else:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'mrna_boxplot_MTOC_enrichment.png'
    # log values for plotting
    dd['value'] = dd['value'].apply(np.log2)
    sns_boxplot(dd,my_pal,figname)

    #protein part
    df_p = pd.read_csv(df_filename_p)
    df_sorted = df_p.sort_values(by='MTOC',ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()

    # plot bar plot cytoplasmic mpi
    mpis,err=compute_mpis(df_sorted)
    if is_periph:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_protein_paired_mpis.png'
    else:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'protein_paired_mpis.png'
    bar_profile_median(mpis,proteins,'protein',figname,err)

    # plot boxplot cytoplasmic mpi
    dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'],var_name='Quadrants')
    if is_periph:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_protein_boxplot_MTOC_enrichment.png'
    else:
        figname=check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'protein_boxplot_MTOC_enrichment.png'
    # log values for plotting
    dd['value'] = dd['value'].apply(np.log2)
    sns_boxplot(dd,my_pal,figname)

    # plot dynamic profile cytoplasmic mpi
    df_sorted_m = df_m.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()
    df_sorted_p = df_p.sort_values(by='MTOC', ascending=False).groupby(['Gene', 'timepoint', 'Image'],as_index=False).first()
    for i in range(len(genes)):
        df_sorted_mc = df_sorted_m.copy()
        df_sorted_mc = df_sorted_mc[df_sorted_mc.Gene == genes[i]]
        random_mpis = []
        mpis = []
        for gene, line in df_sorted_mc.groupby(['timepoint']):
            gene_random_mpis = []
            for j in range(100):
                mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
                gene_random_mpis.append(mpi)
            random_mpis.append(gene_random_mpis)
            mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'].values)
            print(gene, p)
            mpis.append(mpi)
        mrna_data = np.zeros((3, 4))
        counter = 0
        for timepoint in timepoints_mrna:
            print(genes[i], '_', timepoint)
            err = np.median(np.abs(np.tile(np.median(random_mpis[counter]), (1, len(random_mpis[counter]))) - random_mpis[counter]))
            upp_env = mpis[counter] + err
            low_env = mpis[counter] - err
            mrna_data[0, counter] = mpis[counter]
            mrna_data[1, counter] = upp_env
            mrna_data[2, counter] = low_env
            counter += 1
        mrna_data = mrna_data / np.mean(mrna_data[0, :])
        if genes[i] in proteins:
            df_sorted_pc = df_sorted_p.copy()
            df_sorted_pc = df_sorted_pc[df_sorted_pc.Gene == genes[i]]
            random_mpis = []
            mpis = []
            for gene, line in df_sorted_pc.groupby(['timepoint']):
                gene_random_mpis = []
                for j in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values, line['Non MTOC'].values)
                    gene_random_mpis.append(mpi)
                random_mpis.append(gene_random_mpis)
                mpi, p = stan.calculate_mpi(line['MTOC'].values, line['Non MTOC'].values)
                mpis.append(mpi)
            counter = 0
            protein_data = np.zeros((3, 4))
            for timepoint in timepoints_protein:
                err = np.median(np.abs(np.tile(np.median(random_mpis[counter]), (1, len(random_mpis[counter]))) - random_mpis[counter]))
                upp_env = mpis[counter] + err
                low_env = mpis[counter] - err
                protein_data[0, counter] = mpis[counter]
                protein_data[1, counter] = upp_env
                protein_data[2, counter] = low_env
                counter += 1
            protein_data = protein_data / np.mean(protein_data[0, :])

        figname = check_dir(path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_' if is_periph else '' + 'MPI_' + genes[i] + '.png'
        if is_periph:
            figname = check_dir(
                path.analysis_dir + 'analysis_MTOC/figures/') + 'periph_MPI_' + genes[i] + '.png'
        else:
            figname = check_dir(path.analysis_dir + 'analysis_MTOC/figures/')+ 'MPI_' + genes[i] + '.png'
        dynamic_profiles(mrna_data, protein_data, timepoints_num_mrna, timepoints_num_protein, genes[i], plot_colors[i], '', '', figname)

if __name__ == "__main__":
    # main(is_periph=True)
    main(is_periph=False)

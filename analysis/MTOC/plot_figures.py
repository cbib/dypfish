#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

'''Plot Figures XXX and YYY for the manuscript bla bla'''

import pandas as pd
import argparse
import sys
from src.plot import *
import src.path as path
import src.statistical_analysis as stan
from src.utils import check_dir, check_dir_only, enable_logger, loadconfig, my_pal

pd.set_option('display.max_rows', 500)


parser = argparse.ArgumentParser()                                               
parser.add_argument("--peripheral", "-p", help = 'boolean flag: perform peripheral computation or not',
                        action = "store_true", default = False)
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)

args = parser.parse_args()
is_periph = args.peripheral
input_dir_name = args.input_dir_name


def compute_mpis(df_sorted,bootstrap_mpi):
    #print(df_sorted)
    mpis = []
    random_mpis = []
    for gene, line in df_sorted.groupby(['Gene']):
        # print(gene)
        # print("#################################")
        # print("non mTOC 1", line['Non MTOC1'].values)
        # print("#################################")
        # print("non mTOC 2", line['Non MTOC2'].values)
        # print("#################################")
        # print("non mTOC 3", line['Non MTOC3'].values)
        # print("#################################")


        gene_random_mpis = []
        #nonMTOC =line['Non MTOC1'].values + line['Non MTOC2'].values+ line['Non MTOC3'].values
        #nonMTOC = np.array(line['Non MTOC1'].values) + np.array(line['Non MTOC2'].values) + np.array(line['Non MTOC3'].values)
        nonMTOC=np.concatenate((np.array(line['Non MTOC1'].values), np.array(line['Non MTOC2'].values),np.array(line['Non MTOC3'].values)), axis=0)

        #print("NON MTOC concatenated ",nonMTOC)
        #sys.exit()
        for i in range(100):
            mpi, p = stan.calculate_random_mpi(line['MTOC'].values,nonMTOC,bootstrap_mpi)
            gene_random_mpis.append(mpi)
        random_mpis.append(gene_random_mpis)
        mpi, p = stan.calculate_mpi(line['MTOC'].values, nonMTOC)
        mpis.append(mpi)
    err = []
    for l in random_mpis:
        err.append(np.std(l))
    return mpis,err

def main():

    configData = loadconfig(input_dir_name)
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    bootstrap_mpi = configData["BOOTSTRAP_MPI"]
    mime_type = configData["PNG_IMAGES_MIME_TYPE"]
    colors=configData["COLORS"]




    check_dir(path.analysis_dir + 'MTOC/figures/')
    if (check_dir_only(path.analysis_dir + 'MTOC/figures/')==False):
        sys.exit(1)

    enable_logger()
    if is_periph:
        print ('Periph mode')
        df_filename_m=check_dir(path.analysis_dir + 'MTOC/dataframe/') + 'periph_global_mtoc_file_mrna.csv'
        df_filename_p=check_dir(path.analysis_dir + 'MTOC/dataframe/') + 'periph_global_mtoc_file_protein.csv'
    else:
        print ('Global mode')
        df_filename_m=check_dir(path.analysis_dir + 'MTOC/dataframe/') + 'global_mtoc_file_mrna.csv'
        df_filename_p=check_dir(path.analysis_dir + 'MTOC/dataframe/') + 'global_mtoc_file_protein.csv'

    #mrna part
    df_m = pd.read_csv(df_filename_m)
    df_sorted=df_m.sort_values(by='MTOC',ascending=False).groupby(['Gene','timepoint','Image'], as_index=False).first()

    # plot bar plot cytoplasmic mpi
    mpis,err=compute_mpis(df_sorted,bootstrap_mpi)

    if is_periph:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'periph_mrna_paired_mpis.png'
    else:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'mrna_paired_mpis.png'

    bar_profile_median(mpis, genes, figname, err, colors)

    # plot boxplot cytoplasmic mpi
    #dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC', 'MTOC leading edge'], var_name='Quadrants')

    #print(len(df_sorted['MTOC'].values))
    dd = pd.melt(df_sorted, id_vars=['Gene'],value_vars=['MTOC','Non MTOC1','Non MTOC2','Non MTOC3','MTOC leading edge'],var_name='Quadrants')

    dd=dd.replace('Non MTOC1', 'Non MTOC')
    dd=dd.replace('Non MTOC2', 'Non MTOC')
    dd=dd.replace('Non MTOC3', 'Non MTOC')

    dd=dd.replace(0.000000, np.nan)
    #print(len(dd['value'].values))

    #dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC','MTOC leading edge'], var_name='Quadrants' )
    print(dd.tail(n=200))
    #sys.exit()

    if is_periph:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'periph_mrna_boxplot_MTOC_enrichment.png'
    else:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'mrna_boxplot_MTOC_enrichment.png'
    # log values for plotting
    dd['value'] = dd['value'].apply(np.log2)


    sns_boxplot(dd,my_pal,figname)

    #protein part
    df_p = pd.read_csv(df_filename_p)
    df_sorted = df_p.sort_values(by='MTOC',ascending=False).groupby(['Gene', 'timepoint', 'Image'], as_index=False).first()
    #df_sorted = df_sorted[df_sorted['timepoint']!="5h"]
    #df_sorted = df_sorted[df_sorted['timepoint'] != "7h"]
    # plot bar plot cytoplasmic mpi
    mpis,err=compute_mpis(df_sorted,bootstrap_mpi)
    if is_periph:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'periph_protein_paired_mpis.png'
    else:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'protein_paired_mpis.png'
    bar_profile_median(mpis,proteins,figname,err, colors)

    # plot boxplot cytoplasmic mpi
    #dd = pd.melt(df_sorted, id_vars=['Gene'], value_vars=['MTOC', 'Non MTOC1', 'MTOC leading edge'], var_name='Quadrants')
    dd = pd.melt(df_sorted, id_vars=['Gene'],value_vars=['MTOC', 'Non MTOC1', 'Non MTOC2', 'Non MTOC3', 'MTOC leading edge'], var_name='Quadrants')

    dd = dd.replace('Non MTOC1', 'Non MTOC')
    dd = dd.replace('Non MTOC2', 'Non MTOC')
    dd = dd.replace('Non MTOC3', 'Non MTOC')
    if is_periph:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'periph_protein_boxplot_MTOC_enrichment.png'
    else:
        figname=check_dir(path.analysis_dir + 'MTOC/figures/') + 'protein_boxplot_MTOC_enrichment.png'
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
            #nonMTOC = line['Non MTOC1'].values + line['Non MTOC2'].values + line['Non MTOC3'].values
            nonMTOC = np.concatenate((np.array(line['Non MTOC1'].values), np.array(line['Non MTOC2'].values),np.array(line['Non MTOC3'].values)), axis=0)
            #nonMTOC=np.mean((np.array(line['Non MTOC1'].values), np.array(line['Non MTOC2'].values),np.array(line['Non MTOC3'].values)))
            for j in range(100):

                mpi, p = stan.calculate_random_mpi(line['MTOC'].values, nonMTOC, bootstrap_mpi)
                gene_random_mpis.append(mpi)
            random_mpis.append(gene_random_mpis)
            mpi, p = stan.calculate_mpi(line['MTOC'].values, nonMTOC)
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
        #mrna_data = mrna_data / np.mean(mrna_data[0, :])

        if genes[i] in proteins:
            df_sorted_pc = df_sorted_p.copy()
            df_sorted_pc = df_sorted_pc[df_sorted_pc.Gene == genes[i]]
            random_mpis = []
            mpis = []
            for gene, line in df_sorted_pc.groupby(['timepoint']):
                gene_random_mpis = []
                #nonMTOC = line['Non MTOC1'].values + line['Non MTOC2'].values + line['Non MTOC3'].values
                nonMTOC = np.concatenate((np.array(line['Non MTOC1'].values), np.array(line['Non MTOC2'].values),np.array(line['Non MTOC3'].values)), axis=0)
                #nonMTOC = np.mean((np.array(line['Non MTOC1'].values), np.array(line['Non MTOC2'].values),np.array(line['Non MTOC3'].values)))

                for j in range(100):
                    mpi, p = stan.calculate_random_mpi(line['MTOC'].values, nonMTOC, bootstrap_mpi)
                    gene_random_mpis.append(mpi)
                random_mpis.append(gene_random_mpis)
                mpi, p = stan.calculate_mpi(line['MTOC'].values, nonMTOC)
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
            #protein_data = protein_data / np.mean(protein_data[0, :])

        figname = check_dir(path.analysis_dir + 'MTOC/figures/') + 'periph_' if is_periph else '' + 'MPI_' + genes[i] + '.png'
        if is_periph:
            figname = check_dir(
                path.analysis_dir + 'MTOC/figures/') + 'periph_MPI_' + genes[i] + '.png'
        else:
            figname = check_dir(path.analysis_dir + 'MTOC/figures/')+ 'MPI_' + genes[i] + '.png'
        dynamic_profiles(mrna_data, protein_data, timepoints_num_mrna, timepoints_num_protein, genes[i], plot_colors[i], '', '', figname)

if __name__ == "__main__":
    # main(is_periph=True)
    main()

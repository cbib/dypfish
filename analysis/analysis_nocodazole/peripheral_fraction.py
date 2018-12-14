#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
import src.plot as plot
import os
import src.constants as constants
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import check_dir

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
# if "log" not in globals():
# logger = Logger.init_logger('REFORMAT_%s'%(cfg.language_code), load_config())
logger.info("Running %s", sys.argv[0])
import hashlib

def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


def peripheral_profile_by_timepoint(file_handler,molecule_type,genes,timepoints,colors,cell_type,image_type=[]):
    if len(image_type)==0:

        for timepoint in timepoints:
            mean_profiles = []
            for gene in genes:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, gene, [timepoint])
                profiles = adsc.compute_peripheral_fraction_profiles(file_handler, image_list)
                mean_profiles.append(np.average(profiles, axis=0))

                # boxplot_name = path.figure_dir + 'boxplot_2h' + gene + '.svg'
                # stan.boxplot_profile(profiles,gene,100,boxplot_name)
            plot_name = path.analysis_dir + 'analysis_nocodazole/figures/'+cell_type+'/peripheral_fraction_' + timepoint + '_timepoint' + str(
                constants.NUM_CONTOURS) + 'contours.svg'

            mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
            figure_title = 'timepoint: ' + timepoint
            plot.profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title,colors, True)
    else:
        for image_t in image_type:
            for timepoint in timepoints:
                mean_profiles = []
                for gene in genes:

                        image_list = helps.preprocess_image_list4(file_handler, molecule_type, gene, [timepoint],image_t)

                        profiles = adsc.compute_peripheral_fraction_profiles(file_handler, image_list)
                        mean_profiles.append(np.average(profiles, axis=0))

                        # boxplot_name = path.figure_dir + 'boxplot_2h' + gene + '.svg'
                        # stan.boxplot_profile(profiles,gene,100,boxplot_name)
                plot_name = path.analysis_dir + 'analysis_nocodazole/figures/' +cell_type + '/peripheral_fraction_' + timepoint + '_timepoint' + '_'+image_t+'_'+ str(constants.NUM_CONTOURS) + 'contours.svg'

                mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
                figure_title = 'timepoint: ' + timepoint
                plot.profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)


def peripheral_profile_by_gene(file_handler,molecule_type,genes,colors,cell_type,gene_root_name,image_type=[]):
    if len(image_type) == 0:
        mean_profiles = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            profiles = adsc.compute_peripheral_fraction_profiles_2D(file_handler, image_list)
            #print(np.average(profiles, axis=0))
            mean_profiles.append(np.average(profiles, axis=0))

            # boxplot_name = path.figure_dir + 'boxplot_2h' + gene + '.svg'
            # stan.boxplot_profile(profiles,gene,100,boxplot_name)
        plot_name = path.analysis_dir + 'analysis_nocodazole/figures/'+cell_type+'/'+gene_root_name+'_peripheral_fraction_all_timepoint' + str(constants.NUM_CONTOURS) +'contours.png'
        mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[0], len(genes), 1)
        figure_title = ''
        print(mean_profiles)

        # df = pd.DataFrame({"arhgdia": mean_profiles[0], 'arhgdia_nocodazole': mean_profiles[1], 'pard3': mean_profiles[2],
        #                    'pard3_nocodazole': mean_profiles[3],"gapdh":mean_profiles[4]},
        #                   )
        # #df=pd.DataFrame(mean_profiles, columns=["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"])
        # df.to_csv("peripheral_fraction_profile.csv")
        plot.profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)
    else:
        mean_profiles_by_conditions=[]
        for image_t in image_type:
            mean_profiles = []
            for gene in genes:
                image_list = helps.preprocess_image_list5(file_handler, molecule_type[0], gene,image_t)
                profiles = adsc.compute_peripheral_fraction_profiles(file_handler, image_list)
                # print(np.average(profiles, axis=0))
                mean_profiles.append(np.average(profiles, axis=0))

                # boxplot_name = path.figure_dir + 'boxplot_2h' + gene + '.svg'
                # stan.boxplot_profile(profiles,gene,100,boxplot_name)
            mean_profiles_by_conditions.append(mean_profiles)
            plot_name = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'contours.svg'
            #mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
            figure_title = 'timepoint: All'
            print(mean_profiles)

            # df = pd.DataFrame({"beta_actin": mean_profiles[0], 'arhgdia': mean_profiles[1], 'gapdh': mean_profiles[2],
            #                    'pard3': mean_profiles[3], 'pkp4': mean_profiles[4], "rab13": mean_profiles[5]},
            #                   )
            # # df=pd.DataFrame(mean_profiles, columns=["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"])
            # df.to_csv("peripheral_fraction_profile_"+image_t+".csv")
            plot.profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)

        norm_LPA=[]
        norm_plus_medium=[]
        for i in range(len(mean_profiles_by_conditions[0])):
            norm_LPA.append(mean_profiles_by_conditions[1][i] / mean_profiles_by_conditions[0][i])
            norm_plus_medium.append(mean_profiles_by_conditions[2][i] / mean_profiles_by_conditions[0][i])
        plot_name = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'LPAcontours.svg'
        figure_title = 'timepoint: All'
        plot.profile(norm_LPA, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)

        plot_name = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'plusmediumcontours.svg'
        figure_title = 'timepoint: All'
        plot.profile(norm_plus_medium, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)


def peripheral_fraction_profile(file_handler,molecule_type,genes,fraction,colors,image_type,gene_root_name,basic_file_handler):
    if len(image_type) == 0:
        fractions = []
        for gene in genes:
            print(gene)
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            fractions.append(adsc.build_histogram_periph_fraction(file_handler, image_list, fraction,path_data,basic_file_handler))
        figname = path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/'+molecule_type[0]+'_peripheral_fraction_'+str(fraction)+'.png'
        plot.fraction_profile(fractions, fraction, genes, figname,colors)
    else:
        for image_t in image_type:
            fractions = []
            for gene in genes:
                image_list = helps.preprocess_image_list5(file_handler, molecule_type[0], gene,image_t)
                fractions.append(adsc.build_histogram_periph_fraction(file_handler, image_list, fraction,path_data,basic_file_handler))
            figname = path.analysis_dir + 'analysis_nocodazole/figures/peripheral_fraction/' +gene_root_name+'_peripheral_fraction_' +image_t+'_'+ str(fraction) + '.png'
            plot.fraction_profile(fractions, fraction, genes, figname, colors)


def peripheral_fractions_profile(file_handler, molecule_type, genes,image_type):

    for i in range(1,11):
        # fig = plt.figures(figsize=(10, 10))
        # plt.ylabel('mRNA enrichment')
        # plt.xlabel('Images')
        # plt.title("peripheral fraction " + str(i))

        c = ["r", "b"]
        gene_counter = 0

        arrs=[]
        for gene in genes:
            # image_list = preprocess_image_list3(file_handler, molecule_type, gene, ['4h'])
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            gene_counter += 1


            arr = adsc.compute_fraction_profile(file_handler, image_list, i)
            arrs.append(arr)

        # fig = plt.figures(figsize=(10, 10))
        #
        # plt.plot(np.arange(len(arrs[0])), arrs[0], color=c[0], linewidth=3, label=genes[0])
        # plt.plot(np.arange(len(arrs[1])), arrs[1], color=c[1], linewidth=3, label=genes[1])
        #
        # # fig = plt.figures(figsize=(10, 10))
        # # plt.plot(arr, color=c[gene_counter], linewidth=3, label=gene)
        # plt.legend(loc='upper right')
        # plt.title("peripheral fraction " + str(i))
        # plt.show()

        #fig = plt.figures(figsize=(10, 10))
        #plt.hist(arr, bins=20, color=c[gene_counter], facecolor='red', alpha=0.75, range=(0, 1), label=gene)

        fig = idsc.plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 2, 1)
        ax1 = fig.add_subplot(1, 2, 2)
        idsc.plt.title("peripheral fraction " + str(i))
        ax.hist(arrs[0], color='r', bins=50, facecolor='red', alpha=0.75, range=(0, 1),
                label=genes[0])
        ax1.hist(arrs[1], bins=50, color='g', facecolor='green', alpha=0.75, range=(0, 1),
                 label=genes[1])
        ax.legend(loc='upper right')
        ax1.legend(loc='upper right')
        #plt.xlabel('Mrna enrichment')
        #plt.ylabel('Image count')
        idsc.plt.show()
        figname=path.figure_dir + "periph_fraction_hist/arhgdia_gapdh_peripheral_fraction_"+str(i)
        fig.savefig(figname)
        idsc.plt.close()


def histogram_peripheral_profile(basic_file_handler,secondary_file_handler, molecule_type, genes, proteins, colors, path_data,cell_type,image_type):
    if len(image_type) == 0:
        periph_fraction = []
        for gene in genes:
            print(gene)

            image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], gene)
            periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
        print(periph_fraction)
        figname = path.analysis_dir + 'analysis_nocodazole/figures/'+cell_type+'/'+molecule_type[
            0] +'_peripheral_fraction.png'
        plot.bar_profile(periph_fraction, genes, figname)


        molecule_type = ['/protein']
        periph_fraction = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], protein)
            periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD,path_data))
        figname = path.analysis_dir + 'analysis_nocodazole/figures/'+cell_type+'/'+molecule_type[0] +'_peripheral_fraction.png'
        plot.bar_profile(periph_fraction, proteins,  figname)
    else:
        for image_t in image_type:
            print(image_t)
            periph_fraction = []
            molecule_type = ['/mrna']
            for gene in genes:
                print(gene)

                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], gene,image_t)
                periph_fraction.append(
                    adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,
                                                 constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))

            figname = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/' + \
                      molecule_type[0] +'_'+image_t+ 'peripheral_fraction.svg'
            plot.bar_profile(periph_fraction, genes, figname)

            molecule_type = ['/protein']
            periph_fraction = []
            for protein in proteins:
                print(protein)
                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], protein,image_t)
                periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
            figname = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/' + \
                      molecule_type[0] + '_' + image_t + '_peripheral_fraction.svg'
            plot.bar_profile(periph_fraction, proteins, figname)


def peripheral_dynamic_profile(basic_file_handler,secondary_file_handler, genes,proteins, colors, mrna_tp, protein_tp, path_data,cell_type,image_type):
    if len(image_type) == 0:

        for i in range(len(genes)):
            molecule_type = ['/mrna']
            tp_start = int(mrna_tp[0].split('h')[0])
            tp_end = int(mrna_tp[len(mrna_tp) - 1].split('h')[0])
            x_mrna = np.arange(tp_start, tp_end, 0.01)
            fig = plt.figure()
            ax = plt.axes()
            ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

            mrna_data = np.zeros((3, len(mrna_tp)))
            counter = 0

            for timepoint in mrna_tp:
                print(genes[i], '_', timepoint)
                image_list = helps.preprocess_image_list3(secondary_file_handler, molecule_type, genes[i], [timepoint])
                periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD,path_data)
                periph_frac_median = np.median(periph_frac)
                print(periph_frac_median)
                err = np.median(np.abs(np.tile(np.median(periph_frac), (1, len(periph_frac))) - periph_frac))
                upp_env = periph_frac_median + err
                low_env = periph_frac_median - err
                mrna_data[0, counter] = periph_frac_median
                mrna_data[1, counter] = upp_env
                mrna_data[2, counter] = low_env
                counter += 1
            mrna_data = mrna_data / np.mean(mrna_data[0, :])
            num_tp=[]
            for tp in mrna_tp:
                num_tp.append(int(tp.split('h')[0]))
            print(num_tp)
            print(mrna_data[0, :])


            spl = interpolate.UnivariateSpline(num_tp, mrna_data[0, :],k=len(num_tp)-1)
            spl_upp = interpolate.UnivariateSpline(num_tp, mrna_data[1, :],k=len(num_tp)-1)
            spl_low = interpolate.UnivariateSpline(num_tp, mrna_data[2, :],k=len(num_tp)-1)
            m_y_new = spl(x_mrna)
            m_y_new_upp = spl_upp(x_mrna)
            m_y_new_down = spl_low(x_mrna)
            solid_mrna, = plt.plot(x_mrna, m_y_new, linestyle="-", color=colors[i], linewidth=2)
            plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=colors[i])
            plt.plot(x_mrna, m_y_new_down, linestyle="-", color=colors[i])
            ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=colors[i], alpha=0.5, interpolate=False)

            if genes[i] in proteins:
                molecule_type = ['/protein']
                tp_start=int(protein_tp[0].split('h')[0])
                tp_end=int(protein_tp[len(protein_tp)-1].split('h')[0])
                counter = 0
                x_protein = np.arange(tp_start, tp_end, 0.01)
                protein_data = np.zeros((3, len(protein_tp)))
                for timepoint in protein_tp:
                    print(genes[i], '_', timepoint)
                    image_list = helps.preprocess_image_list3(secondary_file_handler, ['/protein'], genes[i], [timepoint])
                    periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD,path_data)
                    periph_frac_median = np.median(periph_frac)
                    print("periph frac:" + str(periph_frac_median))
                    err = np.median(np.abs(np.tile(np.median(periph_frac), (1, len(periph_frac))) - periph_frac))
                    upp_env = periph_frac_median + err
                    low_env = periph_frac_median - err
                    protein_data[0, counter] = periph_frac_median
                    protein_data[1, counter] = upp_env
                    protein_data[2, counter] = low_env
                    counter += 1
                protein_data = protein_data / np.mean(protein_data[0, :])
                print(protein_data)
                num_tp = []
                for tp in protein_tp:
                    num_tp.append(int(tp.split('h')[0]))
                spl = interpolate.UnivariateSpline(num_tp, protein_data[0, :],k=len(num_tp)-1)
                spl_upp = interpolate.UnivariateSpline(num_tp, protein_data[1, :],k=len(num_tp)-1)
                spl_low = interpolate.UnivariateSpline(num_tp, protein_data[2, :],k=len(num_tp)-1)
                p_y_new = spl(x_protein)
                p_y_new_upp = spl_upp(x_protein)
                p_y_new_down = spl_low(x_protein)

                dashed_protein, = plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=colors[i],
                                           linewidth=2)
                plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=colors[i])
                plt.plot(x_protein, p_y_new_down, linestyle="--", color=colors[i])
                ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=colors[i], alpha=0.25,
                                interpolate=False)
                plt.legend([dashed_protein, solid_mrna], ['Protein', 'Mrna'])
            else:
                plt.legend([solid_mrna], ['Mrna'])

            ax.set_xlim(0, tp_end+1)
            ax.set_ylim(0, 2.0)
            ax.set_ylabel('Peripheral fraction')
            ax.set_xlabel('Time(hrs)')
            ax.set_title(genes[i])

            plt.savefig(
                path.analysis_dir + 'analysis_nocodazole/figures/'+cell_type+'/peripheral_fraction_' + genes[i] + '.svg')

            plt.show()
    else:
        for image_t in image_type:
            for i in range(len(genes)):
                molecule_type = ['/mrna']
                tp_start = int(mrna_tp[0].split('h')[0])
                tp_end = int(mrna_tp[len(mrna_tp) - 1].split('h')[0])
                x_mrna = np.arange(tp_start, tp_end, 0.01)
                fig = plt.figure()
                ax = plt.axes()
                ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

                mrna_data = np.zeros((3, len(mrna_tp)))
                counter = 0

                for timepoint in mrna_tp:
                    print(genes[i], '_', timepoint)
                    image_list = helps.preprocess_image_list4(secondary_file_handler, molecule_type, genes[i],
                                                              [timepoint],image_t)
                    periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,
                                                               constants.PERIPHERAL_FRACTION_THRESHOLD, path_data)
                    periph_frac_median = np.median(periph_frac)
                    print(periph_frac_median)
                    err = np.median(np.abs(np.tile(np.median(periph_frac), (1, len(periph_frac))) - periph_frac))
                    upp_env = periph_frac_median + err
                    low_env = periph_frac_median - err
                    mrna_data[0, counter] = periph_frac_median
                    mrna_data[1, counter] = upp_env
                    mrna_data[2, counter] = low_env
                    counter += 1
                mrna_data = mrna_data / np.mean(mrna_data[0, :])
                num_tp = []
                for tp in mrna_tp:
                    num_tp.append(int(tp.split('h')[0]))
                print(num_tp)
                print(mrna_data[0, :])

                spl = interpolate.UnivariateSpline(num_tp, mrna_data[0, :], k=len(num_tp) - 1)
                spl_upp = interpolate.UnivariateSpline(num_tp, mrna_data[1, :], k=len(num_tp) - 1)
                spl_low = interpolate.UnivariateSpline(num_tp, mrna_data[2, :], k=len(num_tp) - 1)
                m_y_new = spl(x_mrna)
                m_y_new_upp = spl_upp(x_mrna)
                m_y_new_down = spl_low(x_mrna)
                solid_mrna, = plt.plot(x_mrna, m_y_new, linestyle="-", color=colors[i], linewidth=2)
                plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=colors[i])
                plt.plot(x_mrna, m_y_new_down, linestyle="-", color=colors[i])
                ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=colors[i], alpha=0.5, interpolate=False)

                if genes[i] in proteins:
                    molecule_type = ['/protein']
                    tp_start = int(protein_tp[0].split('h')[0])
                    tp_end = int(protein_tp[len(protein_tp) - 1].split('h')[0])
                    counter = 0
                    x_protein = np.arange(tp_start, tp_end, 0.01)
                    protein_data = np.zeros((3, len(protein_tp)))
                    for timepoint in protein_tp:
                        print(genes[i], '_', timepoint)
                        image_list = helps.preprocess_image_list4(secondary_file_handler, molecule_type, genes[i],
                                                                  [timepoint],image_t)
                        periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler,
                                                                    constants.PERIPHERAL_FRACTION_THRESHOLD,
                                                                   path_data)
                        periph_frac_median = np.median(periph_frac)
                        print("periph frac:" + str(periph_frac_median))
                        err = np.median(np.abs(np.tile(np.median(periph_frac), (1, len(periph_frac))) - periph_frac))
                        upp_env = periph_frac_median + err
                        low_env = periph_frac_median - err
                        protein_data[0, counter] = periph_frac_median
                        protein_data[1, counter] = upp_env
                        protein_data[2, counter] = low_env
                        counter += 1
                    protein_data = protein_data / np.mean(protein_data[0, :])
                    print(protein_data)
                    num_tp = []
                    for tp in protein_tp:
                        num_tp.append(int(tp.split('h')[0]))
                    spl = interpolate.UnivariateSpline(num_tp, protein_data[0, :], k=len(num_tp) - 1)
                    spl_upp = interpolate.UnivariateSpline(num_tp, protein_data[1, :], k=len(num_tp) - 1)
                    spl_low = interpolate.UnivariateSpline(num_tp, protein_data[2, :], k=len(num_tp) - 1)
                    p_y_new = spl(x_protein)
                    p_y_new_upp = spl_upp(x_protein)
                    p_y_new_down = spl_low(x_protein)

                    dashed_protein, = plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=colors[i],
                                               linewidth=2)
                    plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=colors[i])
                    plt.plot(x_protein, p_y_new_down, linestyle="--", color=colors[i])
                    ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=colors[i], alpha=0.25,
                                    interpolate=False)
                    plt.legend([dashed_protein, solid_mrna], ['Protein', 'Mrna'])
                else:
                    plt.legend([solid_mrna], ['Mrna'])

                ax.set_xlim(0, tp_end + 1)
                ax.set_ylim(0, 2.00)
                ax.set_ylabel('Peripheral fraction')
                ax.set_xlabel('Time(hrs)')
                ax.set_title(genes[i])

                plt.savefig(
                    path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/peripheral_fraction_' +image_t+
                    '_'+genes[i] + '.svg')

                #plt.show()


if __name__ == "__main__":


    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    ''' 
    1-You need to create a password.txt file before running to connect via ssh
    '''

    basic_file_basename='basic'
    sec_file_basename = 'secondary'
    path_data = path.raw_data_dir

    # basic_file_basename = 'basics_scratch_data'
    # sec_file_basename = 'secondary_scratch_data'
    # path_data=path.scratch_data_dir


    basic_file_path = path.analysis_data_dir + basic_file_basename+'.h5'
    secondary_file_path = path.analysis_data_dir + sec_file_basename+ '.h5'


    ## Figures and analysis
    molecule_type = ['/mrna']
    colors = ['#F16c1b', '#f1bc1b']
    genes = ["pard3", "pard3_nocodazole"]
    proteins=["pard3", "pard3_nocodazole"]
    gene_root_name = ""
    timepoints = ["3h","5h"]
    cell_type = 'micropatterned'
    image_type = []


    with h5py.File(basic_file_path, "r") as basic_file_handler, h5py.File(secondary_file_path, "r") as secondary_file_handler:

        # # Section to build peripheral profile for each timepoint separately
        #peripheral_profile_by_timepoint(secondary_file_handler,molecule_type,genes,timepoints,colors,cell_type,image_type)
        #
        # # Section to build peripheral profile for all timepoint
        #peripheral_profile_by_gene(secondary_file_handler, molecule_type, genes,colors,cell_type,gene_root_name,image_type)
        # #
        # # # # Section to build peripheral profile fraction 10 and 30
        peripheral_fraction_profile(secondary_file_handler, molecule_type, genes, 10,colors, image_type,gene_root_name,basic_file_handler)

        peripheral_fraction_profile(secondary_file_handler, molecule_type, genes, 30,colors, image_type,gene_root_name,basic_file_handler)
        #
        # # # Section to build histogram and plot for specific fraction
        #peripheral_fractions_profile(secondary_file_handler, molecule_type, genes,image_type)
        # #
        # # # Section to compute bar plot peripheral fraction
        histogram_peripheral_profile(basic_file_handler, secondary_file_handler, molecule_type, genes, proteins, colors, path_data,cell_type,image_type)
        # #
        # # # Section to produce plot interpolation (dynamic profile) of peripheral fraction by timepoint
        #peripheral_dynamic_profile(basic_file_handler, secondary_file_handler, genes, proteins, colors, timepoints, timepoints_protein, path_data,cell_type,image_type)
        #
        #
        #
        #
    molecule_type = ['/protein']
    proteins = ["arhgdia", "arhgdia_nocodazole"]
            #"arhgdia_cytod", "pard3", "pard3_nocodazole", "pard3_cytod"]
    gene_root_name = ""
    timepoints = ["3h","5h"]
    #timepoints_protein = ["3h","5h"]
    cell_type = 'micropatterned'
    image_type = []

    with h5py.File(basic_file_path, "r") as basic_file_handler, h5py.File(secondary_file_path,"r") as secondary_file_handler:
        # # Section to build peripheral profile for each timepoint separately
        # peripheral_profile_by_timepoint(secondary_file_handler,molecule_type,genes,timepoints,colors,cell_type,image_type)
        #
        # # Section to build peripheral profile for all timepoint
        # peripheral_profile_by_gene(secondary_file_handler, molecule_type, genes,colors,cell_type,gene_root_name,image_type)
        # #
        # # # # Section to build peripheral profile fraction 10 and 30
        peripheral_fraction_profile(secondary_file_handler, molecule_type, proteins, 10, colors, image_type,gene_root_name, basic_file_handler)

        peripheral_fraction_profile(secondary_file_handler, molecule_type, proteins, 30,colors, image_type,gene_root_name,basic_file_handler)
        #
        # # # Section to build histogram and plot for specific fraction
        # peripheral_fractions_profile(secondary_file_handler, molecule_type, genes,image_type)
        # #
        # # # Section to compute bar plot peripheral fraction
        #histogram_peripheral_profile(basic_file_handler, secondary_file_handler, molecule_type, genes, proteins, colors, path_data,cell_type,image_type)
        # #
        # # # Section to produce plot interpolation (dynamic profile) of peripheral fraction by timepoint
        # peripheral_dynamic_profile(basic_file_handler, secondary_file_handler, genes, proteins, colors, timepoints, timepoints_protein, path_data,cell_type,image_type)
        #
        #
        #
        #
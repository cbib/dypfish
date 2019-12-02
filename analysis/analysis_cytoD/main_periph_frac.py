#!/usr/bin/python
# encoding: UTF-8

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
import src.constants as constants
import src.acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, check_dir, loadconfig, plot_colors_cytoD


def peripheral_profile(file_handler, molecule_type, genes, colors, cell_type, timepoints=[False]):
    for timepoint in timepoints:
        print(timepoint)
        mean_profiles = []
        for gene in genes:
            if timepoint:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, gene, [timepoint])
            else:
                image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            profiles=adsc.compute_peripheral_fraction_profiles_2D(file_handler, image_list)
            mean_profiles.append(np.average(profiles, axis=0))
        plot_peripheral_profile(mean_profiles, timepoint, genes, colors, cell_type)

def plot_peripheral_profile(mean_profiles, timepoint, genes, colors, cell_type):
    timepoint = timepoint if timepoint else 'All'
    plot_filename = 'peripheral_fraction_' + timepoint + '_timepoint' + str(constants.NUM_CONTOURS) + 'contours.png'
    mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[0], len(genes), 1)
    figure_title = ' peripheral profile ('+ timepoint+')'
    stan.plot_profile(mean_profiles, genes, constants.NUM_CONTOURS,check_dir(path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/')+ plot_filename, figure_title, colors, True)

def peripheral_profile_by_timepoint(file_handler,molecule_type,genes,timepoints,colors,cell_type,image_type=[]):
    if len(image_type)==0:
        for timepoint in timepoints:
            mean_profiles = []
            for gene in genes:
                image_list = helps.preprocess_image_list3(file_handler, molecule_type, gene, [timepoint])
                profiles = adsc.compute_peripheral_fraction_profiles_2D(file_handler, image_list)
                mean_profiles.append(np.average(profiles, axis=0))
            plot_name = path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/peripheral_fraction_' + timepoint + '_timepoint' + str(
                constants.NUM_CONTOURS) + 'contours.svg'
            mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
            figure_title = 'timepoint: ' + timepoint
            stan.plot_profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title,colors, True)
    else:
        for image_t in image_type:
            for timepoint in timepoints:
                mean_profiles = []
                for gene in genes:
                        image_list = helps.preprocess_image_list4(file_handler, molecule_type, gene, [timepoint],image_t)
                        profiles = adsc.compute_peripheral_fraction_profiles(file_handler, image_list)
                        mean_profiles.append(np.average(profiles, axis=0))
                plot_name = path.analysis_dir + 'analysis_cytoD/figures/' +cell_type + '/peripheral_fraction_' + timepoint + '_timepoint' + '_'+image_t+'_'+ str(constants.NUM_CONTOURS) + 'contours.svg'
                mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
                figure_title = 'timepoint: ' + timepoint
                stan.plot_profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)


def peripheral_profile_by_gene(file_handler,molecule_type,genes,colors,cell_type,gene_root_name,image_type=[]):
    if len(image_type) == 0:
        mean_profiles = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            profiles = adsc.compute_peripheral_fraction_profiles_2D(file_handler, image_list)
            mean_profiles.append(np.average(profiles, axis=0))
        plot_name = path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/'+gene_root_name+'_peripheral_fraction_all_timepoint' + str(constants.NUM_CONTOURS) +'contours.png'
        mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[1], len(genes), 1)
        figure_title = 'timepoint: All'
        print(mean_profiles)
        stan.plot_profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)
    else:
        mean_profiles_by_conditions=[]
        for image_t in image_type:
            mean_profiles = []
            for gene in genes:
                image_list = helps.preprocess_image_list5(file_handler, molecule_type[0], gene,image_t)
                profiles = adsc.compute_peripheral_fraction_profiles(file_handler, image_list)
                mean_profiles.append(np.average(profiles, axis=0))
            mean_profiles_by_conditions.append(mean_profiles)
            plot_name = path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'contours.svg'
            figure_title = 'timepoint: All'
            print(mean_profiles)
            stan.plot_profile(mean_profiles, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)

        norm_LPA=[]
        norm_plus_medium=[]
        for i in range(len(mean_profiles_by_conditions[0])):
            norm_LPA.append(mean_profiles_by_conditions[1][i] / mean_profiles_by_conditions[0][i])
            norm_plus_medium.append(mean_profiles_by_conditions[2][i] / mean_profiles_by_conditions[0][i])
        plot_name = path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'LPAcontours.svg'
        figure_title = 'timepoint: All'
        stan.plot_profile(norm_LPA, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)
        plot_name = path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/peripheral_fraction_all_timepoint' +image_t+'_'+ str(
                constants.NUM_CONTOURS) + 'plusmediumcontours.svg'
        figure_title = 'timepoint: All'
        stan.plot_profile(norm_plus_medium, genes, constants.NUM_CONTOURS, plot_name, figure_title, colors, True)

def peripheral_fraction_profile(sec_file_handler,molecule_type,genes,fraction,colors,cell_type,image_type,gene_root_name,basic_file_handler,path_data):
    if len(image_type) == 0:
        fractions = []
        for gene in genes:
            print(gene)
            image_list = helps.preprocess_image_list2(sec_file_handler, molecule_type[0], gene)
            fractions.append(adsc.build_histogram_periph_fraction(sec_file_handler, image_list, fraction,path_data,basic_file_handler))
        print(fractions)
        figname = path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/'+molecule_type[0]+'_peripheral_fraction_'+str(fraction)+'.png'
        stan.plot_fraction_profile(fractions, fraction, genes, figname,colors)
    else:
        for image_t in image_type:
            fractions = []
            for gene in genes:
                print(gene)
                image_list = helps.preprocess_image_list5(sec_file_handler, molecule_type[0], gene,image_t)
                fractions.append(adsc.build_histogram_periph_fraction(sec_file_handler, image_list, fraction,path_data,basic_file_handler))
            figname = path.analysis_dir + 'analysis_cytoD/figures/peripheral_fraction/' + cell_type +'/'+gene_root_name+'_peripheral_fraction_' +image_t+'_'+ str(fraction) + '.png'
            stan.plot_fraction_profile(fractions, fraction, genes, figname, colors)

def peripheral_fractions_profile(file_handler, molecule_type, genes,image_type):
    for i in range(1,11):
        gene_counter = 0
        arrs=[]
        for gene in genes:
            image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
            gene_counter += 1
            arr = adsc.compute_fraction_profile(file_handler, image_list, i)
            arrs.append(arr)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 2, 1)
        ax1 = fig.add_subplot(1, 2, 2)
        plt.title("peripheral fraction " + str(i))
        ax.hist(arrs[0], color='r', bins=50, facecolor='red', alpha=0.75, range=(0, 1),
                label=genes[0])
        ax1.hist(arrs[1], bins=50, color='g', facecolor='green', alpha=0.75, range=(0, 1),
                 label=genes[1])
        ax.legend(loc='upper right')
        ax1.legend(loc='upper right')
        plt.show()
        figname=path.figure_dir + "periph_fraction_hist/arhgdia_gapdh_peripheral_fraction_"+str(i)
        fig.savefig(figname)
        plt.close()


def histogram_peripheral_profile(basic_file_handler,secondary_file_handler, molecule_type, genes, proteins, colors, path_data,cell_type,image_type):
    if len(image_type) == 0:
        periph_fraction = []
        for gene in genes:
            print(gene)
            #image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], gene)
            image_list = helps.preprocess_image_list_1(basic_file_handler, gene)
            periph_fraction.append(adsc.compute_mrna_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
        print(periph_fraction)
        figname = path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/'+molecule_type[
            0] +'_peripheral_fraction.png'
        stan.plot_bar_profile(periph_fraction, genes, 0.3, 'Peripheral fraction', figname,colors)
        molecule_type = ['/protein']
        periph_fraction = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(basic_file_handler, molecule_type[0], protein)
            periph_fraction.append(adsc.compute_protein_periph_fraction(image_list,basic_file_handler, secondary_file_handler,  constants.PERIPHERAL_FRACTION_THRESHOLD,path_data))
        figname = path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/'+molecule_type[0] +'_peripheral_fraction.png'
        stan.plot_bar_profile(periph_fraction, proteins, 0.2, 'Peripheral fraction', figname,colors)
    else:
        for image_t in image_type:
            print(image_t)
            periph_fraction = []
            molecule_type = ['/mrna']
            for gene in genes:
                print(gene)
                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], gene,image_t)
                periph_fraction.append(
                    adsc.compute_periph_fraction(image_list, basic_file_handler, secondary_file_handler,
                                                 constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
            figname = path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/' + \
                      molecule_type[
                          0] +'_'+image_t+ 'peripheral_fraction.svg'
            stan.plot_bar_profile(periph_fraction, genes, 0.3, 'Peripheral fraction', figname, colors)
            molecule_type = ['/protein']
            periph_fraction = []
            for protein in proteins:
                print(protein)
                image_list = helps.preprocess_image_list5(basic_file_handler, molecule_type[0], protein,image_t)
                periph_fraction.append(adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD, path_data))
            figname = path.analysis_dir + 'analysis_nocodazole/figures/' + cell_type + '/' + \
                      molecule_type[0] + '_' + image_t + '_peripheral_fraction.svg'
            stan.plot_bar_profile(periph_fraction, proteins, 0.2, 'Peripheral fraction', figname, colors)

def peripheral_dynamic_profile(basic_file_handler,secondary_file_handler, genes,proteins, colors, mrna_tp, protein_tp, path_data,cell_type,image_type):
    if len(image_type) == 0:
        for i in range(len(genes)):
            molecule_type = ['/mrna']
            tp_start = int(mrna_tp[0].split('h')[0])
            tp_end = int(mrna_tp[len(mrna_tp) - 1].split('h')[0])
            x_mrna = np.arange(tp_start, tp_end, 0.01)
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
                tp_start=int(protein_tp[0].split('h')[0])
                tp_end=int(protein_tp[len(protein_tp)-1].split('h')[0])
                counter = 0
                x_protein = np.arange(tp_start, tp_end, 0.01)
                protein_data = np.zeros((3, len(protein_tp)))
                for timepoint in protein_tp:
                    print(genes[i], '_', timepoint)
                    image_list = helps.preprocess_image_list3(secondary_file_handler, ['/protein'], genes[i], [timepoint])
                    periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD,path_data)
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
                path.analysis_dir + 'analysis_cytoD/figures/'+cell_type+'/peripheral_fraction_' + genes[i] + '.svg')
            plt.show()
    else:
        for image_t in image_type:
            for i in range(len(genes)):
                molecule_type = ['/mrna']
                tp_start = int(mrna_tp[0].split('h')[0])
                tp_end = int(mrna_tp[len(mrna_tp) - 1].split('h')[0])
                x_mrna = np.arange(tp_start, tp_end, 0.01)
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
                        periph_frac = adsc.compute_periph_fraction(image_list,basic_file_handler, secondary_file_handler, constants.PERIPHERAL_FRACTION_THRESHOLD,
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
                    path.analysis_dir + 'analysis_cytoD/figures/' + cell_type + '/peripheral_fraction_' +image_t+
                    '_'+genes[i] + '.svg')

def main():
    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots
    enable_logger()
    configData = loadconfig("cytoD")
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    ## Build peripheral profile plot either for each or for all timepoint
    ## Figures and analysis
    path_data = path.raw_data_dir
    molecule_type = ['/protein']
    image_type = []
    cell_type = cell_type_micropatterned
    with h5py.File(path.basic_file_path, "r") as file_handler,h5py.File(path.secondary_file_path, "r") as secondary_file_handler:
        histogram_peripheral_profile(file_handler, secondary_file_handler, molecule_type, genes, proteins, plot_colors_cytoD, path_data, cell_type, image_type)

if __name__ == "__main__":
    main()

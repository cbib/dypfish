#!/usr/bin/python
# encoding: UTF-8

import matplotlib.pyplot as plt
import h5py
import math
import numpy
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from scipy import interpolate
from src.utils import enable_logger, plot_colors, check_dir

''' 
2-This script is supposed to be ran after compute h_star in dypfish. It produces:
 -bar plot for degree of clustering 
'''
def main():
    # Required descriptors: cell_mask, height_map, zero_level and spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    enable_logger()

    # # produce bar plot for degree of clustering
    # with h5py.File(path.h_star_file_path, "a") as input_file_handler, \
    #         h5py.File(path.mtoc_file_path, "a") as mtoc_file_handler:
    #
    #     # mrna part
    #     molecule_type = ['/mrna']
    #     genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    #     base = math.log(0.5)
    #     mrna_median = []
    #     mrna_err = []
    #     for gene in genes:
    #         image_list = helps.preprocess_image_list2(input_file_handler, molecule_type[0], gene)
    #         dof = adsc.compute_degree_of_clustering(input_file_handler, mtoc_file_handler, image_list)
    #         mrna_median.append(math.log(numpy.median(dof)) - base)
    #         err = numpy.median(numpy.abs(numpy.tile(numpy.median(dof), (1, len(dof))) - dof))
    #         error_median = math.log(numpy.median(dof) + err)
    #         error_median = error_median - math.log(numpy.median(dof)) - base
    #         mrna_err.append(error_median)
    #     fig = plt.figure()
    #     ax = plt.axes()
    #     ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    #     N = len(genes)
    #     ind = numpy.arange(N)
    #     width = 0.35
    #     rects1 = ax.bar(ind, mrna_median, width,
    #                     color=plot_colors,
    #                     yerr=mrna_err,
    #                     error_kw=dict(elinewidth=1, ecolor='black'))
    #     ax.set_xlim(-width, len(ind) + width)
    #     ax.set_ylim(0, 20)
    #     ax.set_ylabel('Degree of clustering(log)')
    #     ax.set_title('Mrna degree of clustering')
    #     ax.set_xticks(ind)
    #     plt.legend([gene for gene in genes], loc='upper right')
    #     ax.legend(rects1, genes)
    #     figname = check_dir(path.analysis_dir + 'analysis_cytoD/figures/') + 'mrna_degree_of_clustering.svg'
    #     fig.savefig(figname)
    #
    #     # protein part
    #     molecule_type = ['/protein']
    #     genes = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    #     base = math.log(0.01)
    #     protein_median = []
    #     protein_err = []
    #     for gene in genes:
    #         image_list = helps.preprocess_image_list2(input_file_handler, molecule_type[0], gene)
    #         dof = adsc.compute_degree_of_clustering(input_file_handler, mtoc_file_handler, image_list)
    #         protein_median.append(math.log(numpy.median(dof)) - base)
    #         err = numpy.median(numpy.abs(numpy.tile(numpy.median(dof), (1, len(dof))) - dof))
    #         protein_err.append(math.log(numpy.median(dof) + err) - math.log(numpy.median(dof)) - base)
    #     fig = plt.figure()
    #     ax = plt.axes()
    #     ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    #     N = len(genes)
    #     ind = numpy.arange(N)
    #     width = 0.35
    #     rects1 = ax.bar(ind, protein_median, width,
    #                     color=plot_colors,
    #                     yerr=protein_err,
    #                     error_kw=dict(elinewidth=1, ecolor='black'))
    #     ax.set_xlim(-width, len(ind) + width)
    #     ax.set_ylim(0, 20)
    #     ax.set_ylabel('Degree of clustering(log)')
    #     ax.set_title('Protein degree of clustering')
    #     ax.set_xticks(ind)
    #     plt.legend([gene for gene in genes], loc='upper right')
    #     ax.legend(rects1, genes)
    #     figname = path.analysis_dir + 'analysis_cytoD/figures/protein_degree_of_clustering.svg'
    #     fig.savefig(figname)
    #
    # # produce plot interpolation of degree of clustering by timepoint
    # with h5py.File(path.h_star_file_path, "a") as input_file_handler, h5py.File(path.mtoc_file_path,
    #                                                                                "a") as mtoc_file_handler:
    #     molecule_type = ['/mrna']
    #     genes = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    #     for i in range(len(genes)):
    #         x_mrna = numpy.arange(2, 5, 0.01)
    #         x_protein = numpy.arange(2, 7, 0.01)
    #         test = numpy.zeros((3, 4))
    #         counter = 0
    #         timepoints = ["2h", "3h", "4h", "5h"]
    #         for timepoint in timepoints:
    #             image_list = helps.preprocess_image_list3(input_file_handler, molecule_type, genes[i], [timepoint])
    #             dof = adsc.compute_degree_of_clustering(input_file_handler, mtoc_file_handler, image_list)
    #             dof_median = numpy.median(dof)
    #             err = numpy.median(numpy.abs(numpy.tile(numpy.median(dof), (1, len(dof))) - dof))
    #             upp_env = dof_median + err
    #             low_env = dof_median - err
    #             test[0, counter] = dof_median
    #             test[1, counter] = upp_env
    #             test[2, counter] = low_env
    #             counter += 1
    #         test = test / numpy.mean(test[0, :])
    #         spl = interpolate.UnivariateSpline([2, 3, 4, 5], test[0, :])
    #         spl_upp = interpolate.UnivariateSpline([2, 3, 4, 5], test[1, :])
    #         spl_low = interpolate.UnivariateSpline([2, 3, 4, 5], test[2, :])
    #         m_y_new = spl(x_mrna)
    #         m_y_new_upp = spl_upp(x_mrna)
    #         m_y_new_down = spl_low(x_mrna)
    #         counter = 0
    #
    #         timepoints = ["2h", "3h", "5h", "7h"]
    #         for timepoint in timepoints:
    #             image_list = helps.preprocess_image_list3(input_file_handler, ['/protein'], genes[i], [timepoint])
    #             dof = adsc.compute_degree_of_clustering(input_file_handler, mtoc_file_handler, image_list)
    #             dof_median = numpy.median(dof)
    #             err = numpy.median(numpy.abs(numpy.tile(numpy.median(dof), (1, len(dof))) - dof))
    #             upp_env = dof_median + err
    #             low_env = dof_median - err
    #             test[0, counter] = dof_median
    #             test[1, counter] = upp_env
    #             test[2, counter] = low_env
    #             counter += 1
    #         test = test / numpy.mean(test[0, :])
    #         spl = interpolate.UnivariateSpline([2, 3, 5, 7], test[0, :])
    #         spl_upp = interpolate.UnivariateSpline([2, 3, 5, 7], test[1, :])
    #         spl_low = interpolate.UnivariateSpline([2, 3, 5, 7], test[2, :])
    #         p_y_new = spl(x_protein)
    #         p_y_new_upp = spl_upp(x_protein)
    #         p_y_new_down = spl_low(x_protein)
    #         fig = plt.figure()
    #         ax = plt.axes()
    #         ax.set_xlim(0, 8)
    #         ax.set_ylim(-2, 8)
    #         dashed_protein, = plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=plot_colors[i],
    #                                    linewidth=2)
    #         plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Mrna", color=plot_colors[i])
    #         plt.plot(x_protein, p_y_new_down, linestyle="--", color=plot_colors[i])
    #         ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=plot_colors[i], alpha=0.25,
    #                         interpolate=False)
    #         solid_mrna, = plt.plot(x_mrna, m_y_new, linestyle="-", color=plot_colors[i], linewidth=2)
    #         plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=plot_colors[i])
    #         plt.plot(x_mrna, m_y_new_down, linestyle="-", color=plot_colors[i])
    #         ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=plot_colors[i], alpha=0.5, interpolate=False)
    #         ax.set_ylabel('Degree of clustering(*)')
    #         ax.set_xlabel('Time(hrs)')
    #         ax.set_title(genes[i])
    #         plt.legend([dashed_protein, solid_mrna], ['Protein', 'Mrna'])
    #         plt.savefig(path.analysis_dir + 'analysis_degree_of_clustering/figures/dof_' + genes[i] + '.svg',
    #                     format='svg')

if __name__ == "__main__":
    main()

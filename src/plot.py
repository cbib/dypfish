from __future__ import print_function

import matplotlib
try:
    import Tkinter
except:
    matplotlib.use('Agg')
    # required on headless Linux servers
import matplotlib.pyplot as plt

from scipy.stats import *
from scipy import interpolate
from utils import plot_colors,check_dir
import src.helpers as helps
import seaborn as sns
from helpers import *


def bar_profile_median(median,genes,molecule_type, figname, err):
    fig = plt.figure()
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    N = len(genes)
    ind = np.arange(N)
    width = 0.35
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, 15)
    ax.set_xticks([])
    fig.savefig(figname)
    plt.close()

def sns_boxplot(dd,my_pal,figname):
    box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants', palette=my_pal)
    box.set_xlabel("", fontsize=15)
    box.set_ylabel("", fontsize=15)
    box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    box.legend_.remove()
    plt.yticks(fontsize=15)
    plt.savefig(figname, format='png')

def bar_profile(data,genes,figname):
    #print(data)
    plt.figure()
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    dataMedians = []
    dataStd = []
    y_lim=0
    y_std = 0
    for l in data:
        max_v=np.median(l)
        max_std = np.std(l)
        if max_v > y_lim:
            y_lim=max_v
            y_std=max_std
        dataMedians.append(np.median(l))
        dataStd.append(np.std(l))
    #print(dataStd)
    ind = np.arange(N)
    width = 0.35
    ax.bar(ind, dataMedians, width,color=plot_colors,
                    yerr=dataStd,
                    error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim + 2* y_std)
    ax.set_xticks(ind)
    ax.set_xticklabels(["" for i in range(0, N)])
    plt.savefig(figname,format='png')

def fraction_profile(fractions,fraction,genes,figname,colors):
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(genes)
    fractionsMeans = []
    fractionsStd = []
    y_lim = 0
    y_std = 0
    for l in fractions:
        max_v = np.mean(l)
        max_std = np.std(l)
        if max_v > y_lim:
            y_lim = max_v
            y_std = max_std
        fractionsMeans.append(np.mean(l))
        fractionsStd.append(np.std(l))
    ind = np.arange(N)
    width = 0.35
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim + 2 * y_std)
    ax.set_ylim(0, 0.2)

    xTickMarks = ["" for i in range(0, len(genes))]
    ax.set_xticks(ind)
    plt.savefig(figname)
    plt.close()

def histogram_noise_measured(nm_arhgdia,nm_arhgdia_cultured):
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    ind = np.arange(1,3)
    width = 0.25
    #ax.set_xlim(-width * 2, len(ind) + width)
    ax.set_ylim(0, 0.5)
    ax.set_title('')
    xTickMarks = ["", ""]
    ax.set_xticks(ind)
    ax.bar(ind, [nm_arhgdia,nm_arhgdia_cultured], width, color=plot_colors)
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/Volume-corrected_noise_measure.png", format='png')

def noise_measured_dynamic_profile(nms, gene, color):
    x_mrna = np.arange(2, 5, 0.01)
    spl = interpolate.UnivariateSpline([2, 3, 4, 5], nms)
    m_y_new = spl(x_mrna)
    fig = plt.figure()
    ax = plt.axes()
    ax.set_xlim(1, 6)
    if gene=="arhgdia":
        y_lim=0.4
    elif gene == "b-actin":
        y_lim=0.15
    elif gene == "gapdh":
        y_lim=0.7
    elif gene == "pard3":
        y_lim=0.35
    elif gene == "pkp4":
        y_lim=0.6
    else:
        y_lim=0.15

    ax.set_ylim(0, y_lim)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=30)
    plt.plot(x_mrna, m_y_new, linestyle="-", color=color, linewidth=2)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ax.set_ylabel('Volume corrected noise measure')
    ax.set_xlabel('Time(hrs)')
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.set_title(gene)
    plt.savefig(path.analysis_dir + "/analysis_spots_density/figures/nm_profile_" + gene+".png", format='png')

def profile(profiles, genes, slice_number, figname,figtitle,colors,save=False):
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=30)
    plt.xticks([w for w in range(0,slice_number+2,10)])
    for i in range(len(genes)):
        print(profiles[i])
        plt.plot(np.arange(slice_number), profiles[i], color=colors[i], linewidth=3, label=genes)
    if save:
        plt.savefig(figname)
    plt.close()


# compare descriptor profile for mrna/protein over time
def dynamic_profiles(mrna_data,protein_data,gene,plot_color,xlabel,ylabel,figpath):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(30)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(30)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.set_xlim(1, 8)
    x_mrna = np.arange(2, 5, 0.01)
    spl = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[0, :])
    spl_upp = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[1, :])
    spl_low = interpolate.UnivariateSpline([2, 3, 4, 5], mrna_data[2, :])
    m_y_new = spl(x_mrna)
    m_y_new_upp = spl_upp(x_mrna)
    m_y_new_down = spl_low(x_mrna)
    plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=plot_color)
    plt.plot(x_mrna, m_y_new_down, linestyle="-", color=plot_color)
    ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=plot_color, alpha=0.5, interpolate=False)

    x_protein = np.arange(2, 7, 0.01)
    spl = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[0, :])
    spl_upp = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[1, :])
    spl_low = interpolate.UnivariateSpline([2, 3, 5, 7], protein_data[2, :])
    p_y_new = spl(x_protein)
    p_y_new_upp = spl_upp(x_protein)
    p_y_new_down = spl_low(x_protein)
    plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=plot_color)
    plt.plot(x_protein, p_y_new_down, linestyle="--", color=plot_color)
    ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=plot_color, alpha=0.25,
                    interpolate=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(gene)
    plt.savefig(figpath,format='png')

def boxplot(vec_a,vec_b,figpath):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.1)
    bp = ax.boxplot([vec_a, vec_b])
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['medians'], color='black', linewidth=3)
    ax.set_xticklabels(['', ''])
    plt.savefig(figpath, format='png')



def process_data(file_handler, timepoints, molecule_type, gene, calc_function, *args):
    # for each timepoint, compute cytoplasmic spread
    # compute median, low and up envelope
    # normalize data by mean of median timepoint result
    tp_l= len(timepoints)
    data = np.zeros((3,tp_l))
    for i in range(tp_l):
        image_list = helps.preprocess_image_list3(file_handler, molecule_type, gene, [timepoints[i]])
        results = calc_function(image_list, *args)
        results_median = np.median(results)
        err = np.median(np.abs(np.tile(np.median(results), (1, len(results))) - results))
        upp_env = results_median + err
        low_env = results_median - err
        data[0, i] = results_median
        data[1, i] = upp_env
        data[2, i] = low_env
    normalized_data = data / np.mean(data[0, :])
    return normalized_data

def data_extractor(genes,proteins, secondary_file_handler, calc_function, *args):
    mrna_tp= ["2h", "3h", "4h", "5h"]
    protein_tp= ["2h", "3h", "5h", "7h"]
    for i in range(len(genes)):
        mrna_data=process_data(secondary_file_handler,mrna_tp,['/mrna'],genes[i], calc_function, *args)
        if genes[i] in proteins:
            protein_data=process_data(secondary_file_handler, protein_tp, ['/protein'], genes[i], calc_function, *args)

        yield mrna_data, protein_data, i

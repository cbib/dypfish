from __future__ import print_function

import traceback
import matplotlib
try:
    import Tkinter
except:
    matplotlib.use('Agg')
    # required on headless Linux servers
import matplotlib.pyplot as plt

from scipy.stats import *
from scipy import interpolate
from utils import plot_colors,check_dir,plot_colors_chx
import src.helpers as helps
import seaborn as sns
from helpers import *



def linear_regression(
        cell_area,
        num_spots,
        graph_file_path_name):

    # display figures
    fig, (ax2, ax1) = plt.subplots(1, 2, sharey=True, figsize=(15, 5))
    xs1 = cell_area[0]
    ys1 = num_spots[0]
    xs2 =  cell_area[1]
    ys2 = num_spots[1]

    # Create linear regression object
    ax1.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.1)
    ax1.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    fit1 = np.polyfit(xs1, ys1, 1)
    p = np.poly1d(fit1)
    fit_fn1 = np.poly1d(fit1)

    ax1.plot(xs1, ys1, 'yo', xs1, fit_fn1(xs1), '--k', c='#0080ff')
    ax2.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.1)
    ax2.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    fit2 = np.polyfit(xs2, ys2, 1)
    p = np.poly1d(fit2)

    fit_fn2 = np.poly1d(fit2)
    ax2.plot(xs2, ys2, 'yo', xs2, fit_fn2(xs2), '--k', c='#646464')
    ax2.axis([200, 1500, -100, 600])
    ax1.axis([200, 800, -100, 600])

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for axis in ['bottom', 'left']:
        ax1.spines[axis].set_linewidth(3)
        ax2.spines[axis].set_linewidth(3)
    print(graph_file_path_name)
    plt.savefig(graph_file_path_name, format="png")
    plt.close()

def boxplot2(data,graph_file_path_name, ylim):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, direction='inout', length=8, width=2, colors='black')
    bp = ax.boxplot(data)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)
    plt.setp(bp['boxes'], color='black', linewidth=2)
    plt.setp(bp['whiskers'], color='black', linewidth=2)
    plt.setp(bp['medians'], color='black', linewidth=2)
    plt.setp(bp['caps'], color='black', linewidth=2)
    plt.setp(bp['fliers'], color='black', linewidth=2)
    axes = plt.gca()
    axes.set_ylim([0, ylim])
    #plt.xticks([1, 2], ['Micropatterned cells', 'Standard cells'])
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.savefig(graph_file_path_name,format="png")
    plt.close()

def sns_linear_regression(data_1,data_2,color,graph_file_path_name):

    sns.set(style="white", color_codes=True)
    annot_kws = {'prop': {'family': 'monospace', 'weight': 'bold', 'size': 8}}
    res1 = pearsonr(data_1, data_2)
    #print(", ".join(["x" + unichr(u) for u in (0x2070, 0x00B9, 0x00B2, 0x00B3, 0x2074, 0x2075, 0x2076, 0x2077, 0x2078, 0x2079)]))
    sns.set(font_scale=1)
    data3=np.array(data_2)/np.array(data_1)
    data3 = data3[data3 <= 0.6]
    data3=data3-np.median(data3)
    #plt.plot(data3,marker='o')
    #plt.scatter(np.arange(1,len(data3)+1),np.array(data3))
    #plt.show()
    #sns.plt.xlim(350, 850)
    #print("data3:",len(data3))
    #sns.distplot(data3, hist=True, kde=True,bins=50, color='darkblue',hist_kws={'edgecolor': 'black'},kde_kws={'linewidth': 4})
    #plt.ylim((0,7))
    #plt.savefig(graph_file_path_name, format="png")
    #plt.close()
    sns_plot_regression = sns.jointplot(x=data_1, y=data_2, kind='reg',color=color)
    #sns_plot_regression.ax_marg_x.set_xlim(350,850)
    phantom, = sns_plot_regression.ax_joint.plot([], [], linestyle="", alpha=0)
    sns_plot_regression.ax_joint.legend([phantom], ['pearsonr={:f}, R'+ unichr(0x00B2)+'={:f}, p={:.2E}'.format(res1[0],res1[0]**2, res1[1])],**annot_kws)
    sns_plot_regression.savefig(graph_file_path_name, format="png")


def bar_profile_median(median,genes, figname, err, colors):
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
    ax.bar(ind, median, width, color=colors,
           yerr=err,
           error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(np.min(median)-np.min(err)-2, np.max(median)+np.max(err)+1)
    ax.set_xticks([])
    fig.savefig(figname)
    plt.close()

def sns_boxplot(dd,my_pal,figname):
    fig = plt.figure()
    box = sns.boxplot(x='Gene', y='value', data=dd, hue='Quadrants', palette=my_pal)
    box.set_xlabel("", fontsize=15)
    box.set_ylabel("", fontsize=15)
    box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    #box.set(ylim=(-3, 3))
    box.legend_.remove()
    plt.yticks(fontsize=15)
    fig.savefig(figname, format='png')

def bar_profile(data,genes,figname,colors):
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
    ax.bar(ind, dataMedians, width,color=colors,
                    yerr=dataStd,
                    error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim + 2* y_std)
    #ax.set_ylim(0, 0.03)
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
    ax.bar(ind, fractionsMeans, width, color=plot_colors,
           yerr=fractionsStd,
           error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, y_lim + 2 * y_std)
    #ax.set_ylim(0, 0.2)

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
    plt.savefig(path.analysis_dir + "/spots_density/figures/Volume-corrected_noise_measure.png", format='png')

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
    plt.savefig(path.analysis_dir + "/spots_density/figures/nm_profile_" + gene+".png", format='png')

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
        #print(profiles[i])
        #print(np.arange(slice_number))
        plt.plot(np.arange(slice_number), profiles[i], color=colors[i], linewidth=3, label=genes)
    if save:
        plt.savefig(figname)
    plt.close()


# compare descriptor profile for mrna/protein over time
def dynamic_profiles(mrna_data,protein_data,timepoints_num_mrna, timepoints_num_protein, gene,plot_color,xlabel,ylabel,figpath):
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
    #ax.set_ylim(-2, 3)
    x_mrna = np.arange(np.min(timepoints_num_mrna), np.max(timepoints_num_mrna), 0.01)
    spl = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[0, :],k=len(timepoints_num_mrna)-1)
    spl_upp = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[1, :],k=len(timepoints_num_mrna)-1)
    spl_low = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[2, :],k=len(timepoints_num_mrna)-1)
    m_y_new = spl(x_mrna)
    m_y_new_upp = spl_upp(x_mrna)
    m_y_new_down = spl_low(x_mrna)
    plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=plot_color)
    plt.plot(x_mrna, m_y_new_down, linestyle="-", color=plot_color)
    ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=plot_color, alpha=0.5, interpolate=False)

    x_protein = np.arange(np.min(timepoints_num_protein), np.max(timepoints_num_protein), 0.01)
    spl = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[0, :],k=len(timepoints_num_mrna)-1)
    spl_upp = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[1, :],k=len(timepoints_num_mrna)-1)
    spl_low = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[2, :],k=len(timepoints_num_mrna)-1)
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
        try:
            mrna_data=process_data(secondary_file_handler,mrna_tp,['/mrna'],genes[i], calc_function, *args)
            if genes[i] in proteins:
                protein_data=process_data(secondary_file_handler, protein_tp, ['/protein'], genes[i], calc_function, *args)
            yield mrna_data, protein_data, i
        except Exception as e:
            print("Got exception : %s" % str(e))
            tr = traceback.format_exc()
            print(tr)
            raise


def data_extractor_generic(genes, proteins, timepoints_m, timepoints_p ,secondary_file_handler, calc_function, *args):
    for i in range(len(genes)):
        try:
            mrna_data=process_data(secondary_file_handler,timepoints_m,['/mrna'],genes[i], calc_function, *args)
            if genes[i] in proteins:
                protein_data = process_data(secondary_file_handler, timepoints_p, ['/protein'], genes[i], calc_function,
                                            *args)
            yield mrna_data,protein_data, i
        except Exception as e:
            print("Got exception : %s" % str(e))
            tr = traceback.format_exc()
            print(tr)
            raise


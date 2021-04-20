# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

# These functions are specific for plotting the results produced by the analysis scripts

import matplotlib
from loguru import logger

import constants
from scipy.interpolate import interp1d
from helpers import create_dir_if_needed_for_filepath
from path import global_root_dir
from mpi_calculator import DensityStats
import math
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate
from scipy.stats import pearsonr
import pandas as pd
import helpers
import pathlib
import numpy as np

pd.set_option('display.max_rows', 1000)


def sns_boxplot(dd, my_pal, figname, x="Gene", y="value", hue='Quadrants'):
    fig = plt.figure()
    box = sns.boxplot(x=x, y=y, data=dd, palette=my_pal)
    box.set_xlabel("", fontsize=15)
    box.set_ylabel("", fontsize=15)
    box.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    box.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    # box.set(ylim=(-3, 3))
    # box.legend_.remove()
    plt.yticks(fontsize=15)
    fig.savefig(figname, format='png')
    plt.close()


def bar_profile_median_timepoints(df : pd.DataFrame, palette, tgt_fp, gene, fixed_yscale=0):
    """
    Plots a barplot for 'd_of_c' column for 2 molecules at each timepoint for gene
    Dataframe df contains:
       - 'd_of_c' column to be plotted
       - 'Molecule' column containing either 'mrna' or 'protein'
       - 'Timepoint' column providing the timepoints where they are measured
       - 'error' column for sem error from the median
       - 'CI' column for the confidence interval
    palette is a dict with keys 'mrna' and 'protein'
    fixed_yscale are for setting the y_scale upper limit
    """
    bar_width = 0.35
    fig, ax = plt.subplots()
    labels = np.sort(df["Timepoint"].unique())
    index = np.arange(len(labels))
    df[['lower','higher']] = pd.DataFrame(df['CI'].to_list(), columns=['lower', 'higher'], index= df.index)
    dfs = df.sort_values('Timepoint')
    # Fixing readability - scaling up for errors and low values (due to shift)
    dfs.loc[(dfs.error < 0.2), "error"] = dfs['error'] * 2
    dfs = dfs.round({'lower': 1, 'higher': 1})

    mrna = dfs[dfs["Molecule"] == "mrna"]
    protein = dfs[dfs["Molecule"] == "protein"]
    mask = mrna['d_of_c'] != 0
    x_new = np.linspace(index[mask].min(), index[mask].max(), 20)
    f_cubic = interp1d(index[mask], mrna['d_of_c'][mask], kind='cubic')
    ax.plot(x_new, f_cubic(x_new), zorder=5, color=helpers.colorscale("A5073E", 0.9))
    mask = protein['d_of_c'] != 0
    x_new = np.linspace(np.min([index[mask]+bar_width]), np.max([index[mask]+bar_width]), 20)
    f_cubic = interp1d(index[mask]+bar_width, protein['d_of_c'][mask], kind='cubic')
    ax.plot(x_new, f_cubic(x_new), zorder=5, color=helpers.colorscale("A5073E", 1.2))

    ax.bar(index, mrna["d_of_c"], bar_width, yerr=mrna["error"], color=palette["mrna"], error_kw=dict(elinewidth=6, ecolor='black'))
    ax.bar(index + bar_width, protein["d_of_c"], bar_width, yerr=protein["error"], color=palette["protein"], error_kw=dict(elinewidth=6, ecolor='black'))

    plt.ylim([0, fixed_yscale])
    ax.set_xlabel("", fontsize=15)
    ax.set_ylabel("", fontsize=15)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    all_timepoints = np.sort(list(set(mrna['Timepoint']) | set(protein['Timepoint'])))
    plt.xticks(index + bar_width/2, all_timepoints, fontsize=20)

    # create the table for the CI
    CI = pd.DataFrame({"mRNA": mrna[["lower", "higher"]].values.tolist(),
                        "protein": protein[["lower", "higher"]].values.tolist()}, index=all_timepoints)
    CI = CI.T
    plt.table(cellText=CI.values,
              rowLabels=CI.index.values,
              colLabels=list(str(' ') * len(all_timepoints)), # empty since we already have the xticks
              loc='bottom',
              bbox=[0.0, -0.3, 1, .28])
    plt.subplots_adjust(bottom=0.28)


    fig.savefig(tgt_fp, format='png')
    plt.close()
    logger.info("Generated image at {}", tgt_fp)


def sns_barplot(dd, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type", err="err"):
    fig = plt.figure()
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ax.spines['left'].set_linewidth(3)
    plt.yticks(fontsize=20)
    idx1 = np.arange(dd["Timepoint"].nunique())
    width = 0.35
    idx2 = [x + width for x in idx1]
    x_labels = sorted(dd["Timepoint"].unique())
    mrna_values = list(dd[dd["Molecule_type"] == "mrna"][y])
    mrna_err = list(dd[dd["Molecule_type"] == "mrna"][err])
    protein_values = list(dd[dd["Molecule_type"] == "protein"][y])
    protein_err = list(dd[dd["Molecule_type"] == "protein"][err])
    protein_err = [x if (float(x) < 0.0000000000001) else 0.05 for x in protein_err]
    mrna_err = [x if (float(x) < 0.0000000000001) else 0.05 for x in mrna_err]
    assert (len(mrna_values) == len(protein_values))
    assert (len(mrna_values) == len(mrna_err))
    ax.bar(idx1, mrna_values, width, color=my_pal["mrna"],
           yerr=mrna_err, error_kw=dict(elinewidth=1, ecolor='black'))
    ax.bar(idx2, protein_values, width, color=my_pal["protein"],
           yerr=protein_err, error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(idx1) + width)
    #    ax.set_ylim(-1, 1)
    ax.set_xticks(idx1 + width / 2)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    ax.set_xticklabels(x_labels)
    fig.savefig(tgt_fp)
    plt.close()


def sns_violinplot(dd, my_pal, tgt_fp, x="Gene", y='value', hue=None, no_legend=True, rotation=0):
    fig = plt.figure(figsize=(10, 10))
    # medianprops = dict(linestyle='--', color='firebrick')
    # sns.boxplot(x=x, y=y, data=dd, hue=hue, showfliers=False, showcaps=True, whis=(25.0, 75.0),showbox=False, meanline=True, showmeans=True, linewidth=0.8, medianprops=medianprops)
    ax = sns.violinplot(x=x, y=y, data=dd, hue=hue, palette=my_pal)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    # ax.legend(["",""])
    # box.set(ylim=(0, 8))
    plot_xlabels = constants.analysis_config['MRNA_GENES_LABEL']
    ax.set_xticklabels(plot_xlabels, rotation=rotation)
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=20)
    if no_legend:
        ax.legend_.remove()
    fig.savefig(tgt_fp, format='png')
    plt.close()


def bar_profile(data, genes, tgt_fp):
    plot_colors = constants.analysis_config['PLOT_COLORS']
    plt.figure(figsize=(8, 8))
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    dataMedians = []
    dataStdErr = []
    for l in data:
        dataMedians.append(np.median(l))
        dataStdErr.append(np.std(l)/math.sqrt(len(l)))
    ind = np.arange(len(genes))
    width = 0.35
    ax.bar(ind, dataMedians, width, color=plot_colors,
           yerr=dataStdErr,
           error_kw=dict(elinewidth=1, ecolor='black'))
    ax.set_xlim(-width, len(ind) + width)
    ax.set_xticks(ind)
    ax.set_xticklabels(["" for i in range(0, len(genes))])
    create_dir_if_needed_for_filepath(figname)
    plt.savefig(tgt_fp, format='png')
    plt.close()


def bar_profile_median(genes, medians, err, CI, molecule_type, fixed_yscale, tgt_fp):
    """
    Plot a barplot for each gene with height given by medians, error bars defined by err
    and confidence intervals by CI; CI is a dictionary with keys = genes
    """
    plot_colors = constants.analysis_config['PLOT_COLORS']
    fig, ax = plt.subplots()

    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    ax.spines['left'].set_linewidth(3)
    plt.yticks(fontsize=20)
    N = len(genes)
    ind = np.arange(N)
    width = 0.35
    ax.bar(ind, medians, width, color=plot_colors, yerr=err, error_kw=dict(elinewidth=6, ecolor='black'))

    # add confidence intervals
    # TODO here add the table for CI in the same way as in timepoints

    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, fixed_yscale) # same as for the timepoints
    ax.set_xticks([])

    fig.savefig(tgt_fp)
    plt.close()
    logger.info("Generated image at {}", tgt_fp)


def bar_profile_simple(data, tgt_fp, plot_colors):
    plt.figure()
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    N = len(data)
    ind = np.arange(N)
    width = 0.35
    ax.bar(ind, data, width, color=plot_colors)
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, 0.7)
    ax.set_xticks(ind)
    ax.set_xticklabels(["" for i in range(0, N)])
    plt.savefig(tgt_fp, format='png')
    plt.close()


def plot_MPI(density_stats: DensityStats, molecule_type, tgt_fp):
    labels = density_stats.make_labels()
    mpis, errs = density_stats.mpi()
    bar_profile_median(mpis, labels, errs, tgt_fp)


def plot_boxplot_DOC(mrna_data, protein_data, gene, tgt_fp,plot_colors):
    """
    The timepoint_list has to have the same order as in the mpi() function
    """

    print(mrna_data)
    # plot_colors = constants.analysis_config['PLOT_COLORS']
    # mrna_density_stats.update_group_key(['Timepoint'])
    # protein_density_stats.update_group_key(['Timepoint'])
    #
    # for color_num, gene in enumerate(molecule_list):
    #     # Collect mrna and protein density statistics for this gene only
    #     dd = {'Molecule_type': [], 'Timepoint': [], 'MPI': [], 'err': []}
    #     gene_stats = mrna_density_stats.subset_stats("Gene", gene)
    #     prot_stats = protein_density_stats.subset_stats("Gene", gene)
    #     mrna_mpis, mrna_errs = gene_stats.mpi()
    #     prot_mpis, prot_errs = prot_stats.mpi()
    #     dd["MPI"] = mrna_mpis + prot_mpis
    #     dd["err"] = mrna_errs + prot_errs
    #     dd["Molecule_type"] = ["mrna"] * len(mrna_timepoint_list) + ["protein"] * len(protein_timepoint_list)
    #     dd["Timepoint"] = sorted(mrna_timepoint_list) + sorted(protein_timepoint_list)
    #     for tp in np.setdiff1d(mrna_timepoint_list, protein_timepoint_list):
    #         dd["Timepoint"].append(tp);
    #         dd["Molecule_type"].append("protein")
    #         dd["MPI"].append(0);
    #         dd["err"].append(0)
    #     for tp in np.setdiff1d(protein_timepoint_list, mrna_timepoint_list):
    #         dd["Timepoint"].append(tp);
    #         dd["Molecule_type"].append("mrna")
    #         dd["MPI"].append(0);
    #         dd["err"].append(0)
    #
    #     df = pd.DataFrame(dd)
    #     df.sort_values("Timepoint", axis=0, ascending=True, inplace=True)
    #
    #     tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_DYNAMIC_MPI'].format(gene=gene)
    #     tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
    #                           tgt_image_name)
    #
    #     my_pal = {"mrna": str(plot_colors[color_num]),
    #               "protein": str(helpers.color_variant(plot_colors[color_num], +80))}
    #     helpers.create_dir_if_needed_for_filepath(tgt_fp)
    #     sns_barplot(df, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type", err="err")
    #     # sns_barplot_simple(df, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type")
    #     logger.info("Generated image at {}", tgt_fp)

# plot_mtoc_enrichment never used ??
def plot_mtoc_enrichment(density_stats: DensityStats, molecule_type, limit_threshold, log=False):
    # melt dataframe and relabel all non MTOC quadrants
    value_vars = density_stats.quadrant_labels + [density_stats.mtoc_quadrant_label]
    dd = pd.melt(density_stats.df, id_vars=['Gene'], value_vars=value_vars, var_name='Quadrants')
    for label in density_stats.quadrant_labels:
        dd = dd.replace(label, 'Non MTOC')
    dd = dd.replace(0.000000, np.nan)  # TODO why this is here? there should be no np.nan
    # apply log scale
    if (log):
        dd['value'] = dd['value'].apply(np.log2)
    # Choose color palette
    my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#1a8cff"}

    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    ## remove outliers
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]

    helpers.create_dir_if_needed_for_filepath(tgt_fp)
    sns_violinplot(dd, my_pal, tgt_fp, rotation=45)
    logger.info("Generated image at {}", tgt_fp)

# plot_hist_ratio never used ??
def plot_hist_ratio(density_stats: DensityStats, molecule_type, limit_threshold, groupby=['Gene']):
    df = density_stats.df
    df['MTOC ratio'] = density_stats.ratios()
    dd = pd.melt(df, id_vars=groupby, value_vars=['MTOC ratio'], var_name='Quadrants')
    dd = dd.replace(0.000000, np.nan)
    tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT_RATIO'].format(molecule_type=molecule_type)
    tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                          tgt_image_name)
    my_pal = {"MTOC ratio": "#66b2ff"}
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]  # dd[dd["value"] < limit_threshold]

    helpers.create_dir_if_needed_for_filepath(tgt_fp)
    sns_violinplot(dd, my_pal, tgt_fp, x=groupby[0])
    logger.info("Generated image at {}", tgt_fp)


def compute_violin_plot_ratio(density_stats: DensityStats, molecule_type, tgt_fp, limit_threshold=6, groupby=['Gene'], term=""):
    df = density_stats.df
    if term != "":
        df[groupby[0]] = df.apply(lambda row: term if term[2:len(term)] in row[groupby[0]] else row[groupby[0]], axis=1)
    groups = [tp for tp, line in df.groupby(groupby, sort=False)]
    df["MTOC_ratio"] = df.apply(lambda row: row['MTOC'] / row.filter(regex=("Non MTOC.*")).mean() if row.filter(
        regex=("Non MTOC.*")).mean() != 0 else 0.0, axis=1)
    dd = pd.melt(df, id_vars=groupby[0], value_vars=["MTOC_ratio"], var_name='Quadrants')
    dd = dd.replace(0.000000, np.nan)
    dd.dropna(inplace=True)

    my_pal = {}
    for i, color in enumerate(constants.analysis_config['PLOT_COLORS'][0:len(groups)]):
        my_pal[groups[i]] = color
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]  # dd[dd["value"] < limit_threshold]

    sns_violinplot(dd, my_pal, tgt_fp, x=groupby[0], no_legend=False, rotation=45)


def compute_categorical_violin_plot_ratio(density_stats: DensityStats, molecule_type, tgt_fp, limit_threshold=6,
                                          groupby=['Gene'], term="", gene=""):
    df = density_stats.df
    if term != "":
        df["Label"] = df.apply(
            lambda row: term if term[2:len(term) - 1] in row[groupby[len(groupby) - 1]] else 'Control',
            axis=1)

    groups = [tp for tp, line in df.groupby(["Label"], sort=False)]
    df["MTOC_ratio"] = df.apply(lambda row: row['MTOC'] / row.filter(regex=("Non MTOC.*")).mean() if row.filter(
        regex=("Non MTOC.*")).mean() != 0 else 0.0, axis=1)
    dd = pd.melt(df, id_vars=groupby[len(groupby) - 1], value_vars=["MTOC_ratio"], var_name='Quadrants')
    if term != "":
        dd["Quadrants"] = df["Label"].values
    if gene != "":
        dd["Gene"] = [gene for i in range(len(df["Label"].values))]
    dd = dd.replace(0.000000, np.nan)
    dd.dropna(inplace=True)

    my_pal = {}
    for i, color in enumerate(constants.analysis_config['PLOT_COLORS'][0:len(groups)]):
        my_pal[groups[i]] = color
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]  # dd[dd["value"] < limit_threshold]
    sns_violinplot(dd, my_pal, tgt_fp, x="Gene", hue="Quadrants")


def compute_violin_plot_enrichment(density_stats: DensityStats, molecule_type, tgt_fp, limit_threshold=6, log=False,
                                   groupby="Gene"):
    df = density_stats.df
    ## melt dataframe and group together all non MTOC quadrant

    value_vars = [density_stats.mtoc_quadrant_label] + density_stats.quadrant_labels
    dd = pd.melt(df, id_vars=[groupby], value_vars=value_vars,
                 var_name='Quadrants')
    for label in density_stats.quadrant_labels:
        dd = dd.replace(label, 'Non MTOC')
    dd = dd.replace(0.000000, np.nan)
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]
    dd.dropna(inplace=True)
    ## apply log scale
    if (log):
        dd['value'] = dd['value'].apply(np.log2)
    ## Choose color palette
    my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#1a8cff"}

    ## remove outliers
    sns_violinplot(dd, my_pal, tgt_fp, x=groupby, hue=dd['Quadrants'], rotation=45)


def compute_categorical_violin_plot_enrichment(density_stats: DensityStats, molecule_type, tgt_fp, limit_threshold=8, log=False,
                                               groupby="Gene", term="", gene=""):
    df = density_stats.df
    value_vars = [density_stats.mtoc_quadrant_label] + density_stats.quadrant_labels

    ## melt dataframe and group together all non MTOC quadrant
    groups = [tp for tp, line in df.groupby(groupby, sort=False)]
    if term != "":
        df[groupby] = df.apply(lambda row: term if term[2:len(term)] in row[groupby] else row["Gene"], axis=1)

    dd = pd.melt(df, id_vars=[groupby], value_vars=value_vars,
                 var_name='Quadrants')
    for label in density_stats.quadrant_labels:
        dd = dd.replace(label, 'Non MTOC')
    if gene != "":
        dd["Gene"] = [gene for i in range(len(df["Label"].values))]
    dd = dd.replace(0.000000, np.nan)
    outliers = helpers.detect_outliers(np.array(dd["value"]), limit_threshold)
    dd = dd[~np.isin(dd["value"], outliers)]  # dd[dd["value"] < limit_threshold]
    dd.dropna(inplace=True)
    ## apply log scale
    if (log):
        dd['value'] = dd['value'].apply(np.log2)
    ## Choose color palette
    my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#1a8cff"}

    ## remove outliers
    sns_violinplot(dd, my_pal, tgt_fp, x=groupby, hue=dd['Quadrants'])


# compare descriptor profile for mrna/protein over time
def dynamic_profiles(mrna_data, protein_data, gene, xlabel, ylabel, figpath, plot_colors):
    timepoints_num_mrna = constants.dataset_config['TIMEPOINTS_NUM_MRNA']
    timepoints_num_protein = constants.dataset_config['TIMEPOINTS_NUM_PROTEIN']
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
    # ax.set_ylim(-2, 8)
    x_mrna = np.arange(np.min(timepoints_num_mrna), np.max(timepoints_num_mrna), 0.01)
    spl = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[0, :], k=len(timepoints_num_mrna) - 1)
    spl_upp = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[1, :], k=len(timepoints_num_mrna) - 1)
    spl_low = interpolate.UnivariateSpline(timepoints_num_mrna, mrna_data[2, :], k=len(timepoints_num_mrna) - 1)
    m_y_new = spl(x_mrna)
    m_y_new_upp = spl_upp(x_mrna)
    m_y_new_down = spl_low(x_mrna)
    plt.plot(x_mrna, m_y_new, linestyle="-", color=plot_colors)
    plt.plot(x_mrna, m_y_new_upp, linestyle="-", color=plot_colors)
    plt.plot(x_mrna, m_y_new_down, linestyle="-", color=plot_colors)
    ax.fill_between(x_mrna, m_y_new_upp, m_y_new_down, facecolor=plot_colors, alpha=0.5, interpolate=False)

    x_protein = np.arange(np.min(timepoints_num_protein), np.max(timepoints_num_protein), 0.01)
    spl = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[0, :], k=len(timepoints_num_mrna) - 1)
    spl_upp = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[1, :], k=len(timepoints_num_mrna) - 1)
    spl_low = interpolate.UnivariateSpline(timepoints_num_protein, protein_data[2, :], k=len(timepoints_num_mrna) - 1)
    p_y_new = spl(x_protein)
    p_y_new_upp = spl_upp(x_protein)
    p_y_new_down = spl_low(x_protein)
    plt.plot(x_protein, p_y_new, linestyle="--", label="Protein", color=plot_colors)
    plt.plot(x_protein, p_y_new_upp, linestyle="--", label="Protein", color=plot_colors)
    plt.plot(x_protein, p_y_new_down, linestyle="--", label="Protein", color=plot_colors)
    ax.fill_between(x_protein, p_y_new_upp, p_y_new_down, facecolor=plot_colors, alpha=0.25, interpolate=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_title(gene)
    plt.savefig(figpath, format='png')


#TODO never user => probably replaced by boxplot_MPI
def plot_dynamic_MPI(mrna_df, prot_df, genes, tgt_fp):
    plot_colors = constants.analysis_config['PLOT_COLORS']
    for i, gene in enumerate(genes):
        data_mrna = np.zeros((3, len(constants.dataset_config['TIMEPOINTS_MRNA'])))
        df_mrna = mrna_df[mrna_df["Gene"] == gene]
        cpt = 0
        for tp, line in df_mrna.groupby(['Timepoint']):
            mpi, err = helpers.compute_mpis(df_mrna[df_mrna["Timepoint"] == tp],
                                            constants.analysis_config['BOOTSTRAP_MPI'])
            # err_median = np.median(np.abs(np.tile(np.median(err), (1, len(err))) - err))
            upp_env = mpi + np.std(err)
            low_env = mpi - np.std(err)
            data_mrna[0, cpt] = mpi
            data_mrna[1, cpt] = upp_env
            data_mrna[2, cpt] = low_env
            cpt += 1
        data_prot = np.zeros((3, len(constants.dataset_config['TIMEPOINTS_PROTEIN'])))
        df_protein = prot_df[prot_df["Gene"] == gene]
        cpt = 0
        for tp, line in df_protein.groupby(['Timepoint']):
            mpi, err = helpers.compute_mpis_2(df_protein[df_protein["Timepoint"] == tp],
                                              constants.analysis_config['BOOTSTRAP_MPI'])
            upp_env = mpi + np.std(err)
            low_env = mpi - np.std(err)
            data_prot[0, cpt] = mpi
            data_prot[1, cpt] = upp_env
            data_prot[2, cpt] = low_env
            cpt += 1
        dynamic_profiles(data_mrna, data_prot, gene, 'Time(hrs)', 'MTOC polarity index', tgt_fp, plot_colors[i])


def sns_linear_regression(data_1, data_2, color, graph_file_path_name):
    sns.set(style="white", color_codes=True)
    annot_kws = {'prop': {'family': 'monospace', 'weight': 'bold', 'size': 8}}
    res1 = pearsonr(data_1, data_2)
    # print(", ".join(["x" + unichr(u) for u in (0x2070, 0x00B9, 0x00B2, 0x00B3, 0x2074, 0x2075, 0x2076, 0x2077, 0x2078, 0x2079)]))
    sns.set(font_scale=1)
    data3 = np.array(data_2) / np.array(data_1)
    data3 = data3[data3 <= 0.6]
    data3 = data3 - np.median(data3)
    sns_plot_regression = sns.jointplot(x=data_1, y=data_2, kind='reg', color=color)
    # sns_plot_regression.ax_marg_x.set_xlim(350,850)
    phantom, = sns_plot_regression.ax_joint.plot([], [], linestyle="", alpha=0)
    sns_plot_regression.ax_joint.legend([phantom], [
        'pearsonr={:f}, R' + chr(0x00B2) + '={:f}, p={:.2E}'.format(res1[0], res1[0] ** 2, res1[1])], **annot_kws)
    sns_plot_regression.savefig(graph_file_path_name, format="png")


def profile(profiles, genes, num_contours, figname):
    plot_colors = constants.analysis_config['PLOT_COLORS']
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    ax.yaxis.grid(which="major", color='black', linestyle='-', linewidth=0.25)
    ax.tick_params(right=False, top=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    plt.xticks(fontsize=30)
    plt.xticks([w for w in range(0, num_contours + 2, 10)])
    for i in range(len(genes)):
        plt.plot(np.arange(num_contours), profiles[i], color=plot_colors[i], linewidth=3, label=genes)
    plt.savefig(figname)
    plt.close()


def histogram_noise_measured(nm_arhgdia, nm_arhgdia_cultured, figname):
    plot_colors = constants.analysis_config['PLOT_COLORS']
    ax = plt.axes()
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=30)
    ind = np.arange(1, 3)
    width = 0.25
    # ax.set_xlim(-width * 2, len(ind) + width)
    ax.set_ylim(0, 0.8)
    ax.set_title('')
    xTickMarks = ["", ""]
    ax.set_xticks(ind)
    ax.bar(ind, [nm_arhgdia, nm_arhgdia_cultured], width, color=plot_colors)
    plt.savefig(figname)
    plt.close()


def plot_figure(total_mads_arhgdia, total_mads_arhgdia_cultured, figname):
    plt.figure()
    ax = plt.axes()
    ax.tick_params(right=False, top=False, bottom=False, direction='inout', length=8, width=3, colors='black')
    for axis in ['left']:
        ax.spines[axis].set_linewidth(3)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    plt.plot(np.mean(total_mads_arhgdia, axis=0), color='blue')
    plt.plot(np.mean(total_mads_arhgdia_cultured, axis=0), color='black')
    ax.set_xlim(0, len(total_mads_arhgdia[0]))
    plt.savefig(figname)
    plt.close()

def spline_graph(grid_mat, figname, band_n=100):
    ax1 = plt.subplot()
    x_mrna = np.arange(0, band_n, 0.5)
    fact = np.max(np.array(grid_mat)) / 10
    data = (np.array(grid_mat).flatten() / fact).astype(int) + 1
    spl = interpolate.UnivariateSpline(np.arange(0, band_n, 1), data)
    m_y_new = spl(x_mrna)
    plt.plot(x_mrna, m_y_new)
    ax1.set_ylim((0, 12))
    ax1.set_xlim((0, band_n - 1))
    plt.savefig(figname)
    plt.close()

def heatmap(grid_mat, figname, band_n=100):
    ax0 = plt.subplot()
    ax0.set_yticks([])
    major_ticks = np.arange(0, int(band_n) + 1, 1)
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.set_xticks(major_ticks)
    plt.imshow(grid_mat, cmap='coolwarm', aspect=band_n / 3)
    plt.ylim((0, 0.4))
    plt.savefig(figname)
    plt.close()


def compute_heatmap(ranking, gene, tgt_fp, size=4, xtickslabel=['2h', '3h', '5h', '7h'], ytickslabel = ['2h', '3h', '4h', '5h']):
    im = np.flipud(np.kron(ranking, np.ones((10, 10))))
    plt.imshow(im, extent=[0, size, 0, size], cmap='GnBu', interpolation='nearest')
    ax = plt.axes()
    ax.set_ylabel("mRNA  - Time (hrs)")
    ax.set_xlabel("Protein  - Time (hrs)")
    myxticklabels = xtickslabel
    ax.xaxis.set(ticks=np.arange(0.5, size + 0.5, 1), ticklabels=myxticklabels)
    myyticklabels = ytickslabel
    ax.yaxis.set(ticks=np.arange(0.5, size + 0.5, 1), ticklabels=myyticklabels)
    ax.set_title(gene)
    plt.savefig(tgt_fp)
    plt.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from dataclasses import dataclass
from typing import List

import pandas as pd
from loguru import logger

import constants
import mpi_calculator
import plot
from helpers import open_repo
from image_set import ImageSet
from path import global_root_dir
from mpi_calculator import DensityStats

# global config
OUTLIERS_THRESHOLD = 3


def build_labels(quadrants_num):
    return [f"Non MTOC{i}" for i in range(1, quadrants_num)]


def plot_boxplot_MPI(mrna_density_stats: DensityStats, protein_density_stats: DensityStats,
                     molecule_list, mrna_timepoint_list, protein_timepoint_list, tgt_fp):
    """
    The timepoint_list has to have the same order as in the mpi() function
    """
    plot_colors = constants.analysis_config['PLOT_COLORS']
    mrna_density_stats.update_group_key(['Timepoint'])
    protein_density_stats.update_group_key(['Timepoint'])

    for color_num, gene in enumerate(molecule_list):
        # Collect mrna and protein density statistics for this gene only
        dd = {'Molecule_type': [], 'Timepoint': [], 'MPI': [], 'err': []}
        gene_stats = mrna_density_stats.subset_stats("Gene", gene)
        prot_stats = protein_density_stats.subset_stats("Gene", gene)
        mrna_mpis, mrna_errs = gene_stats.mpi()
        prot_mpis, prot_errs = prot_stats.mpi()
        dd["MPI"] = mrna_mpis + prot_mpis
        dd["err"] = mrna_errs + prot_errs
        dd["Molecule_type"] = ["mrna"] * len(mrna_timepoint_list) + ["protein"] * len(protein_timepoint_list)
        dd["Timepoint"] = sorted(mrna_timepoint_list) + sorted(protein_timepoint_list)
        for tp in np.setdiff1d(mrna_timepoint_list, protein_timepoint_list):
            dd["Timepoint"].append(tp);
            dd["Molecule_type"].append("protein")
            dd["MPI"].append(0);
            dd["err"].append(0)
        for tp in np.setdiff1d(protein_timepoint_list, mrna_timepoint_list):
            dd["Timepoint"].append(tp);
            dd["Molecule_type"].append("mrna")
            dd["MPI"].append(0);
            dd["err"].append(0)

        df = pd.DataFrame(dd)
        df.sort_values("Timepoint", axis=0, ascending=True, inplace=True)

        my_pal = {"mrna": str(plot_colors[color_num]),
                  "protein": str(helpers.color_variant(plot_colors[color_num], +80))}
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_DYNAMIC_MPI'].format(gene=gene)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        helpers.create_dir_if_needed_for_filepath(tgt_fp)
        plot.sns_barplot(df, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type", err="err")
        # sns_barplot_simple(df, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type")
        logger.info("Generated image at {}", tgt_fp)


def compute_density_per_quadrant(repo, molecule_type, groupby, quadrants_num, quadrant_labels,
                                 molecule_list, time_points, mpi_sample_size):
    density_per_quadrant = []
    for molecule in molecule_list:
        for timepoint in time_points:
            image_set = ImageSet(repo, [f"{molecule_type}/{molecule}/{timepoint}/"])
            if image_set.__sizeof__() < 5:
                logger.warning("Image set is small for {}", molecule)
            dict_gene = image_set.compute_peripheral_normalised_quadrant_densities(
                quadrants_num=quadrants_num,
                mtoc_quadrant_label='MTOC',
                quadrant_labels=quadrant_labels
            )
            # dict_gene["MTOC leading edge"] = image_set.mtoc_is_in_leading_edge()
            dict_gene["Gene"] = molecule
            dict_gene["Timepoint"] = timepoint
            density_per_quadrant.append(pd.DataFrame(dict_gene))

    return DensityStats(
        df=pd.concat(density_per_quadrant),
        group_key=groupby,
        mpi_sample_size=mpi_sample_size,
        quadrant_labels=quadrant_labels,
        mtoc_quadrant_label='MTOC'
    )


configurations = [
    ["src/analysis/mtoc/config_periph_original.json", "", "", ""],
    ["src/analysis/mtoc/config_periph_prrc2c.json", "arhgdia", "prrc2c_depleted", "Timepoint"]
]

num_protein_quadrants = {
    "src/analysis/mtoc/config_periph_original.json": 3
}

# TODO use this info when building the DensityStats objects

group_keys = {
    "src/analysis/mtoc/config_periph_prrc2c.json": ['Gene', 'Timepoint']
}

# TODO recuperer commentaires indiquant nom de figure généré, le passer en argument des fonctions de plot
if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        ## mRNA
        genes_list = constants.dataset_config['MRNA_GENES']
        mrna_time_points = constants.dataset_config['TIMEPOINTS_MRNA']
        mpi_sample_size = constants.analysis_config['MPI_SUB_SAMPLE_SIZE']
        quadrant_labels = build_labels(quadrants_num=4)
        groupby = group_keys.get(conf[0], ['Gene'])
        mRNA_df = compute_density_per_quadrant(
            repo=repo,
            molecule_type="mrna",
            quadrants_num=4,
            quadrant_labels=quadrant_labels,
            molecule_list=genes_list,
            time_points=mrna_time_points,
            groupby=groupby,
            mpi_sample_size=mpi_sample_size
        )
        # plot.plot_hist_ratio(mRNA_df, 'mrna', limit_threshold=OUTLIERS_THRESHOLD)
        if conf[2] == "":
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT_RATIO'].format(
                molecule_type='mrna')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.compute_violin_plot_ratio(mRNA_df, 'mrna', tgt_fp)
            # # # # New Figure 3.C left mRNA MTOC enrichment for original data
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(
                molecule_type='mrna')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.compute_violin_plot_enrichment(mRNA_df, 'mrna', tgt_fp, limit_threshold=OUTLIERS_THRESHOLD)
        else:
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT_RATIO'].format(molecule_type='mrna')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),tgt_image_name)
            plot.compute_categorical_violin_plot_ratio(mRNA_df, 'mrna', tgt_fp, limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], gene=conf[1], groupby=[conf[3]])

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type='mrna')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),tgt_image_name)
            plot.compute_categorical_violin_plot_enrichment(mRNA_df, 'mrna', tgt_fp, limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], groupby=conf[3])

        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MPI'].format(molecule_type='mrna')
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        plot.plot_MPI(mRNA_df, 'mrna', tgt_fp)

        ## proteins
        protein_time_points = constants.dataset_config['TIMEPOINTS_PROTEIN']
        protein_list = constants.dataset_config['PROTEINS']
        protein_quadrant_labels = build_labels(num_protein_quadrants.get(conf[0], 4))
        prot_df = compute_density_per_quadrant(
            repo=repo,
            molecule_type="protein",
            quadrants_num=num_protein_quadrants.get(conf[0], 4),
            quadrant_labels=protein_quadrant_labels,
            molecule_list=protein_list,
            time_points=protein_time_points,
            groupby=groupby,
            mpi_sample_size=mpi_sample_size
        )
        if conf[2] == "":
            # plot.compute_violin_plot_ratio(prot_df, 'protein')
            # # # # New Figure 3.C left mRNA MTOC enrichment for original data
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(
                molecule_type='protein')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.compute_violin_plot_enrichment(prot_df, 'protein', tgt_fp, limit_threshold=OUTLIERS_THRESHOLD)
        else:
            # plot.compute_categorical_violin_plot_ratio(prot_df, 'protein', limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], gene=conf[1], groupby=[conf[3]])
            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type='protein')
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
            plot.compute_categorical_violin_plot_enrichment(prot_df, 'protein', tgt_fp, limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], groupby=conf[3])

        # plot.plot_mtoc_enrichment(prot_df, 'protein', limit_threshold=5)
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MPI'].format(molecule_type='protein')
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        plot.plot_MPI(prot_df, 'protein', tgt_fp)

        ## combined
        if ("original" in conf[0]):
            plot.plot_boxplot_MPI(mRNA_df, prot_df, protein_list, mrna_time_points, protein_time_points)

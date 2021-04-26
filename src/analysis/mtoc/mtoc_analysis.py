#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from dataclasses import dataclass
from typing import List

import pandas as pd
from loguru import logger
import numpy as np
import constants
import mpi_calculator
import plot
from helpers import open_repo, create_dir_if_needed_for_filepath, color_variant
from image_set import ImageSet
from path import global_root_dir
from mpi_calculator import DensityStats

# global config
OUTLIERS_THRESHOLD = 3


def build_labels(quadrants_num):
    return [f"Non MTOC{i}" for i in range(1, quadrants_num)]


def plot_boxplot_MPI(mrna_density_stats: DensityStats, protein_density_stats: DensityStats,
                     molecule_list, mrna_timepoint_list, protein_timepoint_list):
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
                  "protein": str(color_variant(plot_colors[color_num], +80))}
        tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_DYNAMIC_MPI'].format(gene=gene)
        tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                              tgt_image_name)
        create_dir_if_needed_for_filepath(tgt_fp)
        plot.sns_barplot(df, my_pal, tgt_fp, x="Timepoint", y="MPI", hue="Molecule_type", err="err")
        logger.info("Generated image at {}", str(tgt_fp).split("analysis/")[1])


def compute_density_per_quadrant(analysis_repo, molecule_type, groupby, quadrants_num, quadrant_labels,
                                 molecule_list, time_points, mpi_sample_size):
    density_per_quadrant = []
    for molecule in molecule_list:
        for timepoint in time_points:
            image_set = ImageSet(analysis_repo, [f"{molecule_type}/{molecule}/{timepoint}/"])
            if image_set.__sizeof__() < 5:
                logger.warning("Image set is small for {}", molecule)
            dict_gene = image_set.compute_normalised_quadrant_densities(
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
    ["src/analysis/mtoc/config_original.json", "", "", ""],
    ["src/analysis/mtoc/config_nocodazole_arhgdia.json", "arhgdia", "Nocodazole+", "Gene"],
    ["src/analysis/mtoc/config_nocodazole_pard3.json", "pard3", "Nocodazole+", "Gene"],
    ["src/analysis/mtoc/config_prrc2c.json", "arhgdia", "prrc2c_depleted", "Timepoint"]
]

num_protein_quadrants = {
    "src/analysis/mtoc/config_original.json": 3,
}

# TODO use this info when building the DensityStats objects

group_keys = {
    "src/analysis/mtoc/config_prrc2c.json": ['Gene', 'Timepoint']
}

# figure 3B left/right mRNA/protein MTOC cytoplasmic enrichment for original data (violin plot)
# figure 3C mRNA/protein cytoplasmic MPI
# figure 4C top/bottom mRNA/protein MTOC cytoplasmic enrichment for prrc2c (violin plot)
# figure 3D (2 figures) and S3D (dynamic cytoplasmic MPI) (2 figures)


if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = open_repo()
        tp_mrna = constants.dataset_config['TIMEPOINTS_MRNA']
        tp_proteins = constants.dataset_config['TIMEPOINTS_PROTEIN']
        _mpi_sample_size = constants.analysis_config['MPI_SUB_SAMPLE_SIZE']

        _groupby = group_keys.get(conf[0], ['Gene'])
        dfs = []
        for genes, timepoints, _molecule_type, quads, in zip([constants.analysis_config['MRNA_GENES'], constants.analysis_config['PROTEINS']], [tp_mrna, tp_proteins],
                                                            ["mrna", "protein"], [4, num_protein_quadrants.get(conf[0], 4)]):
            _quadrant_labels = build_labels(quadrants_num=quads)
            df = compute_density_per_quadrant(
                analysis_repo=repo,
                molecule_type=_molecule_type,
                quadrants_num=quads,
                quadrant_labels=_quadrant_labels,
                molecule_list=genes,
                time_points=timepoints,
                groupby=_groupby,
                mpi_sample_size=_mpi_sample_size
            )
            dfs.append(df)

            if conf[2] == "":
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT_RATIO'].format(molecule_type=_molecule_type)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.compute_violin_plot_ratio(df, _molecule_type, tgt_fp)

                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type=_molecule_type)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.compute_violin_plot_enrichment(df, _molecule_type, tgt_fp, limit_threshold=OUTLIERS_THRESHOLD)
            else:
                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT_RATIO'].format(molecule_type=_molecule_type)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.compute_categorical_violin_plot_ratio(df, _molecule_type, tgt_fp, limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], gene=conf[1], groupby=[conf[3]])

                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type=_molecule_type)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.compute_categorical_violin_plot_enrichment(df, _molecule_type, tgt_fp, limit_threshold=OUTLIERS_THRESHOLD, term=conf[2], groupby=conf[3])

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MPI'].format(molecule_type=_molecule_type)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.plot_MPI(df, _molecule_type, tgt_fp)

        if "original" in conf[0]:
            plot_boxplot_MPI(dfs[0], dfs[1], constants.analysis_config['PROTEINS'], tp_mrna, tp_proteins)

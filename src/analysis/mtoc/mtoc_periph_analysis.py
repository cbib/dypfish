#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
import pandas as pd
from loguru import logger
import constants
import plot
import helpers
from image_set import ImageSet
from path import global_root_dir
from mpi_calculator import DensityStats

# global config
OUTLIERS_THRESHOLD = 3


def build_labels(quadrants_num):
    return [f"Non MTOC{i}" for i in range(0, quadrants_num-1)]

def compute_density_per_quadrant(analysis_repo, molecule_type, groupby_key, quadrants_num, quadrant_labels,
                                 molecule_list, time_points, mpi_sample_size):
    assert len(quadrant_labels) == quadrants_num-1
    density_per_quadrant = []
    dict_gene = {}
    for molecule in molecule_list:
        for timepoint in time_points:
            image_set = ImageSet(analysis_repo, [f"{molecule_type}/{molecule}/{timepoint}/"])
            if image_set.__sizeof__() < 5:
                logger.warning("Image set is small for {}", molecule)
            res = image_set.compute_normalised_quadrant_densities(quadrants_num=quadrants_num, peripheral_flag=False)
            mtoc_quadrants = res[res[:,1]==1][:, 0]
            non_mtoc_quadrants = res[res[:,1]==0][:, 0].reshape(image_set.__sizeof__(), quadrants_num-1)
            dict_gene["Gene"] = [molecule for i in range(image_set.__sizeof__())]
            dict_gene["Timepoint"] = [timepoint for i in range(image_set.__sizeof__())]
            dict_gene["MTOC"] = mtoc_quadrants
            for i in range(0, quadrants_num-1):
                dict_gene["Non MTOC" + str(i)] = non_mtoc_quadrants[:,i]
            density_per_quadrant.append(pd.DataFrame(dict_gene))

    return DensityStats(df=pd.concat(density_per_quadrant), group_key=groupby_key, mpi_sample_size=mpi_sample_size,
                        quadrant_labels=quadrant_labels, mtoc_quadrant_label='MTOC')


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

# figure S3A left/right mRNA/protein MTOC peripheral enrichment for original data (violin plot)
# figure S3B mRNA/protein peripheral MPI
# figure S3C dynamic peripheral MPI (4 figures)

if __name__ == '__main__':
    for conf in configurations:
        conf_full_path = pathlib.Path(global_root_dir, conf[0])
        constants.init_config(analysis_config_js_path=conf_full_path)
        repo = helpers.open_repo()
        tp_mrna = constants.dataset_config['TIMEPOINTS_MRNA']
        tp_proteins = constants.dataset_config['TIMEPOINTS_PROTEIN']
        mrna_genes = constants.analysis_config['MRNA_GENES']
        proteins = constants.analysis_config['PROTEINS']
        mpi_sample_size = constants.analysis_config['MPI_SUB_SAMPLE_SIZE']

        groupby_key = group_keys.get(conf[0], ['Gene'])
        dfs = []
        for genes, timepoints, _molecule_type, quads, in zip([mrna_genes, proteins], [tp_mrna, tp_proteins],
                                                             ["mrna", "protein"],
                                                             [4, num_protein_quadrants.get(conf[0], 4)]):
            quadrant_labels = build_labels(quadrants_num=quads)
            df = compute_density_per_quadrant(analysis_repo=repo, molecule_type=_molecule_type,
                                              quadrants_num=quads, quadrant_labels=quadrant_labels,
                                              molecule_list=genes, time_points=timepoints,
                                              groupby_key=groupby_key, mpi_sample_size=mpi_sample_size)
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
                plot.compute_categorical_violin_plot_ratio(df, _molecule_type, tgt_fp, limit_threshold=OUTLIERS_THRESHOLD,
                                                           term=conf[2], gene=conf[1], groupby_key=[conf[3]])

                tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type=_molecule_type)
                tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
                plot.compute_categorical_violin_plot_enrichment(df, _molecule_type, tgt_fp, limit_threshold=OUTLIERS_THRESHOLD,
                                                                term=conf[2], groupby_key=conf[3])

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MPI'].format(molecule_type=_molecule_type)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.plot_MPI(df, _molecule_type, tgt_fp)

        if "original" in conf[0]:
            plot.plot_boxplot_MPI(dfs[0], dfs[1], constants.analysis_config['PROTEINS'], tp_mrna, tp_proteins)

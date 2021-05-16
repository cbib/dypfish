#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib

import pandas as pd
from loguru import logger

import constants
import helpers
import plot
from image_set import ImageSet
from mpi_calculator import DensityStats
from path import global_root_dir

# global config
OUTLIERS_THRESHOLD = 3


def build_labels(quadrants_num):
    return [f"Non MTOC{i}" for i in range(0, quadrants_num-1)]

def compute_density_per_quadrant(analysis_repo, molecule_type, groupby_key, quadrants_num, quadrant_labels,
                                 molecule_list, time_points, mpi_sample_size, peripheral_flag=False):
    assert len(quadrant_labels) == quadrants_num - 1
    density_per_quadrant = []
    for molecule in molecule_list:
        for timepoint in time_points:
            density_statictics = {}
            image_set = ImageSet(analysis_repo, [f"{molecule_type}/{molecule}/{timepoint}/"])
            if len(image_set.images) < 5:
                logger.warning("Image set is small for {}", molecule)
            res = image_set.compute_normalised_quadrant_densities(quadrants_num=quadrants_num,
                                                                  peripheral_flag=peripheral_flag)
            mtoc_quadrants = res[res[:,1] == 1][:, 0]
            num_images = res.shape[0] // quadrants_num
            non_mtoc_quadrants = res[res[:,1] == 0][:, 0].reshape(num_images, quadrants_num-1)
            density_statictics["Gene"] = [molecule] * num_images
            density_statictics["Timepoint"] = [timepoint] * num_images
            density_statictics["MTOC"] = mtoc_quadrants
            for i in range(0, quadrants_num-1):
                density_statictics["Non MTOC" + str(i)] = non_mtoc_quadrants[:,i]
            density_per_quadrant.append(pd.DataFrame(density_statictics))

    return DensityStats(df=pd.concat(density_per_quadrant), group_key=groupby_key, mpi_sample_size=mpi_sample_size,
                        quadrant_labels=quadrant_labels, mtoc_quadrant_label='MTOC')

configurations = [
    ["src/analysis/mtoc/config_original.json", "Gene", [4, 4]], # it was 3 for proteins
    ["src/analysis/mtoc/config_nocodazole_arhgdia.json", "Gene", [4, 4]],
    ["src/analysis/mtoc/config_nocodazole_pard3.json", "Gene", [4, 4]],
    ["src/analysis/mtoc/config_prrc2c.json", "Timepoint", [4, 4]]
    ["src/analysis/mtoc/config_periph_original.json", "Gene", [4, 4]],
    ["src/analysis/mtoc/config_periph_prrc2c.json", "Timepoint", [4, 4]]
]

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
        repo = helpers.open_repo()
        tp_mrna = constants.dataset_config['TIMEPOINTS_MRNA']
        tp_proteins = constants.dataset_config['TIMEPOINTS_PROTEIN']
        mrna_genes = constants.analysis_config['MRNA_GENES']
        proteins = constants.analysis_config['PROTEINS']
        mpi_sample_size = constants.analysis_config['MPI_SUB_SAMPLE_SIZE']
        peripheral_flag = "periph" in conf[0]

        group_key = group_keys.get(conf[0], ['Gene'])
        dfs = []
        for genes, timepoints, molecule_type, quads, in zip([mrna_genes, proteins], [tp_mrna, tp_proteins],
                                                            ["mrna", "protein"], conf[2]):
            logger.info("Running MPI analysis for {}", conf[0])
            quadrant_labels = build_labels(quadrants_num=quads)
            density_stats = compute_density_per_quadrant(analysis_repo=repo, molecule_type=molecule_type,
                                                             quadrants_num=quads, quadrant_labels=quadrant_labels,
                                                             molecule_list=genes, time_points=timepoints,
                                                             groupby_key=group_key, mpi_sample_size=mpi_sample_size,
                                                             peripheral_flag=peripheral_flag)
            dfs.append(density_stats)

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MTOC_ENRICHMENT'].format(molecule_type=molecule_type)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
            plot.enrichment_violin_plot(density_stats, molecule_type, tgt_fp,
                                        groupby_key=conf[1], limit_threshold=OUTLIERS_THRESHOLD)
            logger.info("Created figure {}", tgt_fp)
            plot.plot_clusters(molecule_type, density_stats.df, peripheral_flag=peripheral_flag)

            tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_MPI'].format(molecule_type=molecule_type)
            tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir),
                                  tgt_image_name)
            plot.plot_MPI(density_stats, molecule_type, tgt_fp)
            logger.info("Created figure {}", tgt_fp)

        if "original" in conf[0]:
            plot.plot_boxplot_MPI(dfs[0], dfs[1], constants.analysis_config['PROTEINS'], tp_mrna, tp_proteins)

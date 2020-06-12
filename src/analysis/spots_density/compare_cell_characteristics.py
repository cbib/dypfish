#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import pathlib
from loguru import logger
import constants
import pandas as pd
from plot import sns_boxplot, sns_linear_regression
from image_set import ImageSet
from path import global_root_dir
from repository import H5RepositoryWithCheckpoint

constants.init_config(
    analysis_config_js_path=pathlib.Path(global_root_dir, "src/analysis/spots_density/config_original.json"))
dataset_root_fp = pathlib.Path(constants.analysis_config['DATASET_CONFIG_PATH'].format(root_dir=global_root_dir)).parent
primary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['PRIMARY_FILE_NAME'])
secondary_fp = pathlib.Path(dataset_root_fp, constants.dataset_config['SECONDARY_FILE_NAME'])
analysis_repo = H5RepositoryWithCheckpoint(repo_path=primary_fp, secondary_repo_path=secondary_fp)

# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
###                                 Figure 1.D bottom left (cell area) Standard vs Micropatterned
# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
logger.info("cell area comparison for arhgdia et arhgdia_cultured 3h data")
df_cell_area = pd.DataFrame(columns=["Gene", "value"])
dict_Cell_area = {"Gene": [], "value": []}
[dict_Cell_area["value"].append(image.compute_cell_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/3h/"]).get_images()]
[dict_Cell_area["Gene"].append("Micropatterned cells") for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/3h/"]).get_images()]
df_cell_area = pd.concat([df_cell_area, pd.DataFrame(dict_Cell_area)])

dict_Cell_area = {"Gene": [], "value": []}
[dict_Cell_area["value"].append(image.compute_cell_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]
[dict_Cell_area["Gene"].append("Standard cells") for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]
df_cell_area = pd.concat([df_cell_area, pd.DataFrame(dict_Cell_area)])

my_pal = {"Micropatterned cells": "#66b2ff", "Standard cells": "#1a8cff"}
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BOXPLOT'].format(model="cell_area")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
sns_boxplot(df_cell_area, my_pal, tgt_fp)

# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
###                                 Figure 1.D bottom left (nucleus area) Standard vs Micropatterned
# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
logger.info("nucleus area comparison for arhgdia et arhgdia_cultured 3h data")
df_nucleus_area = pd.DataFrame(columns=["Gene", "value"])
dict_nucleus_area = {"Gene": [], "value": []}
[dict_nucleus_area["value"].append(image.compute_nucleus_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/3h/"]).get_images()]
[dict_nucleus_area["Gene"].append("Micropatterned cells") for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/3h/"]).get_images()]
df_nucleus_area = pd.concat([df_nucleus_area, pd.DataFrame(dict_nucleus_area)])

dict_nucleus_area = {"Gene": [], "value": []}
[dict_nucleus_area["value"].append(image.compute_nucleus_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]
[dict_nucleus_area["Gene"].append("Standard cells") for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]
df_nucleus_area = pd.concat([df_nucleus_area, pd.DataFrame(dict_nucleus_area)])

my_pal = {"Micropatterned cells": "#66b2ff", "Standard cells": "#1a8cff"}
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_BOXPLOT'].format(model="nucleus_area")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
sns_boxplot(df_nucleus_area, my_pal, tgt_fp)

# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
###                                 Figure 1.D Top right(transcript total count by cell area)Micropatterned cell
# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
dict_transcript_by_cell = {"total_transcript": [], "cell_area": []}
[dict_transcript_by_cell["total_transcript"].append(image.compute_cytoplasmic_total_spots()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/"]).get_images()]
[dict_transcript_by_cell["cell_area"].append(image.compute_cell_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia/"]).get_images()]
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT'].format(cell_type="micropatterned")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
sns_linear_regression(dict_transcript_by_cell["cell_area"], dict_transcript_by_cell["total_transcript"], "#0A3950",
                      tgt_fp)

# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
###                                 Figure 1.D Top left (transcript total count by cell area) Cultured cells
# ################################## ################################## ################################## #################################
# ################################## ################################## ################################## #################################
dict_transcript_by_cell = {"total_transcript": [], "cell_area": []}
[dict_transcript_by_cell["total_transcript"].append(image.compute_cytoplasmic_total_spots()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]
[dict_transcript_by_cell["cell_area"].append(image.compute_cell_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/3h/"]).get_images()]

[dict_transcript_by_cell["total_transcript"].append(image.compute_cytoplasmic_total_spots()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/1h/"]).get_images()]
[dict_transcript_by_cell["cell_area"].append(image.compute_cell_area()) for image in
 ImageSet(analysis_repo, ["mrna/arhgdia_cultured/1h/"]).get_images()]
tgt_image_name = constants.analysis_config['FIGURE_NAME_FORMAT_PLOT'].format(cell_type="standard")
tgt_fp = pathlib.Path(constants.analysis_config['FIGURE_OUTPUT_PATH'].format(root_dir=global_root_dir), tgt_image_name)
sns_linear_regression(dict_transcript_by_cell["cell_area"], dict_transcript_by_cell["total_transcript"], "#1E95BB",
                      tgt_fp)

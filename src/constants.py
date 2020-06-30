#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import json
import path
from loguru import logger

analysis_config = None
dataset_config = None

# basic image descriptors
ZLINES_PATH_SUFFIX = "/z_lines"
Z_LINE_DISTANCE_PATH_SUFFIX= "/z_line_distance"
SPOTS_PATH_SUFFIX = "/spots"
IF_PATH_SUFFIX = "/IF"
CELL_MASK_PATH_SUFFIX = "/cell_mask"
NUCLEUS_MASK_PATH_SUFFIX = "/nucleus_mask"
CYTOPLASM_MASK_PATH_SUFFIX = "/cytoplasm_mask"
PERIPHERAL_MASK_PATH_SUFFIX = "/peripheral_mask"
NUCLEUS_CENTROID_PATH_SUFFIX = "/nucleus_centroid"
MTOC_POSITION_PATH_SUFFIX = '/mtoc_position'
MTOC_LEADING_EDGE_SUFFIX = '/mtoc_is_in_leading_edge'

# secondary basic image descriptors
NUCLEUS_AREA_PATH_SUFFIX = "/nucleus_area"
CELL_AREA_PATH_SUFFIX = "/cell_area"
CELL_MASK_DISTANCE_PATH_SUFFIX = '/cell_mask_distance_map'
SPOTS_PERIPHERAL_DISTANCE_2D_PATH_SUFFIX = '/spots_peripheral_distance_2D'
CYTOPLASMIC_SPOTS_PATH_SUFFIX = '/cytoplasmic_spots'
PERIPHERAL_SPOTS_PATH_SUFFIX = '/peripheral_spots'
CYTOPLASMIC_SPOT_COUNT_PATH_SUFFIX = '/cytoplasmic_spot_count'
PERIPHERAL_SPOT_COUNT_PATH_SUFFIX = '/peripheral_spot_count'
TOTAL_INTENSITY_PATH_SUFFIX = '/total_intensity'
CELL_TOTAL_INTENSITY_PATH_SUFFIX = '/cell_total_intensity'
CYTOPLASMIC_TOTAL_INTENSITY_PATH_SUFFIX = '/cytoplasmic_total_intensity'
PERIPHERAL_TOTAL_INTENSITY_PATH_SUFFIX = '/peripheral_total_intensity'
CYTOPLASMIC_INTENSITIES_PATH_SUFFIX = '/cytoplasmic_intensities'
PERIPHERAL_INTENSITIES_PATH_SUFFIX = '/peripheral_intensities'
QUADRANT_MASK_PATH_SUFFIX = '/quadrant_mask'

# 3d image descriptors
HEIGHT_MAP_PATH_SUFFIX = '/height_map'
ZERO_LEVEL_PATH_SUFFIX = '/zero_level'
NUCLEUS_VOLUME_PATH_SUFFIX = "/nucleus_volume"
CELL_MASK_SLICES_PATH_SUFFIX = '/cell_mask_slices'
CELL_VOLUME_PATH_SUFFIX = "/cell_volume"
SPOTS_PERIPHERAL_DISTANCE_3D_PATH_SUFFIX = '/spots_peripheral_distance_3D'
CLUSTERING_INDICES_PATH_SUFFIX = '/clustering_index'
QUADRANT_DENSITIES_PATH_SUFFIX = '/quadrant_densities'
PERIPHERAL_QUADRANT_DENSITIES_PATH_SUFFIX = '/peripheral_quadrant_densities'
QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX = '/quadrant_and_slice_densities'
PERIPHERAL_QUADRANT_AND_SLICE_DENSITIES_PATH_SUFFIX = '/peripheral_quadrant_and_slice_densities'

def init_config(analysis_config_js_path):
    global analysis_config
    global dataset_config
    with open(analysis_config_js_path, "r") as fh:
        analysis_config = json.load(fh)

    dataset_config_file = analysis_config['DATASET_CONFIG_PATH'].format(root_dir=path.global_root_dir)
    with open(dataset_config_file, "r") as fh:
        dataset_config = json.load(fh)
    logger.info("Configs {} and {} loaded", analysis_config_js_path, dataset_config_file)


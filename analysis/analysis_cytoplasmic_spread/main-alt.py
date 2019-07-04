#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
import h5py
import src.plot as plot
from src import acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors
from src.utils import cell_type_micropatterned


CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX = "cytoplasmic_spread.png"

def cytoplasmic_spread(
    basic_h5_file_handler,
    genes,
    cell_type,
    molecule_type,
    save_into_dir,
    logger=None
    ):  # TODO argument genes vs proteins ?
    """
    * Generates a cytoplasmic spread graph and saves it in a file
    * Returns a tuple embedding 2 values : temporary file path name,
      associated metadata (as a dict)
    * Parameters :
    - basic_h5_file_handler : handler linked to a 'basic.h5' file
    - genes : list of strings (genes IDs)
      Ex. : ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    - cell_type : cell type ID
      Ex. : "micropatterned"
    - molecule_type : a string (molecule type ID)
      Ex. : "mrna" or "protein"
    - save_into_dir : a string (path of the directory where the graph will be saved)
    - logger : external logger
    """

    if logger:
        logger.debug("Generating cytoplasmic spread graph (%s)..." % molecule_type)

    graph_filename = "%s_%s" % (molecule_type, CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX)
    graph_file_path = os.path.join(save_into_dir, graph_filename)
    graph_metadata = {"filename": graph_filename}

    # for each gene in genes list, compute cytoplasmic spread
    cyt_spreads = []

    for gene in genes:
        image_list = helps.preprocess_image_list2(
            basic_h5_file_handler,
            '/' + molecule_type, gene
            )
        cyt_spreads.append(
            adsc.compute_cytoplasmic_spread(
                image_list,
                basic_h5_file_handler,
                path_data=None
                )
            )

    plot.bar_profile(cyt_spreads, genes, graph_file_path)

    if logger:
        logger.debug("Graph generated on path : %s" % graph_file_path)

    return {graph_file_path: graph_metadata}



CYTOPLASMIC_SPREAD_GRAPHS_DYNMIC_PROFILES_PREFIX = "cyt_spread"

def cytoplasmic_spread_dynamic_profiles(
    basic_h5_file_handler,
    genes,
    proteins,
    save_into_dir,
    logger=None
    ):
    """
    * Generates cytoplasmic spread dynamic profiles and saves them in files
    * Returns a list of tuples embedding 2 values :
      temporary file path name, associated metadata (as a dict)
    * Parameters :
    - basic_h5_file_handler : handler linked to a 'basic.h5' file
    - genes : list of strings (genes IDs)
      Ex. : ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    - proteins :  list of strings (proteins IDs)
      Ex. : ["beta_actin", "arhgdia", "gapdh", "pard3"]
    - save_into_dir : a string (path of the directory where the graphs will be saved)
    """

    # This part produces plot interpolation of cytoplasmic spread by timepoint
    if logger:
        logger.debug("Generating cytoplasmic spread dynamic profiles...")

    graphs_details = []

    if logger:
        logger.debug("Extracting data...")

    data_generator = plot.data_extractor(
        genes,
        proteins,
        basic_h5_file_handler,
        adsc.compute_cytoplasmic_spread,
        basic_h5_file_handler,
        path.path_data
        )

    try:
        for mrna_data, protein_data, i in data_generator:
            gene = genes[i]

            if logger:
                logger.debug("Processing gene : %s..." % gene)

            try:
                graph_filename = "%s_%s.png" % (CYTOPLASMIC_SPREAD_GRAPHS_DYNMIC_PROFILES_PREFIX, gene)

                if logger:
                    logger.debug("Generating graph %s..." % graph_filename)

                graph_file_path = os.path.join(save_into_dir, graph_filename)
                graph_metadata = {"filename": graph_filename}
                plot.dynamic_profiles(
                    mrna_data,
                    protein_data,
                    genes[i],
                    plot_colors[i],
                    'Time(hrs)',
                    'Cytoplasmic spread',
                    graph_file_path
                    )

                if logger:
                    logger.debug("Graph generated on path : %s" % graph_file_path)

                graphs_details.append({graph_file_path: graph_metadata})

            except Exception as e:
                if logger:
                    logger.exception(
                        "Could not generate cytoplasmic spread dynamic profile for gene : %s" % gene
                        )

    except Exception as e:
        if logger:
            logger.exception(
                "Got exception ! (maybe raised by data_generator ?)"
                )

    return graphs_details


def main(
    basic_h5_file_path_name,
    save_into_dir
    ):

    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()

    # Compute bar plot cytoplasmic spread
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(path.basic_file_path, "r") as basic_h5_file_handler:
        cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            cell_type=cell_type_micropatterned,
            molecule_type='mrna',
            save_into_dir=save_into_dir
            )
        cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            cell_type=cell_type_micropatterned,
            molecule_type='protein',
            save_into_dir=save_into_dir
            )
        cytoplasmic_spread_dynamic_profiles(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            save_into_dir=save_into_dir
            )



if __name__ == "__main__":
    basic_h5_file_path_name = path.basic_file_path
    save_into_dir = path.analysis_dir + "/analysis_cytoplasmic_spread/figures/"
    main(basic_h5_file_path_name, save_into_dir)

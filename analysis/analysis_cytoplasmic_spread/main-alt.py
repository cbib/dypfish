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


CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX = "cytoplasmic_spread"

def cytoplasmic_spread(
    basic_h5_file_handler,
    genes,
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
    - molecule_type : a string (molecule type ID)
      Ex. : "mrna" or "protein"
    - save_into_dir : a string (path of the directory where the graph will be saved)
    - logger : external logger
    """
    assert os.path.isdir(save_into_dir)

    if logger:
        logger.debug("Generating cytoplasmic spread graph (%s)..." % molecule_type)

    graph_file_name = "%s_%s.png" % (molecule_type, CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX)
    graph_file_path_name = os.path.join(save_into_dir, graph_file_name)
    graph_metadata = {"filename": graph_file_name}

    # for each gene in genes list, compute cytoplasmic spread
    cyt_spreads = []

    for gene in genes:
        image_list = helps.preprocess_image_list2(
            basic_h5_file_handler,
            '/' + molecule_type,
            gene
            )
        cyt_spreads.append(
            adsc.compute_cytoplasmic_spread(
                image_list,
                basic_h5_file_handler,
                path_data=None  # should be the path to raw images data, but in fact this parameter is not used in this case
                )
            )

    plot.bar_profile(cyt_spreads, genes, graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))

    if logger:
        logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}



CYTOPLASMIC_SPREAD_GRAPHS_DYNAMIC_PROFILES_PREFIX = "cyt_spread"

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
    assert os.path.isdir(save_into_dir)

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
        None   # should be the path to raw images data, but in fact this parameter is not used in this case
        )

    try:
        for mrna_data, protein_data, i in data_generator:
            gene = genes[i]
            colour = plot_colors[i]

            if logger:
                logger.debug("Processing gene : %s..." % gene)

            try:
                if logger:
                    logger.debug("Generating graph %s..." % graph_file_name)

                graph_file_name = "%s_%s.png" % (CYTOPLASMIC_SPREAD_GRAPHS_DYNAMIC_PROFILES_PREFIX, gene)
                graph_file_path_name = os.path.join(save_into_dir, graph_file_name)
                graph_metadata = {"filename": graph_file_name}

                plot.dynamic_profiles(
                    mrna_data,
                    protein_data,
                    gene,
                    colour,
                    'Time(hrs)',
                    'Cytoplasmic spread',
                    graph_file_path_name
                    )

                assert(os.path.isfile(graph_file_path_name))

                if logger:
                    logger.debug("Graph generated on path : %s" % graph_file_path_name)

                graphs_details.append({graph_file_path_name: graph_metadata})

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

    resulting_graphs_details = []

    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()

    # Compute bar plot cytoplasmic spread
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(basic_h5_file_path_name, "r") as basic_h5_file_handler:

        graph_details = cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            molecule_type='mrna',
            save_into_dir=save_into_dir
            )
        resulting_graphs_details.append(graph_details)

        graph_details = cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            molecule_type='protein',
            save_into_dir=save_into_dir
            )
        resulting_graphs_details.append(graph_details)

        graphs_details = cytoplasmic_spread_dynamic_profiles(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            save_into_dir=save_into_dir
            )
        resulting_graphs_details += graphs_details
        print(resulting_graphs_details)


if __name__ == "__main__":
    basic_h5_file_path_name = path.basic_file_path
    save_into_dir = os.path.join(path.analysis_dir, "analysis_cytoplasmic_spread/figures/")
    if not os.path.isdir(save_into_dir):
        os.mkdir(save_into_dir)
    main(basic_h5_file_path_name, save_into_dir)

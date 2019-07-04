#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
from collections import OrderedDict
import h5py
import src.plot as plot
from src import acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors


PNG_IMAGES_MIME_TYPE = "image/png"

CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX = "cytoplasmic_spread"

def cytoplasmic_spread(
    raw_images_dir_path_name,
    basic_h5_file_handler,
    genes,
    molecule_type,
    save_into_dir_path_name,
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
    - save_into_dir_path_name : a string (path of the directory where the graph will be saved)
    - logger : external logger
    """
    assert os.path.isdir(save_into_dir_path_name)

    if logger:
        logger.debug("Generating cytoplasmic spread graph (%s)..." % molecule_type)

    graph_file_name = "%s_%s.png" % (molecule_type, CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX)
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
        }

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
                raw_images_dir_path_name
                )
            )

    plot.bar_profile(cyt_spreads, genes, graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))

    if logger:
        logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}



CYTOPLASMIC_SPREAD_DYNAMIC_PROFILES_GRAPHS_PREFIX = "cyt_spread"

def cytoplasmic_spread_dynamic_profiles(
    raw_images_dir_path_name,
    basic_h5_file_handler,
    genes,
    proteins,
    save_into_dir_path_name,
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
    - save_into_dir_path_name : a string (path of the directory where the graphs will be saved)
    """
    assert os.path.isdir(save_into_dir_path_name)

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
        raw_images_dir_path_name
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

                graph_file_name = "%s_%s.png" % (CYTOPLASMIC_SPREAD_DYNAMIC_PROFILES_GRAPHS_PREFIX, gene)
                graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
                graph_metadata = {
                    "file_name": graph_file_name,
                    "mime_type": PNG_IMAGES_MIME_TYPE
                    }

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
    raw_images_dir_path_name,
    basic_h5_file_path_name,
    save_into_dir_path_name
    ):
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    enable_logger()

    resulting_graphs_details_as_list = []

    # Compute bar plot cytoplasmic spread
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(basic_h5_file_path_name, "r") as basic_h5_file_handler:

        graph_details = cytoplasmic_spread(
            raw_images_dir_path_name=raw_images_dir_path_name,
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            molecule_type='mrna',
            save_into_dir_path_name=save_into_dir_path_name
            )
        resulting_graphs_details_as_list.append(graph_details)

        graph_details = cytoplasmic_spread(
            raw_images_dir_path_name=raw_images_dir_path_name,
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            molecule_type='protein',
            save_into_dir_path_name=save_into_dir_path_name
            )
        resulting_graphs_details_as_list.append(graph_details)

        graphs_details = cytoplasmic_spread_dynamic_profiles(
            raw_images_dir_path_name=raw_images_dir_path_name,
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            save_into_dir_path_name=save_into_dir_path_name
            )
        resulting_graphs_details_as_list += graphs_details

    resulting_graphs_details_as_dict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_dict.update(graph_details)

    return resulting_graphs_details_as_dict


if __name__ == "__main__":

    raw_images_dir_path_name = path.path_data
    basic_h5_file_path_name = path.basic_file_path
    save_into_dir_path_name = os.path.join(path.analysis_dir, "analysis_cytoplasmic_spread/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details = main(
        raw_images_dir_path_name,
        basic_h5_file_path_name,
        save_into_dir_path_name
        )

    print("The following graphs were generated in directory %s: " % save_into_dir_path_name)
    for file_path_name in resulting_graphs_details:
        print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])

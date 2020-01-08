#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
import argparse
from collections import OrderedDict
import h5py
import src.plot as plot
from src import acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors, loadconfig


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name



CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX = "cytoplasmic_spread"

def cytoplasmic_spread(
    basic_h5_file_handler,
    genes,
    mime_type,
    molecule_type,
    colors,
    image_width,
    image_height,
    save_into_dir_path_name,
    ext_logger=None
    ):  #
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
    - ext_logger : external logger
    """
    assert os.path.isdir(save_into_dir_path_name)

    if ext_logger:
        ext_logger.debug("Generating cytoplasmic spread graph (%s)..." % molecule_type)

    graph_file_name = "%s_%s.png" % (molecule_type, CYTOPLASMIC_SPREAD_GRAPHS_SUFFIX)
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": mime_type
        }

    # for each gene in genes list, compute cytoplasmic spread
    cyt_spreads = []


    for gene in genes:
        if ext_logger:
            ext_logger.debug("Computing cytoplasmic spread for " + gene)
        image_list = helps.preprocess_image_list2(
            basic_h5_file_handler,
            '/' + molecule_type,
            gene
            )
        cyt_spreads.append(
            adsc.compute_cytoplasmic_spread(
                image_list,
                basic_h5_file_handler,
                image_width,
                image_height
                )
            )

    plot.bar_profile(cyt_spreads, genes, graph_file_path_name,colors)
    assert(os.path.isfile(graph_file_path_name))

    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}



CYTOPLASMIC_SPREAD_DYNAMIC_PROFILES_GRAPHS_PREFIX = "cyt_spread"

def cytoplasmic_spread_dynamic_profiles(
    basic_h5_file_handler,
    genes,
    proteins,
    colors,
    timepoints_num_mrna,
    timepoints_num_protein,
    image_width,
    image_height,
    mime_type,
    save_into_dir_path_name,
    ext_logger=None
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
    if ext_logger:
        ext_logger.debug("Generating cytoplasmic spread dynamic profiles...")

    graphs_details = []

    if ext_logger:
        ext_logger.debug("Extracting data...")

    data_generator = plot.data_extractor(
        genes,
        proteins,
        basic_h5_file_handler,
        adsc.compute_cytoplasmic_spread,
        basic_h5_file_handler,
        image_width,
        image_height
        )

    try:
        for mrna_data, protein_data, i in data_generator:
            gene = genes[i]
            colour = colors[i]

            if ext_logger:
                ext_logger.debug("Processing gene : %s..." % gene)

            try:
                graph_file_name = "%s_%s.png" % (CYTOPLASMIC_SPREAD_DYNAMIC_PROFILES_GRAPHS_PREFIX, gene)
                graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
                graph_metadata = {
                    "file_name": graph_file_name,
                    "mime_type": mime_type
                    }

                if ext_logger:
                    ext_logger.debug("Generating graph %s..." % graph_file_name)

                plot.dynamic_profiles(
                    mrna_data,
                    protein_data,
                    timepoints_num_mrna,
                    timepoints_num_protein,
                    gene,
                    colour,
                    'Time(hrs)',
                    'Cytoplasmic spread',
                    graph_file_path_name
                    )

                assert(os.path.isfile(graph_file_path_name))

                if ext_logger:
                    ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

                graphs_details.append({graph_file_path_name: graph_metadata})

            except Exception as e:
                if ext_logger:
                    ext_logger.exception(
                        "Could not generate cytoplasmic spread dynamic profile for gene : %s" % gene
                        )
                raise

    except Exception as e:
        if ext_logger:
            ext_logger.exception(
                "Got exception ! (maybe raised by data_generator ?)"
                )
        #raise

    return graphs_details


def main(
    save_into_dir_path_name,
    ext_logger=None
    ):
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    resulting_graphs_details_as_list = []

    configData = loadconfig(input_dir_name)

    # Compute bar plot cytoplasmic spread
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    mime_type = configData["PNG_IMAGES_MIME_TYPE"]
    image_width = configData["IMAGE_WIDTH"]
    image_height = configData["IMAGE_HEIGHT"]
    timepoints_num_protein=configData["TIMEPOINTS_NUM_PROTEIN"]
    timepoints_num_mrna=configData["TIMEPOINTS_NUM_MRNA"]

    colors=configData["COLORS"]

    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name, "r") as basic_h5_file_handler:

        graph_details = cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            mime_type=mime_type,
            molecule_type='mrna',
            colors=colors,
            image_width=image_width,
            image_height=image_height,
            save_into_dir_path_name=save_into_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list.append(graph_details)

        graph_details = cytoplasmic_spread(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            mime_type=mime_type,
            molecule_type='protein',
            colors=colors,
            image_width=image_width,
            image_height=image_height,
            save_into_dir_path_name=save_into_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list.append(graph_details)

        graphs_details = cytoplasmic_spread_dynamic_profiles(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            colors=colors,
            timepoints_num_mrna=timepoints_num_mrna,
            timepoints_num_protein=timepoints_num_protein,
            image_width=image_width,
            image_height=image_height,
            mime_type=mime_type,
            save_into_dir_path_name=save_into_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list += graphs_details

    resulting_graphs_details_as_odict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_odict.update(graph_details)

    return resulting_graphs_details_as_odict


if __name__ == "__main__":

    logger=enable_logger()

    save_into_dir_path_name = os.path.join(path.analysis_dir, "cytoplasmic_spread/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details = main(
        save_into_dir_path_name,logger
        )

    print("\nThe following graphs were generated in directory %s :\n" % save_into_dir_path_name)
    for file_path_name in resulting_graphs_details:
        print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])
    print("\n")

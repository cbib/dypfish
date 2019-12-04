#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
import argparse
from collections import OrderedDict
import h5py
import src.plot as plot
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors,loadconfig


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


CYTOPLASMIC_TOTAL_COUNT_GRAPHS_SUFFIX = "total_cytoplasmic_transcript"

def cytoplasmic_total_count(
    basic_h5_file_handler,
    genes,
    molecule_type,
    save_into_dir_path_name,
    mime_type,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    if ext_logger:
        ext_logger.debug("Generating cytoplasmic total count graphs (%s)..." % molecule_type)

    cytoplasmic_total = []
    for gene in genes:
        if ext_logger:
            ext_logger.debug("Processing gene %s..." % gene)

        image_list = helps.preprocess_image_list2(
            basic_h5_file_handler,
            '/' + molecule_type,
            gene
            )

        cytoplasmic_total.append(
            adsc.compute_cytoplasmic_total(
                image_list,
                basic_h5_file_handler,
                path_data=raw_images_dir_path_name
                )
            )

    graph_file_name = "%s_%s.png" % (molecule_type, CYTOPLASMIC_TOTAL_COUNT_GRAPHS_SUFFIX)
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": mime_type
        }

    plot.bar_profile(cytoplasmic_total, genes, graph_file_path_name)

    assert(os.path.isfile(graph_file_path_name))

    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


CYTOPLASMIC_TOTAL_COUNT_DYNAMIC_PROFILES_GRAPHS_PREFIX = "cyt_total"

def cytoplasmic_total_dynamic_profiles(
    basic_h5_file_handler,
    genes,
    proteins,
    timepoints_mrna,
    timepoints_protein,
    timepoints_num_mrna,
    timepoints_num_protein,
    save_into_dir_path_name,
    mime_type,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    if ext_logger:
        ext_logger.debug("Generating cytoplasmic total count dynamic profiles.")

    graphs_details = []

    data_generator = plot.data_extractor_generic(
        genes,
        proteins,
        timepoints_mrna,
        timepoints_protein,
        basic_h5_file_handler,
        adsc.compute_cytoplasmic_total,
        basic_h5_file_handler,
        raw_images_dir_path_name
        )

    try:

        for mrna_data, protein_data, i in data_generator:
            gene = genes[i]
            colour = plot_colors[i]

            if ext_logger:
                ext_logger.debug("Processing gene : %s..." % gene)

            try:
                graph_file_name = "%s_%s.png" % (CYTOPLASMIC_TOTAL_COUNT_DYNAMIC_PROFILES_GRAPHS_PREFIX, gene)
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
                    gene,
                    timepoints_num_mrna,
                    timepoints_num_protein,
                    colour,
                    'Time(hrs)',
                    'Cytoplasmic total',
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
        # raise

    return graphs_details



def main(
    save_into_dir_path_name,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):

    resulting_graphs_details_as_list = []

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts
    configData = loadconfig(input_dir_name)

    # Compute bar plot cytoplasmic spread
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    mime_type=configData["PNG_IMAGES_MIME_TYPE"]
    basic_file_name = configData["BASIC_FILE_NAME"]


    with h5py.File(path.analysis_data_dir+input_dir_name+'/'+basic_file_name, "r") as basic_h5_file_handler:

        graph_details = cytoplasmic_total_count(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            molecule_type='mrna',
            save_into_dir_path_name=save_into_dir_path_name,
            mime_type=mime_type,
            raw_images_dir_path_name=raw_images_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list.append(graph_details)

        graph_details = cytoplasmic_total_count(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            molecule_type='protein',
            save_into_dir_path_name=save_into_dir_path_name,
            mime_type=mime_type,
            raw_images_dir_path_name=raw_images_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list.append(graph_details)

        graphs_details = cytoplasmic_total_dynamic_profiles(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            timepoints_mrna=timepoints_mrna,
            timepoints_protein=timepoints_protein,
            timepoints_num_mrna=timepoints_num_mrna,
            timepoints_num_protein=timepoints_num_protein,
            save_into_dir_path_name=save_into_dir_path_name,
            mime_type=mime_type,
            raw_images_dir_path_name=raw_images_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list += graphs_details

    resulting_graphs_details_as_odict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_odict.update(graph_details)

    return resulting_graphs_details_as_odict


if __name__ == "__main__":

    enable_logger()

    save_into_dir_path_name = os.path.join(path.analysis_dir, "analysis_cytoplasmic_total_count/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details = main(
        save_into_dir_path_name
        )

    print("\nThe following graphs were generated in directory %s :\n" % save_into_dir_path_name)
    for file_path_name in resulting_graphs_details:
        print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])
    print("\n")

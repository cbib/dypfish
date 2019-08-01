#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
from collections import OrderedDict
import h5py
import src.plot as plot
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, plot_colors


PNG_IMAGES_MIME_TYPE = "image/png"

CYTOPLASMIC_TOTAL_COUNT_GRAPHS_SUFFIX = "total_cytoplasmic_transcript"

def cytoplasmic_total_count(
    basic_h5_file_handler,
    genes,
    molecule_type,
    save_into_dir_path_name,
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
        "mime_type": PNG_IMAGES_MIME_TYPE
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
    save_into_dir_path_name,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    if ext_logger:
        ext_logger.debug("Generating cytoplasmic total count dynamic profiles.")

    graphs_details = []

    data_generator = plot.data_extractor(
        genes,
        proteins,
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
                    "mime_type": PNG_IMAGES_MIME_TYPE
                    }

                if ext_logger:
                    ext_logger.debug("Generating graph %s..." % graph_file_name)

                plot.dynamic_profiles(
                    mrna_data,
                    protein_data,
                    gene,
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
    basic_h5_file_path_name,
    save_into_dir_path_name,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):

    resulting_graphs_details_as_list = []

    # Required descriptors: spots, IF, cell mask an height_map
    # Compute bar plot cytoplasmic total transcripts
    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]

    with h5py.File(basic_h5_file_path_name, "r") as basic_h5_file_handler:

        # graph_details = cytoplasmic_total_count(
        #     basic_h5_file_handler=basic_h5_file_handler,
        #     genes=genes,
        #     molecule_type='mrna',
        #     save_into_dir_path_name=save_into_dir_path_name,
        #     raw_images_dir_path_name=raw_images_dir_path_name,
        #     ext_logger=ext_logger
        #     )
        # resulting_graphs_details_as_list.append(graph_details)

        graph_details = cytoplasmic_total_count(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=proteins,
            molecule_type='protein',
            save_into_dir_path_name=save_into_dir_path_name,
            raw_images_dir_path_name=raw_images_dir_path_name,
            ext_logger=ext_logger
            )
        resulting_graphs_details_as_list.append(graph_details)

        graphs_details = cytoplasmic_total_dynamic_profiles(
            basic_h5_file_handler=basic_h5_file_handler,
            genes=genes,
            proteins=proteins,
            save_into_dir_path_name=save_into_dir_path_name,
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

    basic_h5_file_path_name = path.basic_file_path
    save_into_dir_path_name = os.path.join(path.analysis_dir, "analysis_cytoplasmic_total_count/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details = main(
        basic_h5_file_path_name,
        save_into_dir_path_name
        )

    print("\nThe following graphs were generated in directory %s :\n" % save_into_dir_path_name)
    for file_path_name in resulting_graphs_details:
        print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])
    print("\n")

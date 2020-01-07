#!/usr/bin/python
# encoding: UTF-8

import os
from collections import OrderedDict
import numpy as np
import h5py
import argparse
import src.plot as plot
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, plot_colors, check_dir,loadconfig

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name


PNG_IMAGES_MIME_TYPE = "image/png"

def peripheral_profile(
    mean_profiles,
    timepoint,
    genes,
    molecule_type,
    num_contours,
    colors,
    save_into_dir_path_name,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    timepoint = timepoint if timepoint else 'All'
    graph_file_name = molecule_type[0]+ '_peripheral_fraction_profile.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
        }
    # normalize by Gapdh profile (here index 2 in mean profiles array)
    mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
    figure_title = ' peripheral profile ('+ timepoint+')'

    plot.profile(
        mean_profiles,
        genes,
        num_contours,
        graph_file_path_name,
        figure_title,
        colors,
        True
        )

    assert(os.path.isfile(graph_file_path_name))

    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}



def peripheral_profiles(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    molecule_type,
    genes,
    num_contours,
    colors,
    compute_peripheral_fraction_profiles,
    timepoints,
    save_into_dir_path_name,
    ext_logger=None
    ):

    assert os.path.isdir(save_into_dir_path_name)

    graphs_details = []

    if not timepoints:
        timepoints=[False]
    for timepoint in timepoints:
        mean_profiles = []
        for gene in genes:
            if timepoint:
                image_list = helps.preprocess_image_list3(
                    secondary_h5_file_handler,
                    molecule_type,
                    gene,
                    [timepoint]
                    )
            else:
                image_list = helps.preprocess_image_list2(
                secondary_h5_file_handler,
                molecule_type[0],
                gene
                )
            profiles=compute_peripheral_fraction_profiles(
                basic_h5_file_handler,
                secondary_h5_file_handler,
                image_list
                )
            mean_profiles.append(np.average(profiles, axis=0))

        graph_details = peripheral_profile(
            mean_profiles=mean_profiles,
            timepoint=timepoint,
            genes=genes,
            molecule_type=molecule_type,
            num_contours=num_contours,
            colors=colors,
            save_into_dir_path_name=save_into_dir_path_name,
            ext_logger=ext_logger
            )

        graphs_details.append(graph_details)

    return graphs_details


def histogram_peripheral_profile(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    genes,
    proteins,
    colors,
    periph_fraction_cst,
    save_into_dir_path_name,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    graphs_details = []

    try:

        molecule_type = ['mrna']

        if ext_logger:
            ext_logger.debug("generating histogram_peripheral_profile profile for molecule : %s..." % molecule_type[0])

        slashed_molecule_type = [("/%s" % m) for m in molecule_type]
        periph_fraction = []
        for gene in genes:
            image_list = helps.preprocess_image_list2(
                basic_h5_file_handler,
                slashed_molecule_type[0],
                gene
                )
            periph_fraction.append(
                adsc.compute_mrna_periph_fraction(
                    image_list,
                    basic_h5_file_handler,
                    secondary_h5_file_handler,
                    periph_fraction_cst,
                    )
                )

        graph_file_name = molecule_type[0] +'_peripheral_fraction_'+ str(periph_fraction_cst)+'.png'
        graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
        graph_metadata = {
            "file_name": graph_file_name,
            "mime_type": PNG_IMAGES_MIME_TYPE
            }

        plot.bar_profile(periph_fraction, genes, graph_file_path_name)
        assert(os.path.isfile(graph_file_path_name))

        if ext_logger:
            ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

        graphs_details.append({graph_file_path_name: graph_metadata})


    except Exception as e:
        if ext_logger:
            ext_logger.exception(
                "Could not generate histogram_peripheral_profile profile for molecule : %s" % molecule_type[0]
                )
        raise


    try:

        molecule_type = ['protein']

        if ext_logger:
            ext_logger.debug("generating histogram_peripheral_profile profile for molecule : %s..." % molecule_type[0])

        slashed_molecule_type = [("/%s" % m) for m in molecule_type]
        periph_fraction = []
        for protein in proteins:
            image_list = helps.preprocess_image_list2(
                basic_h5_file_handler,
                slashed_molecule_type[0],
                protein
                )
            periph_fraction.append(
                adsc.compute_protein_periph_fraction(
                    image_list,
                    basic_h5_file_handler,
                    secondary_h5_file_handler,
                    periph_fraction_cst,
                    )
                )

        graph_file_name = molecule_type[0] +'_peripheral_fraction.png'
        graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
        graph_metadata = {
            "file_name": graph_file_name,
            "mime_type": PNG_IMAGES_MIME_TYPE
            }

        plot.bar_profile(periph_fraction, proteins, graph_file_path_name)

        assert(os.path.isfile(graph_file_path_name))

        if ext_logger:
            ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

        graphs_details.append({graph_file_path_name: graph_metadata})

    except Exception as e:
        if ext_logger:
            ext_logger.exception(
                "Could not generate histogram_peripheral_profile profile for molecule : %s" % molecule_type[0]
                )
        raise

    return graphs_details



def peripheral_fraction_dynamic_profile(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    genes,
    proteins,
    timepoints_num_mrna,
    timepoints_num_protein,
    colors,
    periph_fraction_cst,
    mrna_tp,
    protein_tp,
    save_into_dir_path_name,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    graphs_details = []

    data_generator = plot.data_extractor_generic(
        genes,
        proteins,
        mrna_tp,
        protein_tp,
        secondary_h5_file_handler,
        adsc.compute_generic_periph_fraction,
        basic_h5_file_handler,
        secondary_h5_file_handler,
        periph_fraction_cst,
        )

    try:

        for mrna_data, protein_data, i in data_generator:
            gene = genes[i]
            plot_color = plot_colors[i]

            if ext_logger:
                ext_logger.debug("generating peripheral fraction dynamic profile for gene : %s..." % gene)

            try:

                graph_file_name = 'peripheral_fraction_' + genes[i] + '.png'
                graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
                graph_metadata = {
                    "file_name": graph_file_name,
                    "mime_type": PNG_IMAGES_MIME_TYPE
                    }

                plot.dynamic_profiles(
                    mrna_data,
                    protein_data,
                    timepoints_num_mrna,
                    timepoints_num_protein,
                    gene,
                    plot_color,
                    'Time(hrs)',
                    'Peripheral fraction',
                    graph_file_path_name
                    )

                assert(os.path.isfile(graph_file_path_name))

                if ext_logger:
                    ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

                graphs_details.append({graph_file_path_name: graph_metadata})

            except Exception as e:
                if ext_logger:
                    ext_logger.exception(
                        "Could not generate peripheral fraction dynamic profile for gene : %s" % gene
                        )
                raise

    except Exception as e:
        if ext_logger:
            ext_logger.exception(
                "Got exception ! (maybe raised by data_generator ?)"
                )
        raise

    return graphs_details



def main(
    save_into_dir_path_name
    ):
    assert os.path.isdir(save_into_dir_path_name)

    resulting_graphs_details_as_list = []
    errors = []

    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots
    ## Build peripheral profile plot either for each or for all timepoint

    configData = loadconfig(input_dir_name)
    genes = configData["GENES"]
    proteins = configData["PROTEINS"]
    timepoints_mrna = configData["TIMEPOINTS_MRNA"]
    timepoints_protein = configData["TIMEPOINTS_PROTEIN"]
    timepoints_num_mrna = configData["TIMEPOINTS_NUM_MRNA"]
    timepoints_num_protein = configData["TIMEPOINTS_NUM_PROTEIN"]
    periph_fraction_cst = configData["PERIPHERAL_FRACTION_THRESHOLD"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    num_contours = configData["NUM_CONTOURS"]

    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name, "r") as basic_h5_file_handler, \
         h5py.File(path.data_dir+input_dir_name+'/'+secondary_file_name, "r") as secondary_h5_file_handler:

        # Section to build peripheral profile graph (figure 2B in dypFISH Paper, only mRNA figure is shown in article)
        try:
            graphs_details = peripheral_profiles(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                molecule_type=['mrna'],
                genes=genes,
                num_contours=num_contours,
                colors=plot_colors,
                compute_peripheral_fraction_profiles=adsc.compute_mrna_peripheral_fraction_profiles_3D,
                timepoints=None,
                save_into_dir_path_name=save_into_dir_path_name
            )
            resulting_graphs_details_as_list += graphs_details

        except Exception as e:
            errors.append("Could not generate peripheral profiles graphs.")
            raise

        try:
            graphs_details = peripheral_profiles(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                molecule_type=['protein'],
                genes=proteins,
                num_contours=num_contours,
                colors=plot_colors,
                compute_peripheral_fraction_profiles=adsc.compute_protein_peripheral_fraction_profiles_3D,
                timepoints=None,
                save_into_dir_path_name=save_into_dir_path_name
            )
            resulting_graphs_details_as_list += graphs_details

        except Exception as e:
            errors.append("Could not generate peripheral profiles graphs.")
            raise



        # Section to compute bar plot peripheral fraction (change PERIPHERAL_FRACTION_THRESHOLD constant in related config.json)
        try:
            graphs_details = histogram_peripheral_profile(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                genes=genes,
                proteins=proteins,
                colors=plot_colors,
                periph_fraction_cst=periph_fraction_cst,
                save_into_dir_path_name=save_into_dir_path_name
                )
            resulting_graphs_details_as_list += graphs_details

        except Exception as e:
            errors.append("Could not generate histogram peripheral profile graphs.")
            raise

        # Section to produce plot interpolation (dynamic profile) of peripheral fraction by timepoint
        try:
            graphs_details = peripheral_fraction_dynamic_profile(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                genes=genes,
                proteins=proteins,
                timepoints_num_mrna=timepoints_num_mrna,
                timepoints_num_protein=timepoints_num_protein,
                colors=plot_colors,
                periph_fraction_cst=periph_fraction_cst,
                mrna_tp=timepoints_mrna,
                protein_tp=timepoints_protein,
                save_into_dir_path_name=save_into_dir_path_name,
            )
            resulting_graphs_details_as_list += graphs_details
        except Exception as e:
            errors.append("Could not generate peripheral fraction dynamic profile graphs.")
            raise


    resulting_graphs_details_as_odict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_odict.update(graph_details)

    return resulting_graphs_details_as_odict, errors


if __name__ == "__main__":

    enable_logger()

    save_into_dir_path_name = os.path.join(path.analysis_dir, "peripheral_fraction_profile/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details, errors = main(
        save_into_dir_path_name=save_into_dir_path_name
        )

    if resulting_graphs_details:
        print("\nThe following graphs were generated in directory %s :\n" % save_into_dir_path_name)
        for file_path_name in resulting_graphs_details:
            print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])
        print("\n")

    if errors:
        print("\nThe following errors were encountered :\n")
        for error_msg in errors:
            print("-> %s" % error_msg)
        print("\n")

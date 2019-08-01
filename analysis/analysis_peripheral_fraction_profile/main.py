#!/usr/bin/python
# encoding: UTF-8

import os
from collections import OrderedDict
import numpy as np
import h5py
import src.plot as plot
import src.constants as constants
import src.acquisition_descriptors as adsc
import src.path as path
import src.helpers as helps
from src.utils import enable_logger, cell_type_micropatterned, plot_colors, check_dir

PNG_IMAGES_MIME_TYPE = "image/png"

def peripheral_profile(
    mean_profiles,
    timepoint,
    genes,
    colors,
    save_into_dir_path_name,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    timepoint = timepoint if timepoint else 'All'
    graph_file_name = 'peripheral_fraction_' + timepoint + '_timepoint' + str(constants.NUM_CONTOURS) + 'contours.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)
    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
        }

    mean_profiles = mean_profiles / np.matlib.repmat(mean_profiles[2], len(genes), 1)
    figure_title = ' peripheral profile ('+ timepoint+')'

    plot.profile(
        mean_profiles,
        genes,
        constants.NUM_CONTOURS,
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
    secondary_h5_file_handler,
    molecule_type,
    genes,
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
        print(timepoint)
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
                secondary_h5_file_handler,
                image_list
                )
            mean_profiles.append(np.average(profiles, axis=0))

        graph_details = peripheral_profile(
            mean_profiles=mean_profiles,
            timepoint=timepoint,
            genes=genes,
            colors=colors,
            save_into_dir_path_name=save_into_dir_path_name,
            ext_logger=ext_logger
            )

        graphs_details.append(graph_details)

    return graphs_details


def peripheral_fraction_profile(
    secondary_h5_file_handler,
    molecule_type,
    genes,
    fraction,
    colors,
    basic_h5_file_handler,
    save_into_dir_path_name,
    raw_images_dir_path_name=None,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    fractions = []
    for gene in genes:
        print(gene)
        image_list = helps.preprocess_image_list2(
            secondary_h5_file_handler,
            molecule_type[0],
            gene
            )
        fractions.append(
            adsc.build_histogram_periph_fraction(
                secondary_h5_file_handler,
                image_list,
                fraction,
                raw_images_dir_path_name,
                basic_h5_file_handler
                )
            )
    print(fractions)
    graph_file_name = 'peripheral_fraction_'+str(fraction)+'.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name,graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
        }

    #plot.fraction_profile(fractions, fraction, genes, fig_file_path_name, colors)
    plot.bar_profile(fractions, genes, graph_file_path_name)

    assert(os.path.isfile(graph_file_path_name))

    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def histogram_peripheral_profile(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    genes,
    proteins,
    colors,
    raw_images_dir_path_name,
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
            print(gene)
            image_list = helps.preprocess_image_list2(
                basic_h5_file_handler,
                slashed_molecule_type[0],
                gene
                )
            periph_fraction.append(
                adsc.compute_periph_fraction(
                    image_list,
                    basic_h5_file_handler,
                    secondary_h5_file_handler,
                    constants.PERIPHERAL_FRACTION_THRESHOLD,
                    raw_images_dir_path_name
                    )
                )

        graph_file_name = molecule_type[0] +'_peripheral_fraction.png'
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
                adsc.compute_periph_fraction(
                    image_list,
                    basic_h5_file_handler,
                    secondary_h5_file_handler,
                    constants.PERIPHERAL_FRACTION_THRESHOLD,
                    raw_images_dir_path_name
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
    colors,
    mrna_tp,
    protein_tp,
    raw_images_dir_path_name,
    save_into_dir_path_name,
    ext_logger=None
    ):
    assert os.path.isdir(save_into_dir_path_name)

    graphs_details = []

    data_generator = plot.data_extractor(
        genes,
        proteins,
        secondary_h5_file_handler,
        adsc.compute_periph_fraction,
        basic_h5_file_handler,
        secondary_h5_file_handler,
        constants.PERIPHERAL_FRACTION_THRESHOLD,
        raw_images_dir_path_name
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
    basic_h5_file_handler,
    secondary_h5_file_handler,
    raw_images_dir_path_name,
    save_into_dir_path_name
    ):
    assert os.path.isdir(save_into_dir_path_name)

    resulting_graphs_details_as_list = []
    errors = []

    # Required descriptors: spots_peripheral_distance, height_map, zero_level and spots

    ## Build peripheral profile plot either for each or for all timepoint
    molecule_type = ['mrna']
    slashed_molecule_type = [("/%s" % m) for m in molecule_type]

    genes = ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"]
    proteins = ["beta_actin", "arhgdia", "gapdh", "pard3"]
    timepoints = ["2h", "3h", "4h", "5h"]
    timepoints_protein = ["2h", "3h", "5h", "7h"]

    with h5py.File(basic_h5_file_path_name, "r") as basic_h5_file_handler, \
         h5py.File(secondary_h5_path_name, "r") as secondary_h5_file_handler:

        # try:
        #     graphs_details = peripheral_profiles(
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         molecule_type=slashed_molecule_type,
        #         genes=genes,
        #         colors=plot_colors,
        #         compute_peripheral_fraction_profiles=adsc.compute_peripheral_fraction_profiles_3D,
        #         timepoints=None,
        #         save_into_dir_path_name=save_into_dir_path_name
        #         )
        #     resulting_graphs_details_as_list += graphs_details
        #
        # except Exception as e:
        #     errors.append("Could not generate peripheral profiles graphs.")

        # Section to build peripheral profile fraction 10 and 30
        try:
            graph_details = peripheral_fraction_profile(
                secondary_h5_file_handler=secondary_h5_file_handler,
                molecule_type=slashed_molecule_type,
                genes=genes,
                fraction=10,
                colors=plot_colors,
                basic_h5_file_handler=basic_h5_file_handler,
                save_into_dir_path_name=save_into_dir_path_name
                )
            resulting_graphs_details_as_list.append(graph_details)

        except Exception as e:
            errors.append("Could not generate peripheral fraction profile graph.")

        # Section to compute bar plot peripheral fraction
        # try:
        #     graphs_details = histogram_peripheral_profile(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         proteins=proteins,
        #         colors=plot_colors,
        #         raw_images_dir_path_name=raw_images_dir_path_name,
        #         save_into_dir_path_name=save_into_dir_path_name
        #         )
        #     resulting_graphs_details_as_list += graphs_details

        except Exception as e:
            errors.append("Could not generate histogram peripheral profile graphs.")

        # Section to produce plot interpolation (dynamic profile) of peripheral fraction by timepoint
        # try:
        #     graphs_details = peripheral_fraction_dynamic_profile(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         proteins=proteins,
        #         colors=plot_colors,
        #         mrna_tp=timepoints,
        #         protein_tp=timepoints_protein,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         raw_images_dir_path_name=raw_images_dir_path_name
        #         )
        #     resulting_graphs_details_as_list += graphs_details
        # except Exception as e:
        #     errors.append("Could not generate peripheral fraction dynamic profile graphs.")


    resulting_graphs_details_as_odict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_odict.update(graph_details)

    return resulting_graphs_details_as_odict, errors


if __name__ == "__main__":

    enable_logger()

    basic_h5_file_path_name = path.basic_file_path
    secondary_h5_path_name = path.secondary_file_path
    raw_images_dir_path_name = path.path_data

    save_into_dir_path_name = os.path.join(path.analysis_dir, "analysis_peripheral_fraction_profile/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details, errors = main(
        basic_h5_file_handler=basic_h5_file_path_name,
        secondary_h5_file_handler=secondary_h5_path_name,
        raw_images_dir_path_name=raw_images_dir_path_name,
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

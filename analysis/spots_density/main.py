#!/usr/bin/python
# encoding: UTF-8
# author: Benjamin Dartigues

import os
from collections import OrderedDict
import h5py
import math
import argparse
import numpy as np
import pandas as pd
import src.plot as plot
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
from src.utils import check_dir, loadconfig

import src.helpers as helps
from src.utils import enable_logger, plot_colors


parser = argparse.ArgumentParser()
parser.add_argument("--input_dir_name", "-i", help='input dir where to find h5 files and configuration file', type=str)
args = parser.parse_args()
input_dir_name = args.input_dir_name

PNG_IMAGES_MIME_TYPE = "image/png"

SPOTS_DENSITY_GRAPHS_SUFFIX = "spots_density"


def get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints
    ):
    cell_areas = []
    nucleus_areas = []
    transcripts_count = []
    for timepoint in timepoints:

        for gene in genes:

            transcript_count=[]
            if timepoint:
                image_list= helps.preprocess_image_list3(
                    basic_h5_file_handler,
                    molecule_type,
                    gene,
                    [timepoint]
                    )
            else:
                image_list = helps.preprocess_image_list2(
                    basic_h5_file_handler,
                    molecule_type[0],
                    gene
                )
            transcripts_count+=[len(idsc.get_spots(basic_h5_file_handler, x)) for x in image_list]
            cell_area=[]
            nucleus_area=[]
            for x in image_list:
                nucleus_n = idsc.count_nucleus(basic_h5_file_handler, x)
                cell_area.append(idsc.get_cell_area(secondary_h5_file_handler, x) / nucleus_n)
                nucleus_area.append(idsc.get_nucleus_area(secondary_h5_file_handler, x) / nucleus_n)
            cell_areas.append(cell_area)
            nucleus_areas.append(nucleus_area)




    return nucleus_areas,cell_areas,transcripts_count


def get_data_wo_nucleus(basic_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        size_coeff
    ):
    cell_areas = []
    nucleus_areas = []
    transcripts_count = []
    for timepoint in timepoints:
        for gene in genes:
            transcript_count=[]
            if timepoint:
                image_list= helps.preprocess_image_list3(
                    basic_h5_file_handler,
                    molecule_type,
                    gene,
                    [timepoint]
                    )
            else:
                image_list = helps.preprocess_image_list2(
                    basic_h5_file_handler,
                    molecule_type[0],
                    gene
                )
            transcript_count=[]
            cyt_transcript_count=[]
            cell_area=[]
            cell_area_wo_nucleus = []
            nucleus_area=[]

            image_cleaned_list=[]
            for x in image_list:


                if (len(get_cytoplasmic_spots(basic_h5_file_handler, x))>=20):

                    image_cleaned_list.append(x)
                    nucleus_n = idsc.count_nucleus(basic_h5_file_handler, x)
                    transcript_count.append(len(idsc.get_spots(basic_h5_file_handler, x))/nucleus_n)
                    #cell_area.append(idsc.get_cell_area(secondary_h5_file_handler, x) / nucleus_n)
                    cell_area_wo_nucleus.append(get_cell_area_wo_nucleus(basic_h5_file_handler, x, size_coeff) / nucleus_n)
                    cyt_transcript_count.append(len(get_cytoplasmic_spots(basic_h5_file_handler, x)))
                    nucleus_mask = idsc.get_nucleus_mask(basic_h5_file_handler, x)

                    nuc_area=nucleus_mask.sum() * math.pow((1 / size_coeff), 2)
                    nucleus_area.append(nuc_area / nucleus_n)

            transcripts_count.append(cyt_transcript_count)
            cell_areas.append(cell_area_wo_nucleus)
            nucleus_areas.append(nucleus_area)

            df = pd.DataFrame(index=image_cleaned_list)
            df["transcript_count"] = cyt_transcript_count
            df["area"] = cell_area_wo_nucleus
            df["density"] = np.array(cyt_transcript_count)/np.array(cell_area_wo_nucleus)
            tmp=np.array(cyt_transcript_count)/np.array(cell_area_wo_nucleus)
            tmp2=tmp-np.median(tmp)
            df["ecart"] = tmp2
            graph_file_name = gene+'.csv'
            graph_file_path_name = os.path.join(path.analysis_dir+"spots_density/figures/", graph_file_name)
            df.to_csv(graph_file_path_name)


    return nucleus_areas,cell_areas,transcripts_count

def get_cytoplasmic_spots(basic_h5_file_handler,image):
    spots=idsc.get_spots(basic_h5_file_handler, image)
    nucleus_mask = idsc.get_nucleus_mask(basic_h5_file_handler, image)
    cyt_spots=[]
    for spot in spots:
        if nucleus_mask[spot[0],spot[1]]==0:

            cyt_spots.append(spot)

    return cyt_spots

def get_cell_area_wo_nucleus(basic_h5_file_handler, image, size_coeff):
        cell_mask=idsc.get_cell_mask(basic_h5_file_handler,image)
        nucleus_mask=idsc.get_nucleus_mask(basic_h5_file_handler,image)
        cell_mask[nucleus_mask==1]=0
        cell_area= cell_mask.sum() * math.pow((1 / size_coeff), 2)
        return cell_area


def compare_nucleus_area( basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas, transcripts_count = get_data(basic_h5_file_handler,
                                             secondary_h5_file_handler,
                                             genes,
                                             molecule_type,
                                             timepoints)
    graph_file_name = 'nucleus_size_micropatterned_vs_standard.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.boxplot2(nucleus_areas, graph_file_path_name, 900)
    assert (os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_cell_area( basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas, transcripts_count = get_data(basic_h5_file_handler,
                                             secondary_h5_file_handler,
                                             genes,
                                             molecule_type,
                                             timepoints)
    graph_file_name = 'cell_size_micropatterned_vs_standard.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.boxplot2(cell_areas, graph_file_path_name, 2500)
    assert (os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}

def compare_spots_density(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)
    #print(save_into_dir_path_name)
    graph_file_name = 'mic_vs_cult.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.linear_regression(cell_areas,transcripts_count,graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_micropatterned_wo_nucleus(
        basic_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        size_coeff,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data_wo_nucleus(basic_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        size_coeff)

    graph_file_name = 'cell_area_wo_nucleus_vs_transcript_count_micropatterned_distplot.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }


    plot.sns_linear_regression(cell_areas[0],transcripts_count[0],colors[0],graph_file_path_name)

    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_micropatterned_wo_nucleus_scratch(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data_wo_nucleus(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)

    graph_file_name = 'cell_area_wo_nucleus_vs_transcript_count_scratch_distplot.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }


    plot.sns_linear_regression(cell_areas[2],transcripts_count[2],colors[2],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}



def compare_spots_density_micropatterned(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)


    graph_file_name = 'cell_area_vs_transcript_count_micropatterned.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.sns_linear_regression(cell_areas[0],transcripts_count[0],colors[0],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_standard(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)

    graph_file_name = 'cell_area_vs_transcript_count_standard.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }

    plot.sns_linear_regression(cell_areas[1],transcripts_count[1],colors[1],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_standard_wo_nucleus(
        basic_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        size_coeff,
        save_into_dir_path_name,
        ext_logger=None
    ):
    nucleus_areas,cell_areas,transcripts_count=get_data_wo_nucleus(basic_h5_file_handler,
        genes,
        molecule_type,
        timepoints,
        size_coeff)
    print(nucleus_areas)

    graph_file_name = 'cell_area_wo_nucleus_vs_transcript_count_standard_distplot.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }


    plot.sns_linear_regression(cell_areas[1],transcripts_count[1],colors[1],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_histogram_micropatterned(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):

    nucleus_areas,cell_areas,transcripts_count=get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)
    graph_file_name = 'cell_area_vs_transcript_count_micropatterned_histogram.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.sns_histogram(cell_areas[0],transcripts_count[0],colors[0],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}


def compare_spots_density_histogram_cultured(
        basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        colors,
        timepoints,
        save_into_dir_path_name,
        ext_logger=None
    ):

    nucleus_areas,cell_areas,transcripts_count=get_data(basic_h5_file_handler,
        secondary_h5_file_handler,
        genes,
        molecule_type,
        timepoints)

    graph_file_name = 'cell_area_vs_transcript_count_standard_histogram.png'
    graph_file_path_name = os.path.join(save_into_dir_path_name, graph_file_name)

    graph_metadata = {
        "file_name": graph_file_name,
        "mime_type": PNG_IMAGES_MIME_TYPE
    }
    plot.sns_histogram(cell_areas[1],transcripts_count[1],colors[1],graph_file_path_name)
    assert(os.path.isfile(graph_file_path_name))
    if ext_logger:
        ext_logger.debug("Graph generated on path : %s" % graph_file_path_name)

    return {graph_file_path_name: graph_metadata}




def main(
    save_into_dir_path_name,
    ext_logger=None
    ):

    assert os.path.isdir(save_into_dir_path_name)
    resulting_graphs_details_as_list = []
    errors = []

    # Required descriptors: spots,  cell mask, height_map
    # Compute cell volume vs total transcripts between arhgdia micropatterned cells and arhgdia cultured cells

    configData = loadconfig(input_dir_name)
    genes = configData["GENES"]
    basic_file_name = configData["BASIC_FILE_NAME"]
    secondary_file_name = configData["SECONDARY_FILE_NAME"]
    size_coeff = configData["SIZE_COEFFICIENT"]

    timepoints = [""]
    colors = ['#0080ff', '#646464']
    molecule_type = ['mrna']

    slashed_molecule_type = [("/%s" % m) for m in molecule_type]

    with h5py.File(path.data_dir+input_dir_name+'/'+basic_file_name,  "r") as basic_h5_file_handler, h5py.File(path.data_dir+input_dir_name+'/'+secondary_file_name,  "r") as secondary_h5_file_handler:


        # try:
        #
        #     graph_details = compare_spots_density_histogram_micropatterned(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         molecule_type=molecule_type,
        #         colors=colors,
        #         timepoints=timepoints,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         ext_logger=ext_logger
        #     )
        #     resulting_graphs_details_as_list.append(graph_details)
        #
        #
        # except Exception as e:
        #     errors.append("Could not generate spots density graphs.")


        # try:
        #
        #     graph_details = compare_spots_density_histogram_cultured(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         molecule_type=molecule_type,
        #         colors=colors,
        #         timepoints=timepoints,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         ext_logger=ext_logger
        #     )
        #     resulting_graphs_details_as_list.append(graph_details)
        #
        #
        # except Exception as e:
        #     errors.append("Could not generate spots density graphs.")


        # try:
        #     graph_details = compare_spots_density(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         molecule_type=molecule_type,
        #         timepoints=timepoints,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         ext_logger=ext_logger
        #     )
        #     resulting_graphs_details_as_list.append(graph_details)
        #
        #
        # except Exception as e:
        #     errors.append("Could not generate spots density graphs.")
        #
        #
        # try:
        #     graph_details = compare_spots_density_micropatterned(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         molecule_type=molecule_type,
        #         colors=colors,
        #         timepoints=timepoints,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         ext_logger=ext_logger
        #     )
        #     resulting_graphs_details_as_list.append(graph_details)
        #
        # except Exception as e:
        #     print("Could not generate spots density graphs for micropatterned cells.")
        #     errors.append("Could not generate spots density graphs for micropatterned cells.")
        #
        # try:
        #     graph_details = compare_spots_density_standard(
        #         basic_h5_file_handler=basic_h5_file_handler,
        #         secondary_h5_file_handler=secondary_h5_file_handler,
        #         genes=genes,
        #         colors=colors,
        #         molecule_type=molecule_type,
        #         timepoints=timepoints,
        #         save_into_dir_path_name=save_into_dir_path_name,
        #         ext_logger=ext_logger
        #     )
        #     resulting_graphs_details_as_list.append(graph_details)
        #
        #
        # except Exception as e:
        #     print("Could not generate spots density graphs for standard culture cells.")
        #
        #     errors.append("Could not generate spots density graphs for standard culture cells.")
        #

        try:
            graph_details = compare_spots_density_micropatterned_wo_nucleus(
                basic_h5_file_handler=basic_h5_file_handler,
                genes=genes,
                molecule_type=molecule_type,
                colors=colors,
                timepoints=timepoints,
                size_coeff=size_coeff,
                save_into_dir_path_name=save_into_dir_path_name,
                ext_logger=ext_logger
            )
            resulting_graphs_details_as_list.append(graph_details)

        except Exception as e:
            print("Could not generate spots density graphs for micropatterned cells.")
            errors.append("Could not generate spots density graphs for micropatterned cells.")

        try:
            graph_details = compare_spots_density_standard_wo_nucleus(
                basic_h5_file_handler=basic_h5_file_handler,
                genes=genes,
                colors=colors,
                molecule_type=molecule_type,
                timepoints=timepoints,
                size_coeff=size_coeff,
                save_into_dir_path_name=save_into_dir_path_name,
                ext_logger=ext_logger
            )
            resulting_graphs_details_as_list.append(graph_details)


        except Exception as e:
            errors.append("Could not generate spots density graphs for standard culture cells.")


        try:
            graph_details = compare_cell_area(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                genes=genes,
                molecule_type=molecule_type,
                timepoints=timepoints,
                save_into_dir_path_name=save_into_dir_path_name,
                ext_logger=ext_logger
            )
            resulting_graphs_details_as_list.append(graph_details)


        except Exception as e:
            errors.append("Could not generate cell area comparison graphs.")

        try:
            graph_details = compare_nucleus_area(
                basic_h5_file_handler=basic_h5_file_handler,
                secondary_h5_file_handler=secondary_h5_file_handler,
                genes=genes,
                molecule_type=molecule_type,
                timepoints=timepoints,
                save_into_dir_path_name=save_into_dir_path_name,
                ext_logger=ext_logger
            )
            resulting_graphs_details_as_list.append(graph_details)


        except Exception as e:
            errors.append("Could not generate nucleus area comparison graphs.")

    resulting_graphs_details_as_odict = OrderedDict()
    for graph_details in resulting_graphs_details_as_list:
        resulting_graphs_details_as_odict.update(graph_details)

    return resulting_graphs_details_as_odict,errors




if __name__ == "__main__":

    enable_logger()

    save_into_dir_path_name = os.path.join(path.analysis_dir, "spots_density/figures/")
    if not os.path.isdir(save_into_dir_path_name):
        os.mkdir(save_into_dir_path_name)

    resulting_graphs_details,errors = main(
        save_into_dir_path_name
    )
    if resulting_graphs_details:
        print("\nThe following graphs were generated in directory %s :\n" % save_into_dir_path_name)
        for file_path_name in resulting_graphs_details:
            print("-> %s" % resulting_graphs_details[file_path_name]["file_name"])
        print("\n")
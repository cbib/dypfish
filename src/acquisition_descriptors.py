#!/usr/bin/env python2
# encoding: UTF-8

import matplotlib
try:
    import Tkinter
except:
    matplotlib.use('Agg')
    # required on headless Linux servers
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d

import image_descriptors
from src.utils import enable_logger

logger = enable_logger()

def compute_fraction_profile(file_handler, image_list, isoline):
    gene_array = []
    for image in image_list:
        spots_peripheral_distance = image_descriptors.get_spots_peripheral_distance_2D(file_handler, image)
        gene_array.append(float(len(np.where(spots_peripheral_distance <= (isoline))[0]))/float(len(spots_peripheral_distance)))
    return gene_array


def compute_peripheral_fraction_profiles_2D(file_handler, image_list):
    profiles = []
    for image in image_list:
        spots_peripheral_distance = image_descriptors.get_spots_peripheral_distance_2D(file_handler, image)
        peripheral_profile = np.zeros(100)
        for i in range(0, 100):
            peripheral_profile[i] = float(len(np.where(spots_peripheral_distance <= (i + 1))[0]))/float(len(spots_peripheral_distance))
        profiles.append(peripheral_profile)

    return profiles


def compute_peripheral_fraction_profiles_3D(file_handler, image_list):
    profiles = []
    for image in image_list:
        print(image)
        spots_peripheral_distance = image_descriptors.get_spots_peripheral_distance(file_handler, image)
        print(spots_peripheral_distance)
        peripheral_profile = np.zeros(100)
        for i in range(0, 100):
            peripheral_profile[i] = float(len(np.where(spots_peripheral_distance <= (i + 1))[0]))/float(len(spots_peripheral_distance))
        profiles.append(peripheral_profile)

    return profiles


def compute_periph_fraction(image_list,basic_file_handler,secondary_file_handler, fraction, path_data):
    arr=[]
    for image in image_list:

        if 'mrna' in image:
            print(image)
            spots_peripheral_distance = image_descriptors.get_spots_peripheral_distance_2D(secondary_file_handler, image)
            arr.append(float(len(spots_peripheral_distance[spots_peripheral_distance <= fraction]))/float(len(spots_peripheral_distance)))
        else:
            cell_mask=image_descriptors.get_cell_mask(basic_file_handler,image)
            nucleus_mask=image_descriptors.get_nucleus_mask(basic_file_handler,image)
            cell_mask_distance_map=image_descriptors.get_cell_mask_distance_map(secondary_file_handler,image)
            if "tp" in image:
                molecule = image.split("/")[1]
                gene = image.split("/")[2].split("_")[0]
                timepoint = image.split("/")[2].split("_")[1]
                image_n = image.split("/")[4]
                IF_image_path = path_data + '/' + molecule + '/' + gene + '/' + timepoint + '/' + image_n + '/IF.tif'
                try:
                    with open(IF_image_path):
                        pass
                except IOError:
                    continue
                print(IF_image_path)
                IF= image_descriptors.get_IF_image_z_summed(molecule,gene,timepoint,image_n,path_data)
            else:
                IF= image_descriptors.get_IF(basic_file_handler,image)

            IF=np.multiply(IF,cell_mask)
            IF_periph=IF[(cell_mask_distance_map<=10) & (cell_mask_distance_map >0)]
            IF_summed=IF[nucleus_mask==0].sum()
            IF_summed_periph=IF_periph.sum()
            periph_frac=float(IF_summed_periph)/float(IF_summed)
            arr.append(periph_frac)
    return arr


def build_histogram_periph_fraction(sec_file_handler,image_list,fraction,path_data,basic_file_handler):
    arr=[]
    for image in image_list:
        if 'mrna' in image:
            spots_peripheral_distance = image_descriptors.get_spots_peripheral_distance(sec_file_handler, image)
            arr.append(float(len(spots_peripheral_distance[spots_peripheral_distance <= fraction]))/float(len(spots_peripheral_distance)))
        else:
            nucleus_mask = image_descriptors.get_nucleus_mask(basic_file_handler, image)
            IF= image_descriptors.get_IF(basic_file_handler,image)
            cell_mask_distance_map = image_descriptors.get_cell_mask_distance_map(sec_file_handler, image)
            IF_periph = IF[(cell_mask_distance_map <= fraction) & (cell_mask_distance_map > 0)]
            IF_summed = IF[nucleus_mask == 0].sum()
            IF_summed_periph = IF_periph.sum()
            periph_frac = float(IF_summed_periph) / float(IF_summed)
            arr.append(periph_frac)
    return arr



def compute_degree_of_clustering(image_list,file_handler,mtoc_file_handler):
    h_star_l=[]
    for image in image_list:
        h_star = image_descriptors.get_h_star(file_handler, image)
        d=np.array(h_star[h_star>1]-1).sum()
        mtoc_quad=image_descriptors.get_mtoc_quad(mtoc_file_handler,image)
        if mtoc_quad == 1.0:
            h_star_l.append(d)
    return h_star_l


def compute_cytoplasmic_spread_2D(file_handler, image_list,path_data):
    cyt_spreads=[]
    for image in image_list:
        print(image)
        if 'mrna' in image:
            cyt_spread=image_descriptors.compute_mrna_cytoplasmic_spread_2D(file_handler,image)
        else:
            cyt_spread=image_descriptors.compute_protein_cytoplasmic_spread_2D(file_handler,image,path_data)
        cyt_spreads.append(cyt_spread)
    return cyt_spreads


def compute_cytoplasmic_spread(image_list, file_handler, path_data):
    cyt_spreads=[]
    for image in image_list:
        if 'mrna' in image:
            cyt_spread=image_descriptors.compute_mrna_cytoplasmic_spread(file_handler,image)
        else:
            cyt_spread=image_descriptors.compute_protein_cytoplasmic_spread(file_handler,image,path_data)
        cyt_spreads.append(cyt_spread)
    return cyt_spreads


def compute_cytoplasmic_total(image_list,file_handler, path_data):
    total_cyts=[]
    for image in image_list:
        print(image)
        if 'mrna' in image:
            total_cyt=image_descriptors.compute_mrna_cytoplasmic_total(file_handler,image)
        else:
            total_cyt=image_descriptors.compute_protein_cytoplasmic_total(file_handler,image,path_data)
        total_cyts.append(total_cyt)
    return total_cyts

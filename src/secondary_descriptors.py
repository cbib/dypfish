#!/usr/bin/python
# encoding: UTF-8
import io
import signal
import os

import h5py
from numpy import matlib

import matplotlib
try:
    import Tkinter
except:
    matplotlib.use('Agg')
    # required on headless Linux servers
import matplotlib.pyplot as plt

import path, logger, constants
import helpers as helps
import numpy as np

import math
from image_descriptors import *

def compute_cell_mask_3d(
    basic_h5_file_handler,
    image
    ):
   height_map = get_height_map(basic_h5_file_handler, image)
   if height_map.size == 0:
       return None

   zero_level = get_zero_level(basic_h5_file_handler, image)
   cell_mask = get_cell_mask(basic_h5_file_handler, image)
   height_map = height_map.astype(float)
   height_map[cell_mask == 0] = np.nan
   reversed_height_map = zero_level-height_map+1

   # Create binary cell masks per slice

   cell_masks = np.zeros((512,512,zero_level))

   # build slice mask
   for slice in range(0,zero_level):
       slice_mask = np.array(reversed_height_map, copy = True)
       slice_mask[slice_mask > slice] = np.nan
       slice_mask[slice_mask <= float(slice)] = 1
       slice_mask=np.nan_to_num(slice_mask)
       cell_masks[:,:,slice]


def set_h_star_mrna(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    image
    ):

    H_star = get_h_star(basic_h5_file_handler, image)
    assert (not H_star.size), 'H_star already defined for %r' % image

    # get descriptors
    spots = get_spots(basic_h5_file_handler, image)

    cell_mask_3d = compute_cell_mask_3d(basic_h5_file_handler, image)
    if cell_mask_3d is None:
        print("[set_h_star_mrna] : Skipping hstar")
        return False
    #h_star=helps.clustering_index_point_process_2d(spots, cell_mask_3d)
    h_star=helps.clustering_index_point_process(spots, cell_mask_3d)
    descriptor = image + '/h_star'
    #secondary_h5_file_handler[descriptor][:] = h_star
    secondary_h5_file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)


def set_h_star_protein(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    raw_images_dir_path,
    image
    ):

    def prune_intensities(image, zero_level):
        molecule, gene, timepoint, number = image.split("/")
        IF_image_path = raw_images_dir_path + gene + '/' + molecule + '_' + timepoint + "/image_" + number + '/IF.tif'
        IF = io.imread(IF_image_path, plugin='tifffile')

        vol_block = np.zeros((512, 512, zero_level))
        for c_slice in range(0, zero_level):
            vol_block[:, :, c_slice] = IF[c_slice, :, :]

        return vol_block  # .reshape((512,512,zero_level + constants.VOLUME_OFFSET))

    cell_mask_3d=compute_cell_mask_3d(basic_h5_file_handler,image)
    if cell_mask_3d is None:
        print("[set_h_star_mrna] : Skipping hstar")
        return False
    zero_level=get_zero_level(basic_h5_file_handler,image)
    IF = prune_intensities(image, zero_level)
    IF = IF * cell_mask_3d
    h_star = helps.clustering_index_random_measure(IF, cell_mask_3d)
    descriptor = image + '/h_star'
    # secondary_h5_file_handler[descriptor][:] = h_star
    secondary_h5_file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)


def set_zero_level(
    basic_h5_file_handler,
    image,
    raw_images_dir_path
    ):
    molecule, gene, timepoint, number = image.split("/")
    tubulin_image_path = raw_images_dir_path+gene+'/'+molecule+'_'+timepoint+"/image_"+number+"/tubulin.tif"
    if not os.path.exists(tubulin_image_path):
        tubulin_image_path = raw_images_dir_path + gene + '/' + molecule + '_' + timepoint +"/"+ number + "/tubulin.tif"

    image_stacked = io.imread(tubulin_image_path, plugin='tifffile')
    #print(image_stacked.shape)
    '''
    if len(image_stacked.shape)  == 2: #only one stack of image
        x_size, y_size = image_stacked.shape
        z_size = 1
    else:
        z_size, x_size, y_size = image_stacked.shape
    '''
    if len(image_stacked.shape) == 2:  # TODO find a way to make things generic here by adapting to both 2d and 3d cases
        print ("This image is not a 3D one, removing it from basic.h5 to prevent further compatibility issues")
        del basic_h5_file_handler[image]
        return
    z_size, x_size, y_size = image_stacked.shape

    #print(z_size, x_size, y_size)
    sums = np.zeros((z_size, 1))
    for n in range(z_size):
        slice = image_stacked[n, :, :] if z_size > 1 else  image_stacked[ :, :]
        sums[n] = slice.sum()
    basic_h5_file_handler[image].attrs['zero_level'] = sums.argmax()


def set_3d_spots(
    basic_h5_file_handler,
    image
    ):
    if basic_h5_file_handler[image]['spots'].shape[1] == 3:
        return
    spots_3d = []
    for spot in basic_h5_file_handler[image]['spots']:
        spots_3d.append((spot[0],spot[1], basic_h5_file_handler[image]['height_map'][spot[0],spot[1]]))

    del basic_h5_file_handler[image]['spots']
    basic_h5_file_handler.create_dataset(image+'/spots', data=spots_3d, dtype=np.float32)



'''Cell mask distance map computation'''
# Given the nucleus and cytoplasm mask and x, y slopes of a line, compute the segment
# of the line that falls withing the nucleus and the segment that falls within the cytoplasm
def compute_line_segments(
    nucleus_mask,
    cytoplasm_mask,
    nucleus_centroid,
    x_slope,
    y_slope
    ):
    logger.info('Compute line segments')
    cytoplasm_segment = np.zeros(constants.MAX_CELL_RADIUS)
    nucleus_segment = np.zeros(constants.MAX_CELL_RADIUS)

    # check for each point of the line (length MAX_CELL_RADIUS) if it falls within the nucleus or the cytoplasm
    #print(x_slope)
    #print(y_slope)

    for point in range(constants.MAX_CELL_RADIUS):
        x = int(round(nucleus_centroid[0] + point * x_slope))
        #print(x)
        y = int(round(nucleus_centroid[1] + point * y_slope))
        #print(y)
        if not (x < 0 or x > 511 or y < 0 or y > 511):
            cytoplasm_segment[point] = cytoplasm_mask[y, x]
            nucleus_segment[point] = nucleus_mask[y, x]
        else:
            logger.info('Warning: spot coordinates out of image bounds', '(', x, ',', y, ')')
    #sys.exit()
    return nucleus_segment, cytoplasm_segment


# Given line segments, return points that fall on the nucleus and on the cytoplasm edges
def compute_edge_points(
    nucleus_segment,
    cytoplasm_segment
    ):
    nucleus_points = np.where(nucleus_segment == 1)
    if len(nucleus_points[0]) > 0:
        nucleus_edge_point = nucleus_points[0][-1]
    else:
        nucleus_edge_point = 1

    cytoplasm_points = np.where(cytoplasm_segment == 1)
    if len(cytoplasm_points[0]) > 0:
        cytoplasm_edge_point = cytoplasm_points[0][-1]
    else:
        cytoplasm_edge_point = 1
    return nucleus_edge_point, cytoplasm_edge_point


def compute_cell_mask_distance_map(
    nucleus_mask,
    cytoplasm_mask,
    contour_points
    ):
    cell_mask_distance_map = np.zeros((512, 512), dtype=np.int)
    for index in range(constants.NUM_CONTOURS):
        if index == 0:
            peripheral_mask = nucleus_mask
        else:
            contour_num = constants.NUM_CONTOURS - index
            peripheral_mask = helps.create_mask(contour_points[:, contour_num, 1], contour_points[:, contour_num, 0],
                                                (512, 512))
            peripheral_mask &= cytoplasm_mask
        cell_mask_distance_map[(peripheral_mask == 1)] = index +1

    cell_mask_distance_map[(cytoplasm_mask == 0)] = 0
    return cell_mask_distance_map


# Create a map with concentric isolines subdividing the cytoplasm into slices NUM_CONTOURS
# Each isoline is built by constructing a polygone from 360 points (one point per degree)
def set_cell_mask_distance_map(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    image
    ):

    #descriptor = image + '/cell_mask_distance_map'
    #assert (descriptor not in secondary_h5_file_handler.keys(), 'cell_mask_distance_map already defined for %r' % image)

    cell_mask_distance_map = get_cell_mask_distance_map(basic_h5_file_handler, image)
    assert (not cell_mask_distance_map.size), 'cell_mask_distance_map already defined for %r' % image
    cell_mask = get_cell_mask(basic_h5_file_handler, image).astype(int)
    nucleus_mask = get_nucleus_mask(basic_h5_file_handler, image).astype(int)
    nucleus_centroid = get_nucleus_centroid(basic_h5_file_handler, image).transpose()
    contour_points = np.zeros((360, 100, 2))
    cytoplasm_mask = (cell_mask == 1) & (nucleus_mask == 0)

    #print(nucleus_centroid)

    # for each degree, analyse the line segment between the nucleus and the periphery
    for degree in range(360):
        angle = degree * 2 * math.pi / 360
        x_slope, y_slope = math.sin(angle), math.cos(angle)
        nucleus_segment, cytoplasm_segment = compute_line_segments(nucleus_mask, cytoplasm_mask, nucleus_centroid,
                                                                   x_slope, y_slope)
        nucleus_edge_point, cytoplasm_edge_point = compute_edge_points(nucleus_segment, cytoplasm_segment)
        segment_length = cytoplasm_edge_point - nucleus_edge_point


        for index in range(constants.NUM_CONTOURS):
            point = nucleus_edge_point + segment_length * index / constants.NUM_CONTOURS
            x = int(round(nucleus_centroid[0] + point * x_slope))
            y = int(round(nucleus_centroid[1] + point * y_slope))
            contour_points[degree, index, :] = [x, y]

    # end = time.time()
    # print(end - start)
    # sys.exit()
    cell_mask_distance_map = compute_cell_mask_distance_map(nucleus_mask, cytoplasm_mask, contour_points)
    # cell_mask_distance_map[(cell_mask == 1) & (cell_mask_distance_map == 0)] = 1
    # cell_mask_distance_map[nucleus_mask == 1] = 0

    # plt.imshow(cell_mask_distance_map, cmap='hot')
    # plt.show()
    descriptor = image + '/cell_mask_distance_map'
    secondary_h5_file_handler.create_dataset(descriptor, data=cell_mask_distance_map, dtype=np.int)


'''Cell and nucleus area computation'''
def set_cell_area(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    image
    ):
    """compute cell surface by pixel using cell mask"""
    cell_area = get_cell_area(secondary_h5_file_handler, image)
    assert (cell_area == -1), 'cell_area already defined for %r' % image

    #descriptor = image + '/cell_area'
    #assert (descriptor not in secondary_h5_file_handler.attrs.keys()), 'cell_area already defined for %r' % image

    cell_mask = get_cell_mask(basic_h5_file_handler, image)
    area = cell_mask.sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)  # * by pixel dimensions
    secondary_h5_file_handler[image].attrs['cell_area'] = area

def set_nucleus_area(
    basic_h5_file_handler,
    secondary_h5_file_handler,
    image
    ):
    """compute nucleus surface in pixel using nucleus mask"""
    nucleus_area = get_nucleus_area(secondary_h5_file_handler, image)
    assert (nucleus_area == -1), 'nucleus_area already defined for %r' % image
    nucleus_mask = get_nucleus_mask(basic_h5_file_handler, image)
    area = nucleus_mask.sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)  # * by pixel dimensions
    secondary_h5_file_handler[image].attrs['nucleus_area'] = area


from threading import Thread

class Preprocess(Thread):

    def __init__(
        self,
        basic_h5_file_handler,
        secondary_h5_file_handler,
        raw_images_dir_path,
        img_list,
        ext_logger=None
        ):
        Thread.__init__(self)
        self.basic_h5_file_handler = basic_h5_file_handler
        self.secondary_h5_file_handler = secondary_h5_file_handler
        self.raw_images_dir_path = raw_images_dir_path
        self.img_list = img_list
        self.ext_logger = ext_logger

    def run(self):
        if self.ext_logger:
            self.ext_logger("Running preprocess...")
        try:
            for image in self.img_list:
                print("Computing descriptors for " + image)

                '''UPDATE BASIC.H5'''
                set_zero_level(
                    self.basic_h5_file_handler,
                    image,
                    self.raw_images_dir_path
                    )

                '''PREPROCESS'''
                set_cell_mask_distance_map(
                    self.basic_h5_file_handler,
                    self.secondary_h5_file_handler,
                    image
                    )
                if 'mrna' in image:
                    set_3d_spots(
                        self.basic_h5_file_handler,
                        image
                        )
                    set_spots_peripheral_distance(
                        self.basic_h5_file_handler,
                        self.secondary_h5_file_handler,
                        image
                        )
                    set_spots_peripheral_distance_2D(
                        self.basic_h5_file_handler,
                        self.secondary_h5_file_handler,
                        image
                        )
                set_cell_area(
                    self.basic_h5_file_handler,
                    self.secondary_h5_file_handler,
                    image
                    )
                set_nucleus_area(
                    self.basic_h5_file_handler,
                    self.secondary_h5_file_handler,
                    image
                    )
                '''PREPROCESS'''

                if 'mrna' in image:
                    set_h_star_mrna(
                        self.basic_h5_file_handler,
                        self.secondary_h5_file_handler,
                        image
                        )
                elif 'protein' in image:
                    set_h_star_protein(
                        self.basic_h5_file_handler,
                        self.secondary_h5_file_handler,
                        self.raw_images_dir_path,
                        image
                        )
        except Exception as e:
            if self.ext_logger:
                self.ext_logger.exception(e)


def generate_secondary_descriptors_hd5(
    basic_h5_file_path_name,
    raw_images_dir_path,
    secondary_h5_file_path_name,
    thread_num=4,
    ext_logger=None
    ):

    if ext_logger:
        ext_logger.debug("Running generate_secondary_descriptors_hd5...")
        ext_logger.debug("basic_h5_file_path_name : %s" % basic_h5_file_path_name)
        ext_logger.debug("raw_images_dir_path : %s" % raw_images_dir_path)
        ext_logger.debug("secondary_h5_file_path_name : %s" % secondary_h5_file_path_name)
        ext_logger.debug("thread_num : %s" % thread_num)
        ext_logger.debug("ext_logger : %s" % ext_logger)

    assert os.path.isfile(basic_h5_file_path_name)
    assert os.path.isdir(raw_images_dir_path)
    assert not os.path.isfile(secondary_h5_file_path_name)

    ## Start géneration

    if ext_logger:
        ext_logger.debug("Start géneration...")

    with h5py.File(basic_h5_file_path_name, "a") as basic_h5_file_handler, \
         h5py.File(secondary_h5_file_path_name, "a") as secondary_h5_file_handler:

        for molecule_type in (['mrna'], ['protein']):
            thread_list = []
            image_list = helps.preprocess_image_list(
                basic_h5_file_handler,
                molecule_type
                )
            for sub_list in np.array_split(
                image_list,
                thread_num
                ):
                if ext_logger:
                    ext_logger.debug("---> len(sub_list): %s" % len(sub_list))
                thread_list.append(
                    Preprocess(
                        basic_h5_file_handler,
                        secondary_h5_file_handler,
                        raw_images_dir_path,
                        sub_list,
                        ext_logger=ext_logger
                        )
                    )

            for thread in thread_list:
                thread.start()

            for thread in thread_list:
                thread.join()

    if ext_logger:
        ext_logger.debug("end of generate_secondary_descriptors_hd5 function")

    return True


if __name__ == "__main__":

    basic_h5_file_path_name = path.analysis_data_dir + 'basic.h5'
    raw_images_dir_path = path.raw_data_dir

    secondary_h5_file_path_name = path.analysis_data_dir + 'secondary.h5'
    if os.path.exists(secondary_h5_file_path_name):
        os.remove(secondary_h5_file_path_name)

    THREAD_NUM = 4

    generate_secondary_descriptors_hd5(
        basic_h5_file_path_name=basic_h5_file_path_name,
        raw_images_dir_path=raw_images_dir_path,
        secondary_h5_file_path_name=secondary_h5_file_path_name,
        thread_num=THREAD_NUM
        )

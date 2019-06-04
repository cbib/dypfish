#!/usr/bin/env python2
# encoding: UTF-8

from __future__ import print_function
from __future__ import print_function
from helpers import *
import logging
import sys

logger = logging.getLogger('DYPFISH IMAGE DESCRIPTORS')

####### Primary descriptors ########
def get_spots(file_handler, image):
    descriptor = image + '/spots'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    spots = np.around(np.array(file_handler[descriptor])).astype(np.int)
    return spots

def get_IF(file_handler, image):
    descriptor = image + '/IF'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    IF = np.array(file_handler[descriptor])
    return IF

def get_nucleus_mask(file_handler, image):
    descriptor = image + '/nucleus_mask'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    nucleus_mask = np.array(file_handler[descriptor])
    return nucleus_mask

def get_z_lines_masks(file_handler, image):
    descriptor = image + '/z_lines'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    z_lines=[]
    for zmask_n in range(1,len(file_handler[descriptor])+1):
        z_linemask=np.array(file_handler[descriptor+'/'+str(zmask_n)]).astype(int)
        z_lines.append(z_linemask)
    return z_lines

def get_cell_mask(file_handler, image):
    descriptor = image + '/cell_mask'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    cell_mask = np.array(file_handler[descriptor]).astype(int)
    return cell_mask

def get_z_lines_mask(file_handler, image):
    descriptor = image + '/z_lines_mask'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    z_lines_mask = np.array(file_handler[descriptor]).astype(int)
    return z_lines_mask

# Returns integer coordinates of the nucleus
def get_multiple_nucleus_centroid(file_handler, image):
    nucleus_centroid_l=[]
    for dataset in file_handler[image]:
        if 'nucleus_centroid' in dataset:
             descriptor = image + '/'+dataset
             nucleus_centroid = np.around(file_handler[descriptor]).flatten().astype(np.int)
             nucleus_centroid_l.append(nucleus_centroid)
    return nucleus_centroid_l

# Returns integer coordinates of the nucleus
def get_nucleus_centroid(file_handler, image):
    descriptor = image + '/nucleus_centroid'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    nucleus_centroid = np.around(file_handler[descriptor]).flatten().astype(np.int)
    return nucleus_centroid

# Returns integer coordinates of the MTOC
def get_mtoc_position(file_handler, image):
    descriptor = image + '/mtoc_position'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    mtoc_position = np.around(file_handler[descriptor]).flatten().astype(np.int)
    return mtoc_position

####### Secondary descriptors #######
def get_h_star(file_handler, image):
    descriptor = image + '/h_star'

    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    h_star = np.array(file_handler[descriptor])
    return h_star

def set_h_star_protein(file_handler, output_file_handler, image):

    cell_mask_3d=compute_cell_mask_3d(file_handler,image)
    zero_level=get_zero_level(file_handler,image)
    IF = prune_intensities(image, zero_level)
    IF = IF * cell_mask_3d
    h_star = clustering_index_random_measure(IF, cell_mask_3d)
    descriptor = image + '/h_star'
    output_file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)


def compute_cell_mask_3d(file_handler,image):
   height_map = get_height_map(file_handler, image)
   zero_level = get_zero_level(file_handler, image)
   cell_mask = get_cell_mask(file_handler, image)
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
       cell_masks[:,:,slice] = slice_mask
   return cell_masks


"""compute h_star"""

def get_zero_level(file_handler, image):
    attributes = file_handler[image].attrs.keys()
    if 'zero_level' not in attributes:
        logging.critical('Zero level attribute not computed')
        exit(1)
        return -1
    zero_level = file_handler[image].attrs['zero_level']
    return zero_level

def set_h_star_mrna(file_handler, output_file_handler, image):
    H_star = get_h_star(file_handler, image)
    assert (not H_star.size), 'H_star already defined for %r' % image
    spots = get_spots(file_handler, image)
    cell_mask_3d=compute_cell_mask_3d(file_handler,image)
    h_star=clustering_index_point_process(spots, cell_mask_3d)
    descriptor = image + '/h_star'
    output_file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)

def get_peripheral_mask(file_handler, image):
    descriptor = image + '/peripheral_mask'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    peripheral_mask = np.array(file_handler[descriptor])
    return peripheral_mask

def set_peripheral_mask(file_handler, image):
    """compute nucleus surface in pixel using nucleus mask"""
    peripheral_mask = get_peripheral_mask(file_handler, image)
    assert (not peripheral_mask.size), 'peripheral_mask already defined for %r' % image
    cell_mask_distance_map = get_cell_mask_distance_map(file_handler, image)

def get_cell_mask_distance_map(file_handler, image):
    descriptor = image + '/cell_mask_distance_map'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    cell_mask_distance_map = np.array(file_handler[descriptor])
    return cell_mask_distance_map


# Given the nucleus and cytoplasm mask and x, y slopes of a line, compute the segment
# of the line that falls withing the nucleus and the segment that falls within the cytoplasm
def compute_line_segments(nucleus_mask, cytoplasm_mask, nucleus_centroid, x_slope, y_slope):
    logger.info('Compute line segments')
    cytoplasm_segment = np.zeros(constants.MAX_CELL_RADIUS)
    nucleus_segment = np.zeros(constants.MAX_CELL_RADIUS)

    # check for each point of the line (length MAX_CELL_RADIUS) if it falls within the nucleus or the cytoplasm
    for point in range(constants.MAX_CELL_RADIUS):
        x = int(round(nucleus_centroid[0] + point * x_slope))
        y = int(round(nucleus_centroid[1] + point * y_slope))
        if not (x < 0 or x > 511 or y < 0 or y > 511):
            cytoplasm_segment[point] = cytoplasm_mask[y, x]
            nucleus_segment[point] = nucleus_mask[y, x]
        else:
            logger.info('Warning: spot coordinates out of image bounds', '(', x, ',', y, ')')
    return nucleus_segment, cytoplasm_segment


# Given line segments, return points that fall on the nucleus and on the cytoplasm edges
def compute_edge_points(nucleus_segment, cytoplasm_segment):
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


'''
def f(x, y, nuc_centroid):
    return x*y + nuc_centroid
'''


def set_height_map(file_handler, image, tubulin_image_path):
    cell_mask = get_cell_mask(file_handler, image)
    zero_level = get_zero_level(file_handler, image)
    spots=get_spots(file_handler,image)

    # get tubulin stacked image
    image3D = io.imread(tubulin_image_path, plugin='tifffile')

    if "gapdh" in image:
        percentile = 7
    elif "beta_actin" in image:
        percentile = 6
    elif "cultured/1h" in image:
        percentile = 2
    elif "cultured/3h" in image:
        percentile = 1
    elif "arhgdia" in image:
        percentile = 9
    elif "pard3" in image:
        percentile = 9
    elif "pkp4" in image:
        percentile = 9
    elif "nocodazole" in image:
        percentile = 2
    else:
        percentile = 9

    height_map = np.zeros((512, 512))
    test_counter=0
    for n in range(zero_level, -1, -1):
        fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2,figsize=(20,10))
        slice = np.array(image3D[n, :, :] , copy=True)
        ax1.imshow(slice, cmap='gray')
        bins_percent=[1,10,20,30,40,50,55,60,70,80,85,90,95,100]
        percent = np.percentile(slice, bins_percent)
        threshold=percent[percentile]
        if "beta_actin" in image or "cultured" in image:
            threshold=percent[percentile]+test_counter
        elif "nocodazole" in image:
            threshold=percent[percentile]+24
        else:
            threshold = percent[percentile]
        slice[slice <= threshold] = 0
        slice[slice > threshold] = 1
        slice[cell_mask == 0] = 0
        ax2.imshow(slice, cmap='gray')
        from scipy import ndimage as ndi
        from skimage.restoration import denoise_nl_means
        slice = ndi.binary_fill_holes(slice)
        slice = denoise_nl_means(slice, 5, 6, 0.4, multichannel=True)
        slice[slice < 1] = 0
        height_map[cell_mask == 0] = 0
        if n == zero_level:
            height_map[cell_mask==1]=1
        height_map[slice > 0] = (zero_level + 1) - n
        ax4.imshow(height_map, cmap='hot')
        test_counter+=4
    plt.imshow(height_map)
    plt.show()


def get_peripheral_distance(file_handler, image):
    descriptor = image + '/peripheral_distance'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    peripheral_distance = np.array(file_handler[descriptor])
    return peripheral_distance

def get_spots_peripheral_distance(file_handler, image):
    descriptor = image + '/spots_peripheral_distance'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    spots_peripheral_distance = np.array(file_handler[descriptor])
    return spots_peripheral_distance

def get_spots_peripheral_distance_2D(file_handler, image):
    descriptor = image + '/spots_peripheral_distance_2D'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    spots_peripheral_distance_2D = np.array(file_handler[descriptor])
    return spots_peripheral_distance_2D

def set_scratch_peripheral_distance(file_handler, output_file_handler, image):
    spots_peripheral_distance = get_spots_peripheral_distance(output_file_handler, image)
    assert (not spots_peripheral_distance.size), 'spots_peripheral_distance already defined for %r' % image
    spots = get_spots(file_handler, image)
    if len(spots) == 0:
        return
    periph_dist_map = get_cell_mask_distance_map(output_file_handler, image)
    spots_peripheral_distance = []
    for i in range(len(spots)):
        if is_in_cytoplasm(file_handler, image, [spots[i,1], spots[i,0]]):
            spots_peripheral_distance.append(periph_dist_map[spots[i,1],spots[i,0]])
    descriptor = image + '/spots_peripheral_distance'
    output_file_handler.create_dataset(descriptor, data=spots_peripheral_distance, dtype=np.uint8)

def set_protein_peripheral_distance(file_handler, output_file_handler, image):
    peripheral_distance = get_peripheral_distance(output_file_handler, image)
    assert (not peripheral_distance.size), 'peripheral_distance already defined for %r' % image
    height_map = get_height_map(file_handler, image)
    zero_level = get_zero_level(file_handler, image)
    periph_dist_map = get_cell_mask_distance_map(output_file_handler, image)
    for n in range(zero_level, -1, -1):
        height_map_copy = np.array(height_map, copy=True)
        height_map_copy[height_map >= zero_level + 1 - n] = 1
        slice_area = height_map_copy[height_map_copy == 1].sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)
        isoline_area = compute_isoline_area(file_handler, output_file_handler, image)
        test_index = find_nearest(isoline_area, slice_area)
        prot_in_slice = spots[np.around(spots[:, 2]) == n]


'''spot peripheral distance 2D et 3D computation'''
def set_spots_peripheral_distance_2D(file_handler, output_file_handler, image):
    spots_peripheral_distance_2D = get_spots_peripheral_distance_2D(output_file_handler, image)
    assert (not spots_peripheral_distance_2D.size), 'spots_peripheral_distance 2D already defined for %r' % image
    spots = get_spots(file_handler, image)
    if len(spots) == 0:
        return
    nucleus_mask = get_nucleus_mask(file_handler, image)

    cell_mask = get_cell_mask(file_handler, image)
    periph_dist_map = get_cell_mask_distance_map(output_file_handler, image)
    spots_peripheral_distance_2D = []
    for spot in spots:

        if nucleus_mask[spot[1],spot[0]]==0:
            if cell_mask[spot[1],spot[0]]==1 and periph_dist_map[spot[1],spot[0]] ==0:
                spots_peripheral_distance_2D.append(int(1))
            else:
                spots_peripheral_distance_2D.append(periph_dist_map[spot[1], spot[0]])

    descriptor = image + '/spots_peripheral_distance_2D'
    output_file_handler.create_dataset(descriptor, data=spots_peripheral_distance_2D, dtype=np.uint8)


def set_spots_peripheral_distance(file_handler, output_file_handler, image):
    spots_peripheral_distance = get_spots_peripheral_distance(output_file_handler, image)
    assert (not spots_peripheral_distance.size), 'spots_peripheral_distance already defined for %r' % image
    spots = get_spots(file_handler, image)
    if len(spots)==0:
        return
    height_map = get_height_map(file_handler, image)
    zero_level = get_zero_level(file_handler, image)
    periph_dist_map = get_cell_mask_distance_map(output_file_handler, image)
    spots_peripheral_distance = []
    for n in range(zero_level, -1, -1):
        height_map_copy = np.array(height_map, copy=True)
        height_map_copy[height_map >= zero_level + 1 - n] = 1
        slice_area = height_map_copy[height_map_copy == 1].sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)
        isoline_area = compute_isoline_area(file_handler, output_file_handler, image)
        test_index = find_nearest(isoline_area, slice_area)
        spots_in_slice = spots[np.around(spots[:, 2]) == n]
        for j in range(len(spots_in_slice)):
            if is_in_cytoplasm(file_handler, image, [spots_in_slice[j][1], spots_in_slice[j][0]]):
                old_periph_distance = periph_dist_map[spots_in_slice[j][1], spots_in_slice[j][0]]
                if old_periph_distance < test_index:
                    old_periph_distance = test_index
                if old_periph_distance == 100:
                    test_index = 0
                max_periph_tubulin = 100 - test_index
                new_periph_distance = (old_periph_distance - (100 - max_periph_tubulin)) * (100.0 / max_periph_tubulin)
                if np.around(new_periph_distance) == 0.0:
                    spots_peripheral_distance.append(int(1))
                else:
                    spots_peripheral_distance.append(int(np.around(new_periph_distance)))

    descriptor = image + '/spots_peripheral_distance'
    output_file_handler.create_dataset(descriptor, data=spots_peripheral_distance, dtype=np.uint8)

def get_spots_peripheral_flags(file_handler, image):
    descriptor = image + '/spots_peripheral_flags'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    spots_peripheral_flags = np.array(file_handler[descriptor])
    return spots_peripheral_flags

def get_spots_cytoplasm_flags(file_handler, image):
    descriptor = image + '/spots_cytoplasm_flags'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    spots_cytoplasm_flags = np.array(file_handler[descriptor])
    return spots_cytoplasm_flags


def set_spots_peripheral_flags(file_handler, output_file_handler, image):
    spots = get_spots(file_handler, image)
    zero_level = get_zero_level(file_handler, image)
    num_spots = len(spots)
    spots_peripheral_distance = get_spots_peripheral_distance(output_file_handler, image)
    spots_peripheral_flags = np.zeros(num_spots, dtype=np.int)
    counter=0
    for n in range(zero_level, -1, -1):
        spots_in_slice = spots[np.around(spots[:, 2]) == n]
        for j in range(len(spots_in_slice)):
            if is_in_cytoplasm(file_handler, image, [spots_in_slice[j][1], spots_in_slice[j][0]]):
                if spots_peripheral_distance[counter] > constants.PERIPHERAL_FRACTION_THRESHOLD:
                    spots_peripheral_flags[counter] = 0
                else:
                    spots_peripheral_flags[counter] = 1
                counter += 1
    descriptor = image + '/spots_peripheral_flags'
    output_file_handler.create_dataset(descriptor, data=spots_peripheral_flags, dtype=np.int)

def get_cell_area(file_handler, image):
    descriptor = image + '/cell_area'
    attributes = file_handler[image].attrs.keys()
    if 'cell_area' not in attributes:
        return -1
    return file_handler[image].attrs['cell_area']

def set_cell_area(file_handler, output_file_handler, image):
    """compute cell surface by pixel using cell mask"""
    cell_area = get_cell_area(output_file_handler, image)
    assert (cell_area == -1), 'cell_area already defined for %r' % image
    cell_mask = get_cell_mask(file_handler, image)
    area = cell_mask.sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)
    output_file_handler[image].attrs['cell_area'] = area


def get_nucleus_area(file_handler, image):
    """get nucleus_area from image descriptors file"""
    attributes = file_handler[image].attrs.keys()
    if 'nucleus_area' not in attributes:
        return -1
    nucleus_area = file_handler[image].attrs['nucleus_area']
    return nucleus_area

def set_nucleus_area(file_handler, output_file_handler, image):
    """compute nucleus surface in pixel using nucleus mask"""
    nucleus_area = get_nucleus_area(output_file_handler, image)
    assert (nucleus_area == -1), 'nucleus_area already defined for %r' % image
    nucleus_mask = get_nucleus_mask(file_handler, image)
    area = nucleus_mask.sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2)
    output_file_handler[image].attrs['nucleus_area'] = area

# return the mask for the quandrants
def get_quadrants(file_handler, image):
    descriptor = image + '/quadrants'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    quadrants = np.array(file_handler[descriptor])
    return quadrants

# Calculates the quadrant mask for the MTOC
def set_quadrants(file_handler, image):
    quadrants = get_quadrants(file_handler, image)
    assert (not quadrants.size), 'mtoc_quadrant already defined for %r' % image
    mtoc_position = get_mtoc_position(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)

    # the quadrant of MTOC is defined by two lines 45 degrees to the right
    right_point = rotate_point(nucleus_centroid, mtoc_position, 45)
    s = slope_from_points(nucleus_centroid, right_point)
    corr = np.arctan(s)
    xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
    rotated_xx, rotated_yy = rotate_meshgrid(xx, yy, -corr)
    sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
    sliceno = sliceno.astype(int)
    quadrant_mask = sliceno + cell_mask
    quadrant_mask[quadrant_mask == 5]=4
    quadrant_mask[cell_mask == 0] = 0
    return quadrant_mask

# Calculates the quadrant mask for the MTOC
def search_protein_quadrants(file_handler, mtoc_file_handler,protein,image):
    print(image)
    quadrants = get_quadrants(file_handler, image)
    assert (not quadrants.size), 'mtoc_quadrant already defined for %r' % image
    mtoc_position = get_mtoc_position(file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    IF = get_IF(file_handler,image)
    IF[cell_mask == 0] = 0
    IF[nucleus_mask == 1] = 0
    intensity_by_quad = np.zeros((90, 4, 2))
    mtoc_quads = []
    height_map = height_map.astype(float)
    height_map[(cell_mask == 1) & (height_map == 0)] = 0.5
    height_map[cell_mask == 0] = 0
    for i in range(90):
        # the quadrant of MTOC is defined by two lines 45 degrees to the right
        right_point = rotate_point(nucleus_centroid, mtoc_position, i)
        s = slope_from_points(nucleus_centroid, right_point)
        corr = np.arctan(s)  # angle wrt to x axis
        xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0],
                             np.array(range(0, 512)) - nucleus_centroid[1])
        rotated_xx, rotated_yy = rotate_meshgrid(xx, yy, -corr)
        sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == 5] = 4
        quadrant_mask[cell_mask == 0] = 0
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        mtoc_quads.append(mtoc_quad)
        for quad_num in range(1, 5):
            intensity_by_quad[i, quad_num - 1, 0] = np.sum(IF[quadrant_mask == quad_num])
            if quad_num == mtoc_quad:
                intensity_by_quad[i, mtoc_quad - 1, 1] = 1
            intensity_by_quad[i, quad_num - 1, 0] /= np.sum(height_map[quadrant_mask == quad_num]) * constants.VOLUME_COEFFICIENT
    return intensity_by_quad

# Calculates the quadrant mask for the MTOC
def search_mrna_quadrants(file_handler, image):
    mtoc_position = get_mtoc_position(file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    spots = get_spots(file_handler, image)
    spot_by_quad = np.zeros((90, 4, 2))
    height_map = height_map.astype(float)
    height_map[(cell_mask == 1) & (height_map == 0)]=0.5
    height_map[cell_mask == 0] = 0
    height_map[nucleus_mask == 1] = 0

    for i in range(90):
        # the quadrant of MTOC is defined by two lines 45 degrees to the right
        right_point = rotate_point(nucleus_centroid, mtoc_position, i)
        s = slope_from_points(nucleus_centroid, right_point)
        corr = np.arctan(s) # angle wrt to x axis
        xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
        rotated_xx, rotated_yy = rotate_meshgrid(xx, yy, -corr)
        sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == 5]=4
        quadrant_mask[cell_mask == 0] = 0
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        # assign each spot to the corresponding quadrant
        for spot in spots:
            if nucleus_mask[spot[1], spot[0]]==1:
                continue
            spot_quad = quadrant_mask[spot[1], spot[0]]
            spot_by_quad[i, spot_quad - 1, 0] += 1
        # mark leading edge quadrant
        spot_by_quad[i, mtoc_quad - 1, 1] = 1
        for quad_num in range(1,5):
            spot_by_quad[i, quad_num - 1, 0] /= np.sum(height_map[quadrant_mask == quad_num]) * constants.VOLUME_COEFFICIENT
    return spot_by_quad




# Calculates the quadrant mask for the MTOC
def search_periph_mrna_quadrants(file_handler,second_file_handler, image):
    print(image)
    mtoc_position = get_mtoc_position(file_handler, image)
    cell_mask_dist_map = get_cell_mask_distance_map(second_file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)
    peripheral_binary_mask =(cell_mask_dist_map > 0) & (cell_mask_dist_map <= constants.PERIPHERAL_FRACTION_THRESHOLD).astype(int)
    spots = get_spots(file_handler, image)
    spot_by_quad = np.zeros((90, 4, 2))
    height_map = height_map.astype(float)
    height_map[(cell_mask == 1) & (height_map == 0)] = 0.5
    height_map[cell_mask == 0] = 0
    for i in range(90):
        # the quadrant of MTOC is defined by two lines 45 degrees to the right
        right_point = rotate_point(nucleus_centroid, mtoc_position, i)
        s = slope_from_points(nucleus_centroid, right_point)
        corr = np.arctan(s)
        xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0], np.array(range(0, 512)) - nucleus_centroid[1])
        rotated_xx, rotated_yy = rotate_meshgrid(xx, yy, -corr)
        sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == 5]=4
        quadrant_mask[cell_mask == 0] = 0
        quadrant_mask[peripheral_binary_mask==0]=0
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        for spot in spots:
            if peripheral_binary_mask[spot[1], spot[0]]:
                spot_quad = quadrant_mask[spot[1], spot[0]]
                spot_by_quad[i, spot_quad - 1, 0] += 1
        spot_by_quad[i, mtoc_quad - 1, 1] += 1
        for quad_num in range(1,5):
            spot_by_quad[i, quad_num - 1, 0] /= np.sum(height_map[(peripheral_binary_mask==1) & (quadrant_mask == quad_num)]) * constants.VOLUME_COEFFICIENT
    return spot_by_quad

# Calculates the quadrant mask for the MTOC
def search_periph_protein_quadrants(file_handler,second_file_handler, protein,image,path_data):
    print(image)
    quadrants = get_quadrants(file_handler, image)
    assert (not quadrants.size), 'mtoc_quadrant already defined for %r' % image
    mtoc_position = get_mtoc_position(file_handler, image)
    cell_mask_dist_map = get_cell_mask_distance_map(second_file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    IF = get_IF(file_handler,image)
    IF[cell_mask == 0] = 0
    IF[nucleus_mask == 1] = 0
    intensity_by_quad = np.zeros((90, 4, 2))
    height_map = height_map.astype(float)
    height_map[cell_mask == 0] = 0
    height_map[(cell_mask == 1) & (height_map == 0)] = 0.5
    peripheral_binary_mask =(cell_mask_dist_map > 0) & (cell_mask_dist_map <= constants.PERIPHERAL_FRACTION_THRESHOLD).astype(int)
    for i in range(90):
        # the quadrant of MTOC is defined by two lines 45 degrees to the right
        right_point = rotate_point(nucleus_centroid, mtoc_position, i)
        s = slope_from_points(nucleus_centroid, right_point)
        corr = np.arctan(s)  # angle wrt to x axis
        xx, yy = np.meshgrid(np.array(range(0, 512)) - nucleus_centroid[0],
                             np.array(range(0, 512)) - nucleus_centroid[1])
        rotated_xx, rotated_yy = rotate_meshgrid(xx, yy, -corr)
        sliceno = ((math.pi + np.arctan2(rotated_xx, rotated_yy)) * (4 / (2 * math.pi)))
        sliceno = sliceno.astype(int)
        quadrant_mask = sliceno + cell_mask
        quadrant_mask[quadrant_mask == 5] = 4
        quadrant_mask[cell_mask == 0] = 0
        mtoc_quad = quadrant_mask[mtoc_position[1], mtoc_position[0]]
        for quad_num in range(1, 5):
            intensity_by_quad[i, quad_num - 1, 0] = np.sum(IF[(peripheral_binary_mask==1) & (quadrant_mask == quad_num)])
            if quad_num == mtoc_quad:
                intensity_by_quad[i, mtoc_quad - 1, 1] += 1
            intensity_by_quad[i, quad_num - 1, 0] /= np.sum(height_map[(peripheral_binary_mask ==1) & (quadrant_mask == quad_num)]) * constants.VOLUME_COEFFICIENT
    return intensity_by_quad


def get_mtoc_quad(file_handler, image):
    """get mtoc quadrant from image descriptors of johnatan"""
    attributes = file_handler[image].attrs.keys()
    if 'mtoc_quad' not in attributes:
        return -1
    mtoc_quad = file_handler[image].attrs['mtoc_quad']
    return mtoc_quad

# rewrite so that it works on arrays
def get_quadrant(file_handler, image, position):
    quadrants = get_quadrants(file_handler, image)
    assert quadrants.size, 'mtoc_quadrant not defined for %r' % image
    return quadrants[position[0], position[1]]

def is_in_cytoplasm(file_handler, image, position):
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    return (cell_mask[position[0], position[1]] == 1) and (nucleus_mask[position[0], position[1]] == 0)

def is_in_z_lines(file_handler, image, position):
    z_lines_mask = get_z_lines_mask(file_handler, image)
    return (z_lines_mask[position[0], position[1]] == 1)

def get_height_map(file_handler, image):
    descriptor = image + '/height_map'
    if descriptor not in file_handler:
        return np.empty(shape=(0, 0))
    height_map = np.array(file_handler[descriptor])
    return height_map

def get_cell_volume(file_handler, image):
    descriptor = image + '/cell_volume'
    attributes = file_handler[image].attrs.keys()
    if 'cell_volume' not in attributes:
        return -1
    cell_volume = file_handler[descriptor]
    return cell_volume

def set_cell_volume(file_handler, image):
    cell_volume = get_cell_volume(file_handler, image)
    assert (cell_volume == -1), 'cell_volume already defined for %r' % image
    cell_mask = get_cell_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    height_map = height_map + constants.VOLUME_OFFSET
    file_handler[image].attrs['cell_volume'] = height_map[np.where(cell_mask[:] == 1)].sum()

def get_nucleus_volume(file_handler, image):
    descriptor = image + '/nucleus_volume'
    attributes = file_handler[image].attrs.keys()
    if 'nucleus_volume' not in attributes:
        return -1
    nucleus_volume = np.array(file_handler[descriptor])
    return nucleus_volume

# need height map to compute, still in progress
# returns the nucleus volume in pixels
def set_nucleus_volume(file_handler, image):
    nucleus_volume = get_nucleus_volume(file_handler, image)
    assert (nucleus_volume != -1), 'nucleus_volume already defined for %r' % image
    nucleus_mask = get_nucleus_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    height_map = height_map + constants.VOLUME_OFFSET
    file_handler[image].attrs['nucleus_volume'] = height_map[np.where(nucleus_mask[:] == 1)].sum()


def compute_protein_cytoplasmic_spread(file_handler, image,path_data):
    cell_mask = get_cell_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    IF = get_IF(file_handler,image)
    height_map += 1
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - nucleus_centroid[0]
    ds2 = np.matlib.repmat(np.asmatrix(np.arange(0, 512).reshape(512, 1)), 1, 512) - nucleus_centroid[1]
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)
    height_map = np.multiply(height_map,cell_mask)
    height_map_dist = np.multiply(height_map,dsAll)
    S = height_map_dist.sum() / height_map.sum()
    dist_IF = np.multiply(IF,dsAll)
    val = dist_IF.sum() / (IF.sum() * S)
    return val

def compute_protein_cytoplasmic_spread_2D(file_handler, image,path_data):
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    IF = get_IF(file_handler,image)
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - nucleus_centroid[0]
    ds2 = np.matlib.repmat(np.asmatrix(np.arange(0, 512).reshape(512, 1)), 1, 512) - nucleus_centroid[1]
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)
    plane_map_dist = np.multiply(cell_mask,dsAll)
    S = plane_map_dist.sum() / cell_mask.sum()
    dist_IF = np.multiply(IF,dsAll)
    val = dist_IF.sum() / (IF.sum() * S)
    return val

def compute_mrna_cytoplasmic_total(file_handler, image):
    spots = get_spots(file_handler, image)
    counter=0
    for i in range(len(spots)):
        if is_in_cytoplasm(file_handler, image, [spots[i, 1], spots[i, 0]]):
            counter+=1
    return counter

def compute_protein_cytoplasmic_total(file_handler, image,path_data):
    print(image)
    if "tp" in image:
        molecule = image.split("/")[1]
        gene = image.split("/")[2].split("_")[0]
        type = image.split("/")[2].split("_")[1]
        image_n = image.split("/")[4]
        IF = get_IF_image_z_summed(molecule,gene,type,image_n,path_data)
    else:
        IF = get_IF(file_handler,image)
    cell_mask=get_cell_mask(file_handler, image)
    nucleus_mask= get_nucleus_mask(file_handler,image)
    return IF[(cell_mask==1) & (nucleus_mask==0)].sum()

# Compare cytoplasmic spread cell with 3D cytoplasmic mrna spread
# to evaluate degree of spread
def compute_mrna_cytoplasmic_spread(file_handler, image):

    cell_mask = get_cell_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    height_map += 1
    spots = get_spots(file_handler, image)

    # Compute all possible distance in a matrix [512x512]
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - nucleus_centroid[0]
    ds2 = np.matlib.repmat(np.asmatrix(np.arange(0, 512).reshape(512,1)), 1, 512) - nucleus_centroid[1]
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)

    # Computing spots distance from nucleus centroid
    points_dist_list = []
    counter = 0
    for i in range(len(spots)):
        if is_in_cytoplasm(file_handler,image,[spots[i,1],spots[i,0]]):
            dist=0.0
            for j in range(2):
                if j == 0:
                    dist+=(spots[i, j] - nucleus_centroid[0])**2
                elif j == 1:
                    dist += (spots[i, j] - nucleus_centroid[1]) ** 2
            points_dist_list.append(math.sqrt(dist))
            counter+=1
    points_dists = np.matrix(points_dist_list)
    points_dists = points_dists.reshape((counter, 1))
    height_map = np.multiply(height_map, cell_mask)
    height_map_dist=np.multiply(height_map,dsAll)

    # S : Average distance of a cytoplasmic voxel from the nucleus centroid
    S = height_map_dist.sum() / height_map.sum()

    # val is average 2D distance from the nucleus centroid of cytoplasmic mRNAs
    # normalized by the cytoplasmic cell spread (taking a value 1 when mRNAs are evenly
    # distributed across the cytoplasm).
    normalized_average_2d_distance = np.mean(points_dists) / S
    return normalized_average_2d_distance

# Compare cytoplasmic spread cell with 2D cytoplasmic mrna spread
# to evaluate degree of spread
def compute_mrna_cytoplasmic_spread_2D(file_handler, image):
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    spots = get_spots(file_handler, image)
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - nucleus_centroid[0]
    ds2 = np.matlib.repmat(np.asmatrix(np.arange(0, 512).reshape(512,1)), 1, 512) - nucleus_centroid[1]
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)

    # Computing spots distance from nucleus centroid
    points_dist_list = []
    counter = 0
    for i in range(len(spots)):
        if is_in_cytoplasm(file_handler,image,[spots[i,1],spots[i,0]]):
            dist=0.0
            for j in range(2):
                if j == 0:
                    dist+=(spots[i, j] - nucleus_centroid[0])**2
                elif j == 1:
                    dist += (spots[i, j] - nucleus_centroid[1]) ** 2
            points_dist_list.append(math.sqrt(dist))
            counter+=1
    points_dists = np.matrix(points_dist_list)
    points_dists = points_dists.reshape((counter, 1))
    plane_map_dist=np.multiply(cell_mask,dsAll)

    # S : Average distance of a cytoplasmic voxel from the nucleus centroid
    S = plane_map_dist.sum() / cell_mask.sum()

    # val is average 2D distance from the nucleus centroid of cytoplasmic mRNAs
    # normalized by the cytoplasmic cell spread (taking a value 1 when mRNAs are evenly
    # distributed across the cytoplasm).
    normalized_average_distance = np.mean(points_dists) / S
    return normalized_average_distance


def compute_isoline_area(file_handler, output_file_handler, image):
    isoline_area = []
    periph_dist_map = get_cell_mask_distance_map(output_file_handler, image)
    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    periph_dist_map[(periph_dist_map == 0) & (cell_mask == 1) & (nucleus_mask == 0)] = 1
    nucleus_area = get_nucleus_area(output_file_handler, image)
    for i in range(101):
        tmp_mask = np.array(periph_dist_map, copy=True)
        tmp_mask[tmp_mask <= i] = 0
        tmp_mask[(tmp_mask > i) & (tmp_mask <= 100)] = 1
        isoline_area.append(tmp_mask.sum() * math.pow((1 / constants.SIZE_COEFFICIENT), 2) + nucleus_area)
    return isoline_area


def compute_cell_volume(file_handler,image):
    height_map=get_height_map(file_handler,image)
    return np.sum(height_map) * constants.VOLUME_COEFFICIENT


'''


def set_mrna_cytoplasmic_spread(file_handler, image):
    cell_mask = get_cell_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    nucleus_centroid = get_nucleus_centroid(file_handler, image)
    height_map += 1
    spots = get_spots(file_handler, image)
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - nucleus_centroid[0]
    ds2=ds1.transpose()
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)
    cytoplasm_spots=[]
    for spot in spots:
        if is_in_cytoplasm(file_handler,image,[spot[1],spot[0]]):
            cytoplasm_spots.append(spot)
    spots=np.matrix(cytoplasm_spots)
    spots=spots.reshape((len(cytoplasm_spots),3))
    height_map = height_map * cell_mask
    S = (height_map * dsAll).sum() / height_map.sum()
    points_dists = np.power(spots[:, 0:2] - 256, 2)
    points_dists = np.sqrt(points_dists[1, :].sum())
    val = np.mean(points_dists) / S
    file_handler[image].attrs['cytoplasmic_spread'] = val
    
    
def set_cell_mask_dist_map_meshgrid(file_handler, image):
    import time
    start = time.time()
    print("hello")

    cell_mask = get_cell_mask(file_handler, image)
    nucleus_mask = get_nucleus_mask(file_handler, image)
    cytoplasm_mask = (cell_mask == 1) & (nucleus_mask == 0)
    nucleus_centroid = get_nucleus_centroid(file_handler, image).transpose()

    angle_slopex = (np.sin((np.arange(360)*2*math.pi/360)))
    angle_slopey = (np.cos((np.arange(360)*2*math.pi/360)))
    cell_radiusx = np.arange(400)
    cell_radiusy = np.arange(400)

    line_segmentx = np.kron(angle_slopex, cell_radiusx)
    line_segmenty = np.kron(angle_slopey, cell_radiusy)
    line_segmentx += int(round(nucleus_centroid[0]))
    line_segmenty += int(round(nucleus_centroid[1]))
    test = np.zeros((144000, 2))
    test[:, 0] = line_segmentx
    test[:, 1] = line_segmenty
    test_reshape = test.reshape((360, 400, 2))

    end = time.time()
    start = time.time()
    cyt_segments = np.zeros((144000, 1))
    nuc_segments = np.zeros((144000, 1))

    for i in range(144000):
        x = test[i,0]
        y = test[i,1]
        if not (x < 0 or x > 511 or y < 0 or y > 511):
            cyt_segments[i, :] = cytoplasm_mask[y,x]
            nuc_segments[i, :] = nucleus_mask[y,x]
        else:
            cyt_segments[i, :] = 0
            nuc_segments[i, :] = 0

    plt.imshow(cell_mask, cmap='gray')
    contours = measure.find_contours(nucleus_mask, 0.8)
    for n, contour in enumerate(contours):
        plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    plt.plot(test_reshape[90,:,0], test_reshape[90,:,1], color='red', linewidth=2)
    plt.show()
    end = time.time()
    print(end - start)

    start = time.time()
    print("hello")
    tis = np.zeros((400, 360, 2))
    ti5 = f5(*np.meshgrid(angle_slopey, cell_radiusy, sparse=False))
    ti5x = ti5 + int(round(nucleus_centroid[0]))
    ti5y = ti5 + int(round(nucleus_centroid[1]))
    tis[:, :, 0]=ti5x
    tis[:, :, 1] = ti5y

    print(tuple(tis[0, 0, :]))
    plt.imshow(cytoplasm_mask, cmap='gray')
    contours = measure.find_contours(nucleus_mask, 0.8)
    for n, contour in enumerate(contours):
        plt.plot(contour[:, 1], contour[:, 0], color='red', linewidth=2)
    plt.show()

    print(cytoplasm_mask[tis[0, 0, 1],tis[0, 0, 0]])
    print(nucleus_mask[tis[0, 0, 1], tis[0, 0, 0]])
    end = time.time()
    print(end - start)
    start = time.time()
    print("hello")

    test_contour = np.zeros((36000, 2))
    cell_contourx = np.arange(100)
    cell_contoury = np.arange(100)
    contours_pointx = np.kron(angle_slopex, cell_contourx)
    contours_pointy = np.kron(angle_slopey, cell_contoury)
    contours_pointx += int(round(nucleus_centroid[0]))
    contours_pointy += int(round(nucleus_centroid[1]))
    test_contour[:, 0] = contours_pointx
    test_contour[:, 1] = contours_pointy
    test_contour = test_contour.reshape((360, 100, 2))

    end = time.time()
    print(end - start)
    sys.exit()



def set_protein_cytoplasmic_spread(file_handler, image):
    cell_mask = get_cell_mask(file_handler, image)
    height_map = get_height_map(file_handler, image)
    IF = get_IF(file_handler, image)
    height_map += 1
    ds1 = np.matlib.repmat(range(0, 512), 512, 1) - 256
    ds2 = np.matlib.repmat(np.array(range(0, 512)).transpose(), 1, 512) - 256
    dsAll = np.power(ds1, 2) + np.power(ds2, 2)
    dsAll = np.sqrt(dsAll)
    height_map = height_map * cell_mask
    S = (height_map * dsAll).sum() / height_map.sum()
    I = IF[2].sum()
    tmp = I * dsAll
    val = tmp.sum() / I.sum() * S
    print(val)


'''


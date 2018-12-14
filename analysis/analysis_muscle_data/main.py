#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import interpolate
import src.constants as cst
import src.acquisition_descriptors as adsc
import src.image_descriptors as idsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])


def compute_minimal_distance(segment_summed):
    for i in range(15):
        if segment_summed[i]!=0:
            return i


if __name__ == "__main__":

    # Required descriptors: cell_area (built from cell_mask), spots
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8

    basic_file_path = path.analysis_data_dir + 'basics_muscle_data.h5'
    all_median_profiles=[]
    molecule_type = ['/mrna']
    genes=['actn2','gapdh']
    timepoints=['mature']
    colors = ['#0A3950', '#1E95BB', '#A1BA6D']

    # with h5py.File(basic_file_path, "a") as file_handler:
    #     dist = []
    #     for gene in genes:
    #         image_list=[]
    #         h_star_l = []
    #         for timepoint in timepoints:
    #             image_tab=[]
    #             image_tabl=np.zeros((len(file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]),13))
    #             image_count=0
    #             total_profile=[]
    #             dist = []
    #             for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
    #                 image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
    #                 print(image)
    #                 cell_mask=idsc.get_cell_mask(file_handler, image)
    #                 z_lines_mask=idsc.get_z_lines_mask(file_handler, image)
    #                 nucleus_mask=idsc.get_nucleus_mask(file_handler,image)
    #                 spots=idsc.get_spots(file_handler,image)
    #                 cell_area = cell_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
    #                 z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
    #                 cytoplasm_area=cell_area-z_lines_area
    #                 #for each spot compute distance to nearest zlines
    #                 straight_line_x = np.arange(z_lines_mask.shape[1])
    #                 straight_line_y = []
    #                 h_line = np.zeros(z_lines_mask.shape[1])
    #                 for x in range(z_lines_mask.shape[1]):
    #                     straight_line_y.append(z_lines_mask.shape[0]/2)
    #                 dist_count = 0
    #                 z_line_area = False
    #                 for x in range(z_lines_mask.shape[1]):
    #                     if z_lines_mask[straight_line_y[x], straight_line_x[x]] == 0:
    #                         z_line_area = False
    #                         if z_line_area == False:
    #                             dist_count += 1
    #                             z_line_area = True
    #                     else:
    #                         if dist_count < 20 and dist_count > 10:
    #                             dist.append(dist_count)
    #                         z_line_area = True
    #                         dist_count = 0
    #                     #h_line[x] = z_lines_mask[straight_line_y[x], straight_line_x[x]]
    #                 #print(np.median(dist))
    #
    #
    #                 image_count += 1
    #
    #     median_distance_between_z_line=int(np.median(np.array(dist)))

    # #compute degree of clustering
    # with h5py.File(basic_file_path, "a") as file_handler:
    #     base = math.log(0.5)
    #     mrna_median = []
    #     mrna_err = []
    #     for gene in genes:
    #         image_list=[]
    #         h_star_l = []
    #         for timepoint in timepoints:
    #             for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
    #                 image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
    #                 print(image)
    #                 z_lines_mask=idsc.get_z_lines_mask(file_handler, image)
    #                 spots=idsc.get_spots(file_handler,image)
    #                 h_star=helps.clustering_index_point_process_2d(spots,z_lines_mask)
    #                 print(h_star)
    #                 descriptor = image + '/h_star'
    #                 # output_file_handler[descriptor][:] = h_star
    #                 file_handler.create_dataset(descriptor, data=h_star, dtype=np.float32)
    #                 d = np.array(h_star[h_star > 1] - 1).sum()
    #                 h_star_l.append(d)
    #         mrna_median.append(math.log(np.median(h_star_l)) - base)
    #         # mrna_median.append(numpy.median(dof))
    #         err = np.median(np.abs(np.tile(np.median(h_star_l), (1, len(h_star_l))) - h_star_l))
    #         error_median = math.log(np.median(h_star_l) + err)
    #         error_median = error_median - math.log(np.median(h_star_l)) - base
    #         mrna_err.append(error_median)
    #     fig = plt.figures()
    #     ax = plt.axes()
    #     ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    #     N = len(genes)
    #     ind = np.arange(N)
    #     width = 0.35
    #     # colors = ['blue', 'lightblue', 'lightgreen', 'orange', 'red', 'yellow']
    #     colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']
    #
    #     ## the bars
    #     rects1 = ax.bar(ind, mrna_median, width,
    #                     color=colors,
    #                     yerr=mrna_err,
    #                     error_kw=dict(elinewidth=1, ecolor='black'))
    #     # axes and labels
    #     ax.set_xlim(-width, len(ind) + width)
    #     ax.set_ylim(0, 20)
    #     ax.set_ylabel('Degree of clustering(log)')
    #     ax.set_title('Mrna degree of clustering')
    #     xTickMarks = ["" for i in range(0, 2)]
    #     ax.set_xticks(ind)
    #     xtickNames = ax.set_xticklabels(xTickMarks)
    #     plt.legend([gene for gene in genes], loc='upper right')
    #     ax.legend(rects1, genes)
    #     figname = path.analysis_dir + 'analysis_muscle_data/figures/mrna_degree_of_clustering.svg'
    #     fig.savefig(figname)
    #     plt.show()
    #     print(h_star_l)

    # with h5py.File(basic_file_path, "a") as file_handler:
    #     for gene in genes:
    #         image_list=[]
    #         h_star_l = []
    #         for timepoint in timepoints:
    #             image_tab=[]
    #             image_tabl=np.zeros((len(file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]),13))
    #             image_count=0
    #             total_profile=[]
    #             dist = []
    #             for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
    #                 image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
    #                 print(image)
    #                 cell_mask=idsc.get_cell_mask(file_handler, image)
    #                 z_lines_mask=idsc.get_z_lines_mask(file_handler, image)
    #                 # cytoplasm_mask=z_lines_mask[z_lines_mask==1]=2
    #                 # cytoplasm_mask=z_lines_mask[z_lines_mask == 0] = 1
    #                 # cytoplasm_mask=z_lines_mask[z_lines_mask == 2] = 0
    #                 # cytoplasm_mask=z_lines_mask[cell_mask==0]=0
    #                 nucleus_mask=idsc.get_nucleus_mask(file_handler,image)
    #                 spots=idsc.get_spots(file_handler,image)
    #                 cell_area = cell_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
    #                 z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
    #                 cytoplasm_area=cell_area-z_lines_area
    #                 # h_star=helps.clustering_index_point_process_2d(spots,z_lines_mask)
    #                 # print(h_star)
    #                 # d = np.array(h_star[h_star > 1] - 1).sum()
    #                 # h_star_l.append(d)
    #
    #                 #for each spot compute distance to nearest zlines
    #
    #                 spots_count=0
    #                 image_spots_min_distance=np.zeros(len(spots))
    #                 for spot in spots:
    #                     if z_lines_mask[spot[1], spot[0]] == 0:
    #                         total_segment = np.zeros((360, median_distance_between_z_line))
    #                         for degree in range(360):
    #                             z_line_segment = np.zeros(median_distance_between_z_line)
    #
    #                             line=np.zeros((median_distance_between_z_line,2))
    #                             angle = degree * 2 * math.pi / 360
    #                             x_slope,y_slope = math.sin(angle), math.cos(angle)
    #                             for point in range(median_distance_between_z_line):
    #                                 x = int(round(spot[0] + point * x_slope))
    #                                 # print(x)
    #                                 y = int(round(spot[1] + point * y_slope))
    #
    #                                 line[point,0]=x
    #                                 line[point,1]=y
    #                                 #print(x,y)
    #                                 #print(z_lines_mask.shape)
    #
    #                                 if (x >= 0 and x < z_lines_mask.shape[1] and y >= 0 and y < z_lines_mask.shape[0]):
    #                                     z_line_segment[point] = z_lines_mask[y, x]
    #                                     total_segment[degree,point]=z_lines_mask[y, x]
    #                                 else:
    #                                     z_line_segment[point] = 0
    #                                     total_segment[degree,point]=0
    #
    #                             #print(z_line_segment)
    #                             #plt.figures(figsize=(20, 20))
    #                             #plt.imshow(z_lines_mask)
    #                             # straight_line_x=np.arange(1100)
    #                             # straight_line_y = []
    #                             # h_line=np.zeros(1100)
    #                             # for x in range(1100):
    #                             #     straight_line_y.append(100)
    #                             # for x in range(1100):
    #                             #     h_line[x]=z_lines_mask[straight_line_y[x], straight_line_x[x]]
    #                             #print(h_line)
    #                             #plt.plot(straight_line_x,straight_line_y,color='yellow')
    #                             #plt.show()
    #
    #
    #                         distance=compute_minimal_distance(np.sum(total_segment,axis=0))
    #                         image_spots_min_distance[spots_count]=distance
    #                     else:
    #                         image_spots_min_distance[spots_count] =0
    #                     spots_count += 1
    #
    #                 #print(image_spots_min_distance)
    #                 #plt.hist(image_spots_min_distance,range=(1,15),bins=15)
    #                 #plt.show()
    #                 z_line_distance_profile = np.zeros(median_distance_between_z_line)
    #                 for i in range(median_distance_between_z_line):
    #                     #print(len(np.where(image_spots_min_distance == i)[0])/float(len(image_spots_min_distance)))
    #                     z_line_distance_profile[i] =float(len(np.where(image_spots_min_distance == i )[0]))/float(len(image_spots_min_distance))
    #
    #                 #print(z_line_distance_profile)
    #                 total_profile.append(z_line_distance_profile)
    #                 #print(z_line_distance_profile)
    #                 image_count += 1
    #             total_profile=np.array(total_profile).reshape((image_count,median_distance_between_z_line))
    #             #print(total_profile)
    #             print(np.median(total_profile,axis=0))
    #             all_median_profiles.append(np.median(total_profile,axis=0))
    #
    #         print(h_star_l)
    # df=pd.DataFrame(all_median_profiles)
    # df.to_csv("all_median_profiles.csv")




    molecule_type = ['/mrna']
    genes = ['actn2']
    timepoints = ['immature']
    with h5py.File(basic_file_path, "a") as file_handler:
        for gene in genes:
            image_list = []
            for timepoint in timepoints:
                image_tab = []
                image_tabl = np.zeros(
                    (len(file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]), 13))
                image_count = 0
                total_profile = []
                for im in file_handler[molecule_type[0] + '/' + gene + '/' + timepoint]:
                    image = molecule_type[0] + '/' + gene + '/' + timepoint + '/' + im
                    print(image)
                    cell_mask = idsc.get_cell_mask(file_handler, image)
                    z_lines_mask = idsc.get_z_lines_mask(file_handler, image)
                    # cytoplasm_mask=z_lines_mask[z_lines_mask==1]=2
                    # cytoplasm_mask=z_lines_mask[z_lines_mask == 0] = 1
                    # cytoplasm_mask=z_lines_mask[z_lines_mask == 2] = 0
                    # cytoplasm_mask=z_lines_mask[cell_mask==0]=0
                    nucleus_mask = idsc.get_nucleus_mask(file_handler, image)
                    spots = idsc.get_spots(file_handler, image)
                    cell_area = cell_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                    z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                    cytoplasm_area = cell_area - z_lines_area

                    # for each spot compute distance to nearest zlines

                    spots_count = 0
                    image_spots_min_distance = np.zeros(len(spots))
                    for spot in spots:
                        if z_lines_mask[spot[1], spot[0]] == 0:
                            total_segment = np.zeros((360, 15))
                            for degree in range(360):
                                z_line_segment = np.zeros(15)

                                line = np.zeros((15, 2))
                                angle = degree * 2 * math.pi / 360
                                x_slope, y_slope = math.sin(angle), math.cos(angle)
                                for point in range(15):
                                    x = int(round(spot[0] + point * x_slope))
                                    # print(x)
                                    y = int(round(spot[1] + point * y_slope))

                                    line[point, 0] = x
                                    line[point, 1] = y
                                    # print(x,y)
                                    # print(z_lines_mask.shape)

                                    if x >= 0 and x < z_lines_mask.shape[1] and y >= 0 and y < z_lines_mask.shape[0]:
                                        z_line_segment[point] = z_lines_mask[y, x]
                                        total_segment[degree, point] = z_lines_mask[y, x]
                                    else:
                                        z_line_segment[point] = 0
                                        total_segment[degree, point] = 0

                                        # print(z_line_segment)
                                        # plt.figures(figsize=(20, 20))
                                        # plt.imshow(z_lines_mask)
                                        # plt.plot(line[:,0], line[:,1], color='yellow')
                                        #
                                        # plt.show()

                            dist = compute_minimal_distance(np.sum(total_segment, axis=0))
                            image_spots_min_distance[spots_count] = dist
                        else:
                            image_spots_min_distance[spots_count] = 0
                        spots_count += 1

                    # print(image_spots_min_distance)
                    # plt.hist(image_spots_min_distance,range=(1,15),bins=15)
                    # plt.show()
                    z_line_distance_profile = np.zeros(15)
                    for i in range(15):
                        #print(len(np.where(image_spots_min_distance == i)[0]) / float(len(image_spots_min_distance)))
                        z_line_distance_profile[i] = float(
                            len(np.where(image_spots_min_distance == i)[0])) / float(
                            len(image_spots_min_distance))

                    # print(z_line_distance_profile)
                    total_profile.append(z_line_distance_profile)
                    # print(z_line_distance_profile)
                    image_count += 1
                # print(total_profile)
                total_profile = np.array(total_profile).reshape((image_count, 15))
                # print(total_profile)
                print(np.median(total_profile, axis=0))
                all_median_profiles.append(np.median(total_profile,axis=0))

    df=pd.DataFrame(all_median_profiles)
    df.to_csv("all_median_profiles_immature.csv")


    all_median_profiles=[]
    plot_name = path.analysis_dir + 'analysis_muscle_data/figures/z_line_distance' + str(15) + 'contours_immature.png'
    figure_title = 'z line spots distance profile'
    genes=["actn2 mature"]

    df=pd.read_csv("all_median_profiles_immature.csv",index_col=0)

    all_median_profiles.append(df.ix[0].values)
    #all_median_profiles.append(df.ix[1].values)



    # all_cumul_median_profiles=[]
    # print(np.cumsum(all_median_profiles[0]))
    # for i in range(15):
    #     new_distance=all_median_profiles[0][i]
    #     if i==0:
    #         old_dist=0
    #     else:
    #         old_dist=all_median_profiles[0][i-1]
    #     print(dist+all_median_profiles[0][i-1])
    #
    #     all_cumul_median_profiles.append(all_median_profiles[0][i]+dist)

    # print(all_cumul_median_profiles)
    stan.plot_profile(all_median_profiles, genes, 15, plot_name, figure_title, colors, True)
    #stan.plot_profile([np.cumsum(all_median_profiles[0]),np.cumsum(all_median_profiles[1])], genes, 15, plot_name, figure_title, colors, True)

























































                # # print(cell_area)
                    # # print(z_lines_area)
                    # # print(cytoplasm_area)
                    # # z_lines_mask[z_lines_mask==1]=2
                    # # z_lines_mask[z_lines_mask == 0] = 1
                    # # z_lines_mask[z_lines_mask == 2] = 0
                    # # z_lines_mask[cell_mask==0]=0
                    # #
                    # # xs = spots[:, 0]
                    # # ys = spots[:, 1]
                    # # plt.imshow(z_lines_mask)
                    # # # # plt.imshow(cell_mask)
                    # # plt.scatter(xs, ys, color='yellow', marker=".", facecolors='none', linewidths=0.5)
                    # # plt.show()
                    # #
                    # # from scipy.ndimage import morphology
                    # # for i in range(10):
                    # #     z_lines_mask= morphology.binary_erosion(z_lines_mask).astype(z_lines_mask.dtype)
                    # #
                    # #     xs = spots[:, 0]
                    # #     ys = spots[:, 1]
                    # #     plt.imshow(z_lines_mask)
                    # #
                    # #     plt.scatter(xs, ys, color='yellow', marker=".", facecolors='none', linewidths=0.5)
                    # #
                    # #     plt.show()
                    # #     count=0
                    # #     for spot in spots:
                    # #         if z_lines_mask[spot[1], spot[0]] == 1:
                    # #             count+=1
                    # #     print(count)
                    # #print(len(spots))
                    # from skimage.morphology import square
                    # from skimage.morphology import binary_dilation
                    #
                    # tab=[]
                    # z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                    # count = 0
                    # for spot in spots:
                    #     if z_lines_mask[spot[1], spot[0]] == 1:
                    #         count += 1
                    # print((float(count) / z_lines_area))
                    # image_tabl[image_count, 0] = (float(count) / z_lines_area)
                    # for i in range(1,13):
                    #     z_lines_mask = binary_dilation(z_lines_mask,square(2))
                    #     xs = spots[:, 0]
                    #     ys = spots[:, 1]
                    #     plt.figures(figsize=(20, 20))
                    #     plt.imshow(z_lines_mask)
                    #
                    #     plt.scatter(xs, ys, color='yellow', marker=".", facecolors='none', linewidths=0.5)
                    #
                    #     plt.show()
                    #     count = 0
                    #     for spot in spots:
                    #         if z_lines_mask[spot[1], spot[0]] == 1:
                    #             count+=1
                    #     #print(count)
                    #     z_lines_area = z_lines_mask.sum() * math.pow((1 / cst.SIZE_COEFFICIENT), 2)
                    #     tab.append(count/z_lines_area)
                    #     #print((float(count)/z_lines_area)/len(spots))
                    #     image_tabl[image_count,i]=(float(count)/z_lines_area)#/len(spots)
                    #
                    # #plt.plot(tab)
                    # #plt.show()
                    # image_tab.append(tab)
                    # image_count+=1
                    #
                    #
                    #
                    #
                    # # xs = spots[:, 0]
                    # # ys = spots[:, 1]
                    # # plt.imshow(z_lines_mask)
                    # #
                    # # plt.scatter(xs, ys, color='yellow', marker=".", facecolors='none', linewidths=0.5)
                    # #
                    # # spots_in_z_line=0
                    # # for spot in spots:
                    # #     if idsc.is_in_z_lines(file_handler, image, [spot[1], spot[0]]):
                    # #         spots_in_z_line+=1
                    # # spots_n=len(spots)
                    # # print(spots_in_z_line)
                    # # print(len(spots))
                    #
                    # # spots_density_z_lines=spots_in_z_line/z_lines_area
                    # # spots_density_cell=cell_area/spots_n
                    # # spots_density_cytoplasm = cytoplasm_area / (spots_n - spots_in_z_line)
                    # #
                    # # print(spots_density_z_lines)
                    # # print(spots_density_cell)
                    # # print(spots_density_cytoplasm)
                    # # plt.show()
                    # image_list.append(image)
                # sys.exit()
                # print(image_tabl)
                #
                # x = np.arange(0, 13, 0.01)
                # print(len(np.mean(image_tabl,axis=0)))
                # spl = interpolate.UnivariateSpline(np.arange(13), np.mean(image_tabl,axis=0), k=3)
                # new = spl(x)
                # plt.plot(x,new,linestyle="-", linewidth=2)
                # plt.show()
                # plt.boxplot(image_tabl)
                # plt.show()



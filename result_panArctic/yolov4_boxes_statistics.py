#!/usr/bin/env python
# Filename: yolov4_boxes_statistics.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 17 April, 2023
"""


import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import time
import vector_gpd
from shapely.geometry import Point


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter




global_bin_size = 50   # remember to change this one
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.

    # print(global_bin_size)
    # s = str(100 * y*bin_size)
    s = "%.0f"%(100 * y*global_bin_size)

    # The percent symbol needs escaping in latex
    # if matplotlib.rcParams['text.usetex'] is True:
    #     return s + r'$\%$'
    # else:
    #     return s + '%'
    if matplotlib.rcParams['text.usetex'] is True:
        return s
    else:
        return s

def histogram2logfile(value_list,bins,hist_tag=None):
    if hist_tag is not None:
        basic.outputlogMessage('the following is the histogram information of %s'%hist_tag)
    # output hist, min, max, average, accumulate percentage
    np_hist,bin_edges = np.histogram(value_list, bins=bins)
    basic.outputlogMessage("total count: " + str(len(value_list)))
    basic.outputlogMessage("np_hist: " + str(np_hist))
    basic.outputlogMessage("min value: " + str(min(value_list)))
    basic.outputlogMessage("max value: " + str(max(value_list)))
    basic.outputlogMessage("average value: " + str(sum(value_list)/float(len(value_list))))
    basic.outputlogMessage("total value: " + str(sum(value_list)))
    if len(value_list) != np.sum(np_hist):
        basic.outputlogMessage('warning: the count (%d) of input is not equal to the count (%d)'
                               ' in histogram'%(len(value_list),int(np.sum(np_hist))))
    acc_per = np.cumsum(np_hist)/np.sum(np_hist)
    basic.outputlogMessage("accumulate percentage: " + str(acc_per))

def plot_hist(value_list,save_path,bins = None, b_xtick_rotate=False, x_tick=None,x_label=None):
    # plot a histogram
    # bin_count = 12
    if bins is None:
        bins = np.arange(0.5, 12.5+1, 1)
    print('bins:', bins)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    n, bins, patches = ax.hist(value_list, bins=bins, alpha=0.75, ec="black", linewidth='1.5',
                               color='grey', hatch='', rwidth=1,density = True)  # density = True, # label=labels,

    # ax.legend(prop={'size': 12})
    if x_tick is not None:
        # ticks = None, labels = None
        plt.xticks(ticks=x_tick,labels=x_tick)
    ax.tick_params(axis='both', which='both', direction='out', length=7, labelsize=20)  # ,width=50 #,
    # if xlabelrotation is not None:
    #     ax.tick_params(axis='x', labelrotation=90)
    if b_xtick_rotate:
        ax.tick_params(axis='x', labelrotation=90)

    # if ylim is not None:
    #     ax.set_ylim(ylim)

    formatter = FuncFormatter(to_percent)
    # Set the formatter
    plt.gca().yaxis.set_major_formatter(formatter)

    if x_label is not None:
        plt.xlabel(x_label,fontsize=20)
    plt.ylabel('Percentage (%)',fontsize=20)

    plt.gcf().subplots_adjust(bottom=0.15)
    # plt.grid(True)
    # plt.show()
    plt.savefig(save_path)  #
    print('save fig to ', os.path.abspath(save_path))

def calculate_width_height_of_boxes(boxes, box_ids):
    # assume height >= width
    height_list = []
    width_list = []

    count_not_rectangle = 0
    # print(boxes[0])
    for id, box in zip(box_ids, boxes):
        points = list(box.boundary.coords)
        if len(points) != 5:
            print('id:',id, ',points of non-rectangle:',len(points), ',will use its minimum_rotated_rectangle')
            count_not_rectangle += 1
            min_bounding_box = box.minimum_rotated_rectangle
            points = list(min_bounding_box.boundary.coords)

        point1 = Point(points[0])
        point2 = Point(points[1])
        point3 = Point(points[2])

        width = point1.distance(point2)
        height = point2.distance(point3)
        # assume height >= width
        if width > height:
            width, height = height, width
        height_list.append(height)
        width_list.append(width)

        # break   # for test
    print('total polygons: %d: rectangle: %d, non-rectangle: %d'%(len(boxes), len(boxes)-count_not_rectangle ,count_not_rectangle))

    return width_list, height_list


def main():
    boxes_shp = os.path.expanduser('~/Data/labelearth.colorado.edu/data/thawslump_boxes/pan_arctic_thawslump_after_webValidation_thr0.5.shp')

    # polygons = vector_gpd.read_polygons_gpd(boxes_shp,b_fix_invalid_polygon=False)
    polygons = vector_gpd.read_shape_gpd_to_NewPrj(boxes_shp,'EPSG:3413')
    poly_ids = vector_gpd.read_attribute_values_list(boxes_shp,'id')
    print('read %d polygons'%len(polygons))

    box_width_list, box_height_list = calculate_width_height_of_boxes(polygons, poly_ids)

    # calculate length and width
    # print(poly_ids[0],'height:',height, 'width',width)

    width_bins = np.arange(min(box_width_list), max(box_width_list)+global_bin_size, global_bin_size)
    # x_tick=year_unique,
    plot_hist(box_width_list, 'YOLOv4_boxes_width_hist.jpg', bins=width_bins, b_xtick_rotate=True, x_label='Width (m)')
    histogram2logfile(box_width_list, width_bins, hist_tag='Width')
    io_function.move_file_to_dst('processLog.txt', 'YOLOv4_boxes_width_hist_log.txt', overwrite=True)


    height_bins = np.arange(min(box_height_list), max(box_height_list)+global_bin_size, global_bin_size)
    # x_tick=year_unique,
    plot_hist(box_height_list, 'YOLOv4_boxes_height_hist.jpg', bins=height_bins, b_xtick_rotate=True, x_label='Height (m)')
    histogram2logfile(box_height_list, height_bins, hist_tag='Height')
    io_function.move_file_to_dst('processLog.txt', 'YOLOv4_boxes_height_hist_log.txt', overwrite=True)

    pass

if __name__ == '__main__':
    main()
    pass
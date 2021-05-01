#!/usr/bin/env python
# Filename: dem_headwall_extraction 
"""
introduction: extraction RTS headwall from DEM

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 30 April, 2021
"""

import os,sys
from optparse import OptionParser
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import vector_gpd

import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import raster_io

import cv2
import numpy as np
import geopandas as gpd
import pandas as pd

from dem_common import dem_headwall_shp_dir


def extract_headwall_from_slope(idx, total, slope_tif, work_dir, save_dir,slope_threshold, min_area, max_area,max_width,process_num):
    '''
    extract RTS headwall from slope
    :param idx:
    :param slope_tif:
    :param slope_threshold:
    :param area_threshold:
    :return:
    '''

    headwall_shp = os.path.splitext(os.path.basename(io_function.get_name_by_adding_tail(slope_tif,'headwall')))[0] + '.shp'
    save_headwall_shp = os.path.join(save_dir,headwall_shp)
    if os.path.isfile(save_headwall_shp):
        print('%s exists, skip'%save_headwall_shp)
        return save_headwall_shp


    print('(%d/%d) extracting headwall from %s'%(idx,total,slope_tif))

    # binary slope
    slope_bin_path = os.path.join(work_dir, os.path.basename(io_function.get_name_by_adding_tail(slope_tif, 'bin')))
    if os.path.isfile(slope_bin_path):
        print('%s exist'%slope_bin_path)
    else:
        slope_data, nodata = raster_io.read_raster_one_band_np(slope_tif)
        bin_slope = np.zeros_like(slope_data,dtype=np.uint8)
        bin_slope[slope_data > slope_threshold] = 1

        # Dilation or opening
        # https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_morphological_ops/py_morphological_ops.html
        kernel = np.ones((3, 3), np.uint8)  # if kernal is 5 or larger, will remove some narrow parts.
        # bin_slope = cv2.dilate(bin_slope,kernel,iterations = 1)
        bin_slope = cv2.morphologyEx(bin_slope, cv2.MORPH_OPEN, kernel)     # use opening to remove some noise
        # bin_slope = cv2.morphologyEx(bin_slope, cv2.MORPH_CLOSE, kernel)    # closing small holes inside

        # save
        slope_bin = bin_slope*255
        raster_io.save_numpy_array_to_rasterfile(slope_bin,slope_bin_path,slope_tif,nodata=0)   # set nodata as 0

    # to shapefile
    slope_bin_shp = vector_gpd.raster2shapefile(slope_bin_path,connect8=True)
    if slope_bin_shp is None:
        return False

    # only keep small and narrow ones
    polygons = vector_gpd.read_polygons_gpd(slope_bin_shp,b_fix_invalid_polygon=False)

    remain_polygons = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, poly in enumerate(polygons):
        # remove quite large or too small ones
        if poly.area > max_area or poly.area < min_area:
            remove_count += 1
            continue
        remain_polygons.append(poly)

    basic.outputlogMessage('remove %d polygons based on area, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons), save_headwall_shp))

    polyons_noMulti = [vector_gpd.MultiPolygon_to_polygons(idx, poly) for idx, poly in enumerate(remain_polygons)]
    remain_polygons = []
    for polys in polyons_noMulti:
        polys = [poly for poly in polys if poly.area > min_area]  # remove tiny polygon before buffer
        remain_polygons.extend(polys)
    print('convert MultiPolygon to polygons, remain %d' % (len(remain_polygons)))

    if len(remain_polygons) < 1:
        return False

    wkt = map_projection.get_raster_or_vector_srs_info_wkt(slope_tif)
    save_pd = pd.DataFrame({'Polygon':remain_polygons})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,save_headwall_shp)

    # add some shape info
    shape_info_list = [vector_gpd.calculate_polygon_shape_info(item) for item in remain_polygons]
    shapeinfo_all_dict = vector_gpd.list_to_dict(shape_info_list)
    vector_gpd.add_attributes_to_shp(save_headwall_shp, shapeinfo_all_dict)

    # calculate width based on medial axis
    calculate_distance_medial_axis(save_headwall_shp,process_num=process_num)

    return save_headwall_shp

def calculate_distance_medial_axis(input_shp, process_num=4):
    print('calculating polygon width based on medial axis')

    code_dir = os.path.expanduser('~/codes/PycharmProjects/ChangeDet_DL/thawSlumpChangeDet')
    sys.path.insert(0, code_dir)
    # calculate width based on expanding areas
    import cal_retreat_rate
    if cal_retreat_rate.cal_expand_area_distance(input_shp,proc_num=process_num):
        return True

def test_calculate_distance_medial_axis():

    # save polygons without holes
    shp = os.path.join('dem_headwall_shp','slope_sub_headwall.shp')

    polygons = vector_gpd.read_polygons_gpd(shp)
    polygon_nohole = [ vector_gpd.fill_holes_in_a_polygon(item) for item in polygons]

    wkt = map_projection.get_raster_or_vector_srs_info_wkt(shp)
    save_nohole_shp = io_function.get_name_by_adding_tail(shp,'nohole')
    save_pd = pd.DataFrame({'Polygon':polygon_nohole})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,save_nohole_shp)

    #
    calculate_distance_medial_axis(save_nohole_shp, process_num=16)



# def test_extract_headwall_from_slope():
#     print('\n')
#     slope = os.path.expanduser('~/Data/tmp_data/slope_sub.tif')
#     working_dir = './'
#     save_dir = dem_headwall_shp_dir
#     if os.path.isdir(working_dir) is False:
#         io_function.mkdir(working_dir)
#     if os.path.isdir(save_dir) is False:
#         io_function.mkdir(save_dir)
#
#     min_slope = 20
#     min_size = 400
#     max_size = 50000
#     max_width = 50
#     process_num = 10
#
#     extract_headwall_from_slope(0, 1, slope, working_dir, save_dir, min_slope, min_size, max_size,max_width,process_num)


def main(options, args):
    input = args[0]

    if input.endswith('.txt'):
        slope_tifs = io_function.read_list_from_txt(input)
    elif os.path.isdir(input):
        slope_tifs = io_function.get_file_list_by_ext('.tif',input, bsub_folder=True)
    else:
        slope_tifs = [ input]

    working_dir = './'
    save_dir = dem_headwall_shp_dir
    if os.path.isdir(working_dir) is False:
        io_function.mkdir(working_dir)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    failed_tifs = []

    min_slope = 20
    max_size = 50000
    min_size = 400
    for idx, slope in enumerate(slope_tifs):
        if extract_headwall_from_slope(idx, len(slope_tifs), slope,working_dir,save_dir, min_slope,min_size,max_size) is False:
            failed_tifs.append(slope)

    io_function.save_list_to_txt('extract_headwall_failed_tifs.txt',failed_tifs)



if __name__ == '__main__':
    usage = "usage: %prog [options] slopefile or slopefile_list_txt or dir "
    parser = OptionParser(usage=usage, version="1.0 2021-4-30")
    parser.description = 'Introduction: extract RTS headwall from slope derived from ArcticDEM '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")



    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

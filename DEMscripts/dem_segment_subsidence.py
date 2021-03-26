#!/usr/bin/env python
# Filename: watershed_segment 
"""
introduction: segment a grey image

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 1 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import vector_features
import raster_statistic

import cv2
import numpy as np
import pandas as pd

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0,code_dir)
sys.path.insert(0,os.path.join(code_dir,'tools'))   # for some modules in this folder
from tools.grey_image_segment import segment_a_grey_image

def post_processing_subsidence(in_shp):
    polygons = vector_gpd.read_polygons_gpd(in_shp)

    # get shapeinfo
    # poly_shapeinfo_list = []
    save_polyons = []
    for poly in polygons:
        # get INarea, INperimete, WIDTH, HEIGHT, ratio_w_h, hole_count
        # shapeinfo = vector_gpd.calculate_polygon_shape_info(poly)     # error:  'MultiPolygon' object has no attribute 'interiors'
        # poly_shapeinfo_list.append(shapeinfo)
        # if shapeinfo['INarea'] < 40:    # remove the one with area smaller than 40 m^2
        if poly.area < 90:    # remove the one with area smaller than 40 m^2
            continue
        save_polyons.append(poly)

    save_pd = pd.DataFrame({'Polygon': save_polyons})
    wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)
    save_shp = io_function.get_name_by_adding_tail(in_shp,'post')
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,save_shp)




def segment_subsidence_on_dem_diff(dem_diff_tif, save_dir):

    out_pre = os.path.splitext(os.path.basename(dem_diff_tif))[0]

    # read images
    one_band_img, nodata = raster_io.read_raster_one_band_np(dem_diff_tif)

    # segmentation by threshold (may have too many noise)
    # mean = np.nanmean(one_band_img)
    # print("mean value is: %.4f"%mean)
    # one_band_img = one_band_img - mean    # cannot use mean which may affect by some Outliers
    out_labels = np.zeros_like(one_band_img,dtype=np.uint8)
    out_labels[ one_band_img < -2 ] = 1     # end in a lot of noise, change to -2, -1 results in a lot of polygons

    # apply median filter
    out_labels = cv2.medianBlur(out_labels, 3)  # with kernal=3

    # save the label
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    label_path = os.path.join(save_dir, out_pre + '_label.tif')
    raster_io.save_numpy_array_to_rasterfile(out_labels, label_path, dem_diff_tif, nodata=0)

    # convert the label to shapefile
    out_shp = os.path.join(save_dir, out_pre + '.shp')
    command_string = 'gdal_polygonize.py -8 %s -b 1 -f "ESRI Shapefile" %s' % (label_path, out_shp)
    res = os.system(command_string)
    if res != 0:
        sys.exit(1)

    # post-processing
    post_processing_subsidence(out_shp)

def get_mean_from_array(in_array, nodata,range=None):
    data_1d = in_array.flatten()
    data_1d = data_1d[ data_1d != nodata]
    data_1d = data_1d[~np.isnan(data_1d)]  # remove nan value
    value = np.mean(data_1d)
    return value

def get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-0.5, min_area=40, max_area=100000000, process_num=1):

    save_shp = io_function.get_name_by_adding_tail(in_shp, 'post')
    if os.path.isfile(save_shp):
        basic.outputlogMessage('%s exists, skip'%save_shp)
        return save_shp

    shp_pre = os.path.splitext(os.path.basename(in_shp))[0]
    # read polygons and label from segment algorithm, note: some polygons may have the same label
    polygons, demD_mean_list = vector_gpd.read_polygons_attributes_list(in_shp,'demD_mean')
    print('Read %d polygons'%len(polygons))


    demD_height, demD_width, demD_band_num, demD_date_type = raster_io.get_height_width_bandnum_dtype(dem_diff_tif)
    # print(demD_date_type)

    # # read mean elevation difference
    # attributes_path = os.path.join(os.path.dirname(in_shp), shp_pre + '_attributes.txt')
    #
    # # for each seg lable [mean, std, pixel count], if dem_diff_tif is float 32, then in meters, if int16, then in centimeter
    # poly_attributes = io_function.read_dict_from_txt_json(attributes_path)

    # if int16, then it's in centimeter
    if demD_date_type == 'int16':
        dem_diff_thread_m = dem_diff_thread_m*100

    remain_polyons = []
    rm_min_area_count = 0
    rm_diff_thr_count = 0
    for poly, demD_mean in zip(polygons, demD_mean_list):
        if poly.area < min_area:
            rm_min_area_count += 1
            continue
        # mean value: not subsidence
        if demD_mean > dem_diff_thread_m:  #
            rm_diff_thr_count += 1
            continue

        remain_polyons.append(poly)

    print('remove %d polygons based on min_area, %d polygons based on dem_diff_threshold, remain %d ones'%(rm_min_area_count, rm_diff_thr_count,len(remain_polyons)))

    if len(remain_polyons) > 1:
        # we should only merge polygon with similar reduction, but we already remove polygons with mean reduction > threshhold
        # merge touch polygons
        print(timeTools.get_now_time_str(), 'start building adjacent_matrix')
        # adjacent_matrix = vector_features.build_adjacent_map_of_polygons(remain_polyons)
        adjacent_matrix = vector_gpd.build_adjacent_map_of_polygons(remain_polyons, process_num=process_num)
        print(timeTools.get_now_time_str(), 'finish building adjacent_matrix')

        if adjacent_matrix is False:
            return False
        merged_polygons = vector_features.merge_touched_polygons(remain_polyons,adjacent_matrix)
        print(timeTools.get_now_time_str(), 'finish merging touched polygons, get %d ones'%(len(merged_polygons)))

        # remove large ones
        remain_polyons = []
        rm_max_area_count = 0
        for poly in merged_polygons:
            if poly.area > max_area:
                rm_max_area_count += 1
                continue
            remain_polyons.append(poly)

        print('remove %d polygons based on max_area, remain %d'%(rm_max_area_count, len(remain_polyons)))

    # calcualte attributes of remain ones: area, dem_diff: mean, std
    poly_areas = [ poly.area for poly in remain_polyons ]
    poly_ids = [ item+1  for item in range(len(remain_polyons)) ]

    save_pd = pd.DataFrame({'poly_ids':poly_ids,'poly_area':poly_areas,'Polygon': remain_polyons})
    wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)

    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,save_shp)

    raster_statistic.zonal_stats_multiRasters(save_shp,dem_diff_tif,stats=['mean','std'],prefix='demD')

    return save_shp



def segment_subsidence_grey_image(dem_diff_grey_8bit, dem_diff, save_dir,process_num, subsidence_thr_m=-0.5, min_area=40, max_area=100000000):
    '''
    segment subsidence areas based on 8bit dem difference
    :param dem_diff_grey_8bit:
    :param dem_diff:
    :param save_dir:
    :param process_num:
    :param subsidence_thr_m: mean value less than this one consider as subsidence (in meter)
    :param min_area: min size in m^2 (defualt is 40 m^2, 10 pixels on ArcticDEM)
    :param max_area: min size in m^2 (default is 10km by 10 km)
    :return:
    '''

    io_function.is_file_exist(dem_diff_grey_8bit)

    # get initial polygons
    # because the label from segmentation for superpixels are not unique, so we may need to get mean dem diff based on polygons, set org_raster=None
    segment_shp_path = segment_a_grey_image(dem_diff_grey_8bit,save_dir,process_num, org_raster=None)

    # get dem elevation information for each polygon
    raster_statistic.zonal_stats_multiRasters(segment_shp_path, dem_diff, stats=['mean', 'std'], prefix='demD',process_num=process_num)

    # get DEM diff information for each polygon.
    dem_diff_shp = get_dem_subscidence_polygons(segment_shp_path, dem_diff, dem_diff_thread_m=subsidence_thr_m,
                                 min_area=min_area, max_area=max_area, process_num=process_num)

    basic.outputlogMessage('obtain elevation reduction polygons: %s'%dem_diff_shp)
    return True


def test_get_dem_subscidence_polygons():
    in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_parallel_sub/WR_dem_diff_DEM_diff_prj_8bit_sub.shp')
    dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/WR_dem_diff_DEM_diff_prj.tif')
    get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-1, min_area=40, max_area=100000000)



def main(options, args):

    dem_diff_path = args[0]
    io_function.is_file_exist(dem_diff_path)
    save_dir = options.save_dir
    process_num = options.process_num
    # segment_subsidence_on_dem_diff(dem_diff_path,save_dir)

    dem_diff_grey_8bit = options.dem_diff_8bit
    ele_diff_thr = options.ele_diff_thr
    min_area = options.min_area
    max_area = options.max_area

    segment_subsidence_grey_image(dem_diff_grey_8bit, dem_diff_path, save_dir, process_num, subsidence_thr_m=ele_diff_thr,
                                  min_area=min_area, max_area=max_area)

    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] dem_diff_tif "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment subsidence based on DEM difference  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save results")

    parser.add_option("-g", "--dem_diff_8bit",
                      action="store", dest="dem_diff_8bit",
                      help="the normalized elevation difference with 8bit")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes for image segmentation")

    parser.add_option("-t", "--ele_diff_thr",
                      action="store", dest="ele_diff_thr", type=float, default=-0.5,
                      help="the threshold for elevation reduction ")

    parser.add_option("-l", "--min_area",
                      action="store", dest="min_area", type=float, default=40,
                      help="the minimum area for each elevation reduction polygon")

    parser.add_option("-u", "--max_area",
                      action="store", dest="max_area", type=float, default=100000000,  # 10 km by 10 km
                      help="the maximum area for each elevation reduction polygon")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

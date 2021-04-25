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
import re

import dem_common

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0,code_dir)
sys.path.insert(0,os.path.join(code_dir,'tools'))   # for some modules in this folder
from tools.grey_image_segment import segment_a_grey_image
from tools.seg_polygonize_cal_attributes import polygonize_label_images

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

def get_dem_diff_8bit(dem_diff_path):
    # find 8bit one
    tif_8bit = io_function.get_name_by_adding_tail(dem_diff_path, '8bit')
    demD_8bit= os.path.join(dem_common.grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
    if os.path.isfile(demD_8bit) is False:
        basic.outputlogMessage('error, 8bit DEM diff not exists: %s '%demD_8bit)
        return None
    return demD_8bit

def get_save_dir(dem_diff_path):

    grid_id = int(re.findall('grid\d+',os.path.basename(dem_diff_path))[0][4:])
    save_dir = os.path.join(dem_common.grid_dem_diffs_segment_dir, 'segment_result_grid%d'%grid_id)
    return save_dir

def get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-0.5, min_area=40, max_area=100000000, process_num=1):

    save_shp = io_function.get_name_by_adding_tail(in_shp, 'post')
    if os.path.isfile(save_shp):
        basic.outputlogMessage('%s exists, skip'%save_shp)
        return save_shp

    shp_pre = os.path.splitext(os.path.basename(in_shp))[0]
    # read polygons and label from segment algorithm, note: some polygons may have the same label
    polygons, demD_mean_list = vector_gpd.read_polygons_attributes_list(in_shp,'demD_mean')
    print('Read %d polygons'%len(polygons))
    if demD_mean_list is None:
        raise ValueError('demD_mean not in %s, need to remove it and then re-create'%in_shp)


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
        machine_name = os.uname()[1]
        if 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:
            print('Warning, some problem of parallel running in build_adjacent_map_of_polygons on curc, but ok in my laptop and uist, change process_num = 1')
            process_num = 1
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

    wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)

    polyons_noMulti = [ vector_gpd.MultiPolygon_to_polygons(idx,poly) for idx,poly in enumerate(remain_polyons) ]
    remain_polyons = []
    for polys in polyons_noMulti:
        polys = [poly for poly in polys if poly.area > min_area]    # remove tiny polygon before buffer
        remain_polyons.extend(polys)
    print('convert MultiPolygon to polygons, remain %d' % (len(remain_polyons)))

    if len(remain_polyons) < 1:
        return None

    # based on the merged polygons, calculate the mean dem diff, relative dem_diff
    buffer_surrounding = 20  # meters
    surrounding_polygons = vector_gpd.get_surrounding_polygons(remain_polyons,buffer_surrounding)
    surrounding_shp = io_function.get_name_by_adding_tail(in_shp, 'surrounding')
    surr_pd = pd.DataFrame({'Polygon': surrounding_polygons})
    vector_gpd.save_polygons_to_files(surr_pd, 'Polygon', wkt, surrounding_shp)
    raster_statistic.zonal_stats_multiRasters(surrounding_shp, dem_diff_tif, stats=['mean', 'std', 'count'], prefix='demD',process_num=process_num)


    # calcualte attributes of remain ones: area, dem_diff: mean, std
    merged_pd = pd.DataFrame({'Polygon': remain_polyons})
    merged_shp = io_function.get_name_by_adding_tail(in_shp, 'merged')
    vector_gpd.save_polygons_to_files(merged_pd, 'Polygon', wkt, merged_shp)
    raster_statistic.zonal_stats_multiRasters(merged_shp, dem_diff_tif, stats=['mean','std','count'], prefix='demD', process_num=process_num)

    # calculate the relative dem diff
    surr_dem_diff_list = vector_gpd.read_attribute_values_list(surrounding_shp,'demD_mean')
    merge_poly_dem_diff_list = vector_gpd.read_attribute_values_list(merged_shp,'demD_mean')
    if len(surr_dem_diff_list) != len(merge_poly_dem_diff_list):
        raise ValueError('The number of surr_dem_diff_list and merge_poly_dem_diff_list is different')
    relative_dem_diff_list = [  mer - sur for sur, mer in zip(surr_dem_diff_list, merge_poly_dem_diff_list) ]

    merge_poly_demD_std_list = vector_gpd.read_attribute_values_list(merged_shp,'demD_std')
    merge_poly_demD_count_list = vector_gpd.read_attribute_values_list(merged_shp,'demD_count')

    # remove large ones
    save_polyons = []
    save_demD_mean_list = []
    save_demD_std_list = []
    save_demD_count_list = []
    save_rel_diff_list = []
    save_surr_demD_list = []
    rm_rel_dem_diff_count = 0
    rm_min_area_count = 0
    for idx in range(len(remain_polyons)):
        # relative dem diff
        if relative_dem_diff_list[idx] > dem_diff_thread_m:  #
            rm_rel_dem_diff_count += 1
            continue

        # when convert MultiPolygon to Polygon, may create some small polygons
        if remain_polyons[idx].area < min_area:
            rm_min_area_count += 1
            continue


        save_polyons.append(remain_polyons[idx])
        save_demD_mean_list.append(merge_poly_dem_diff_list[idx])
        save_demD_std_list.append(merge_poly_demD_std_list[idx])
        save_demD_count_list.append(merge_poly_demD_count_list[idx])
        save_rel_diff_list.append(relative_dem_diff_list[idx])
        save_surr_demD_list.append(surr_dem_diff_list[idx])

    print('remove %d polygons based on relative rel_demD and %d based on min_area, remain %d' % (rm_rel_dem_diff_count, rm_min_area_count, len(save_polyons)))

    if len(save_polyons) < 1:
        print('Warning, no polygons after remove based on relative demD')
        return None

    poly_ids = [ item+1  for item in range(len(save_polyons)) ]
    poly_areas = [poly.area for poly in save_polyons]

    save_pd = pd.DataFrame({'poly_id':poly_ids, 'poly_area':poly_areas,'demD_mean':save_demD_mean_list, 'demD_std':save_demD_std_list,
                             'demD_count':save_demD_count_list, 'surr_demD':save_surr_demD_list, 'rel_demD':save_rel_diff_list ,'Polygon': save_polyons})

    vector_gpd.save_polygons_to_files(save_pd, 'Polygon', wkt, save_shp)

    # add date difference if there are available
    date_diff_base = os.path.basename(dem_diff_tif).replace('DEM_diff','date_diff')
    date_diff_tif = os.path.join(os.path.dirname(dem_diff_tif) , date_diff_base)
    if os.path.isfile(date_diff_tif):
        raster_statistic.zonal_stats_multiRasters(save_shp, date_diff_tif, stats=['mean', 'std'], prefix='dateD',
                                              process_num=process_num)

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

    out_pre = os.path.splitext(os.path.basename(dem_diff_grey_8bit))[0]
    segment_shp_path = os.path.join(save_dir, out_pre + '.shp')

    # get initial polygons
    # because the label from segmentation for superpixels are not unique, so we may need to get mean dem diff based on polygons, set org_raster=None
    label_path = segment_a_grey_image(dem_diff_grey_8bit,save_dir,process_num, org_raster=None)

    if os.path.isfile(segment_shp_path) and vector_gpd.is_field_name_in_shp(segment_shp_path,'demD_mean'):
        basic.outputlogMessage('%s exists, skip'%segment_shp_path)
    else:

        # remove segment_shp_path if it exist, but don't have demD_mean
        if os.path.isfile(segment_shp_path):
            io_function.delete_shape_file(segment_shp_path)

        # remove nodato (it was copy from the input image)
        command_str = 'gdal_edit.py -unsetnodata ' + label_path
        basic.os_system_exit_code(command_str)

        # convert the label to shapefile # remove -8 (to use 4 connectedness.)
        command_string = 'gdal_polygonize.py %s -b 1 -f "ESRI Shapefile" %s' % (label_path, segment_shp_path)
        res = os.system(command_string)
        if res != 0:
            sys.exit(1)

        # get dem elevation information for each polygon
        raster_statistic.zonal_stats_multiRasters(segment_shp_path, dem_diff, stats=['mean', 'std','count'], prefix='demD',process_num=process_num)

    # get DEM diff information for each polygon.
    dem_diff_shp = get_dem_subscidence_polygons(segment_shp_path, dem_diff, dem_diff_thread_m=subsidence_thr_m,
                                 min_area=min_area, max_area=max_area, process_num=process_num)

    basic.outputlogMessage('obtain elevation reduction polygons: %s'%dem_diff_shp)
    return True

def segment_subsidence_grey_image_v2(dem_diff_grey_8bit, dem_diff, save_dir,process_num, subsidence_thr_m=-0.5, min_area=40, max_area=100000000):
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

    out_pre = os.path.splitext(os.path.basename(dem_diff_grey_8bit))[0]
    segment_shp_path = os.path.join(save_dir, out_pre + '.shp')

    # get initial polygons
    # because the label from segmentation for superpixels are not unique, so we may need to get mean dem diff based on polygons, set org_raster=None
    label_path_list = segment_a_grey_image(dem_diff_grey_8bit,save_dir,process_num, org_raster=None,b_save_patch_label=True)

    patch_shp_list = polygonize_label_images(label_path_list, org_raster=dem_diff, stats=['mean', 'std', 'count'], prefix='demD',
                            process_num=process_num, b_remove_nodata=True)

    # post-processing for each patch shp
    post_patch_shp_list = []
    for idx, shp in enumerate(patch_shp_list):
        # get DEM diff information for each polygon.
        post_shp = get_dem_subscidence_polygons(shp, dem_diff, dem_diff_thread_m=subsidence_thr_m,
                                     min_area=min_area, max_area=max_area, process_num=1)
        if post_shp is not None:
            post_patch_shp_list.append(post_shp)

    # merge shapefile
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    vector_gpd.merge_shape_files(post_patch_shp_list,segment_shp_path)

    # post-processing again
    dem_diff_shp = get_dem_subscidence_polygons(segment_shp_path, dem_diff, dem_diff_thread_m=subsidence_thr_m,
                                            min_area=min_area, max_area=max_area, process_num=1)

    basic.outputlogMessage('obtain elevation reduction polygons: %s'%dem_diff_shp)
    return True


def test_get_dem_subscidence_polygons():
    # in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_parallel_sub/WR_dem_diff_DEM_diff_prj_8bit_sub.shp')
    # dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/WR_dem_diff_DEM_diff_prj.tif')
    # get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-1, min_area=40, max_area=100000000)

    in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_result_grid8868/ala_north_slo_extent_latlon_grid_ids_DEM_diff_grid8868_8bit.shp')
    dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/ala_north_slo_extent_latlon_grid_ids_DEM_diff_grid8868.tif')
    get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-2, min_area=40, max_area=100000000, process_num=10)

def main(options, args):

    dem_diff_path = args[0]
    io_function.is_file_exist(dem_diff_path)
    save_dir = options.save_dir
    process_num = options.process_num
    # segment_subsidence_on_dem_diff(dem_diff_path,save_dir)

    dem_diff_grey_8bit = options.dem_diff_8bit
    if dem_diff_grey_8bit is None:
        dem_diff_grey_8bit  = get_dem_diff_8bit(dem_diff_path)

    if save_dir is None:
        save_dir = get_save_dir(dem_diff_path)

    ele_diff_thr = options.ele_diff_thr
    min_area = options.min_area
    max_area = options.max_area

    segment_subsidence_grey_image(dem_diff_grey_8bit, dem_diff_path, save_dir, process_num, subsidence_thr_m=ele_diff_thr,
                                  min_area=min_area, max_area=max_area)

    # the resut is a litte worse
    # segment_subsidence_grey_image_v2(dem_diff_grey_8bit, dem_diff_path, save_dir, process_num, subsidence_thr_m=ele_diff_thr,
    #                               min_area=min_area, max_area=max_area)

    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] dem_diff_tif "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment subsidence based on DEM difference  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",
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

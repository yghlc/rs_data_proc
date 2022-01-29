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
import geopandas as gpd
import re

import dem_common
from dem_common import process_log_dir,grid_no_subscidence_poly_txt
from dem_diff_to_8bit import one_dem_diff_to_8bit

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0,code_dir)
sys.path.insert(0,os.path.join(code_dir,'tools'))   # for some modules in this folder
from tools.grey_image_segment import segment_a_grey_image
from tools.seg_polygonize_cal_attributes import polygonize_label_images

tile_min_overlap = 1

def save_id_grid_no_subsidence(grid_id):
    # grid_no_headwall_txt
    if os.path.isdir(process_log_dir) is False:
        io_function.mkdir(process_log_dir)

    id_list = []
    if os.path.isfile(grid_no_subscidence_poly_txt):
        id_list = io_function.read_list_from_txt(grid_no_subscidence_poly_txt)    # no need covert to int
    id_str = str(grid_id)
    if id_str in id_list:
        return True
    else:
        # save by adding one line
        with open(grid_no_subscidence_poly_txt,'a') as f_obj:
            f_obj.writelines(str(grid_id) + '\n')
        return True


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
    # if not exist, try to generate it
    if os.path.isfile(demD_8bit) is False:
        if  one_dem_diff_to_8bit(dem_diff_path) is False:
            basic.outputlogMessage('error, 8bit DEM diff not exists: %s '%demD_8bit)
            return None
    return demD_8bit

def get_save_dir(dem_diff_path):

    grid_id = int(re.findall('grid\d+',os.path.basename(dem_diff_path))[0][4:])
    save_dir = os.path.join(dem_common.grid_dem_diffs_segment_dir, 'segment_result_grid%d'%grid_id)
    return save_dir

def merge_polygons_patchBYpatch(patch_DN_range_txt, in_polygons, polygon_DNs, process_num=1):
    patch_DN_range_list = [ int(item) for item in io_function.read_list_from_txt(patch_DN_range_txt)]
    # print(patch_DN_range_list)
    # divide polygons
    patch_polygons = {}
    range_count = len(patch_DN_range_list) - 1
    for idx in range(range_count):
        patch_polygons[idx] = []
    for poly, dn in zip(in_polygons, polygon_DNs):
        for idx in range(range_count):
            if dn > patch_DN_range_list[idx] and dn <= patch_DN_range_list[idx+1]:
                patch_polygons[idx].append(poly)

    # merge polygon patch by patch
    polygons_patch_merge = []
    for idx in range(range_count):
        print(timeTools.get_now_time_str(), 'will merge %d polygons for %d th patch'%(len(patch_polygons[idx]),idx))
        if len(patch_polygons[idx]) < 2:
            polygons_patch_merge.extend(patch_polygons[idx])
            continue
        adjacent_matrix = vector_gpd.build_adjacent_map_of_polygons(patch_polygons[idx], process_num=process_num)
        if adjacent_matrix is False:
            polygons_patch_merge.extend(patch_polygons[idx])
            continue
        merged_polygons = vector_features.merge_touched_polygons(patch_polygons[idx], adjacent_matrix)
        polygons_patch_merge.extend(merged_polygons)

    # merge polygon all
    print(timeTools.get_now_time_str(), 'will merge %d polygons for the entire raster' % (len(polygons_patch_merge)))
    adjacent_matrix = vector_gpd.build_adjacent_map_of_polygons(polygons_patch_merge, process_num=process_num)
    if adjacent_matrix is False:
        return polygons_patch_merge
    last_merged_polygons = vector_features.merge_touched_polygons(polygons_patch_merge, adjacent_matrix)
    return last_merged_polygons

def merge_polygon_rasterize(ref_raster, in_polygons):

    # rasterize to raster
    save_raster = os.path.basename(io_function.get_name_by_adding_tail(ref_raster,'merge'))
    raster_io.burn_polygons_to_a_raster(ref_raster,in_polygons,1,save_raster)

    # set nodata
    if raster_io.set_nodata_to_raster_metadata(save_raster,0) is False:
        raise IOError('Set nodata failed for %s'%save_raster)

    # polygonize
    out_shp = vector_gpd.raster2shapefile(save_raster, connect8=True)
    if out_shp is None:
        raise IOError('polygonzied failed for %s' % save_raster)

    # read polygon
    merge_polygons = vector_gpd.read_polygons_gpd(out_shp,b_fix_invalid_polygon=False)

    print(timeTools.get_now_time_str(),'Get %d merged polygons'%len(merge_polygons))
    return merge_polygons


def filter_merge_polygons(in_shp,merged_shp,wkt, min_area,max_area,dem_diff_tif,dem_diff_thread_m,process_num):

    if os.path.isfile(merged_shp):
        # also check the file is complete
        polys, demD_values = vector_gpd.read_polygons_attributes_list(merged_shp,'demD_mean')
        if len(polys) < 1 or demD_values is None or len(demD_values) < 1:
            basic.outputlogMessage('%s already exists, but not complete, will be overwritten'%merged_shp)
        else:
            basic.outputlogMessage('%s exists, skip'%merged_shp)
            return merged_shp

    # read polygons and label from segment algorithm, note: some polygons may have the same label
    # polygons, demD_mean_list = vector_gpd.read_polygons_attributes_list(in_shp,'demD_mean')
    polygons, attributes = vector_gpd.read_polygons_attributes_list(in_shp,['demD_mean','DN'])
    demD_mean_list = attributes[0]
    DN_list = attributes[1]
    print('Read %d polygons'%len(polygons))
    if demD_mean_list is None:
        raise ValueError('demD_mean not in %s, need to remove it and then re-create'%in_shp)

    # replace None (if exists) as nan
    demD_mean_list = np.array(demD_mean_list, dtype=float)

    # replace nan values as 0
    demD_mean_list = np.nan_to_num(demD_mean_list)

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
    if len(remain_polyons) < 1:
        return None


    # we should only merge polygon with similar reduction, but we already remove polygons with mean reduction > threshhold
    # merge touch polygons
    # print(timeTools.get_now_time_str(), 'start building adjacent_matrix')
    # # adjacent_matrix = vector_features.build_adjacent_map_of_polygons(remain_polyons)
    # machine_name = os.uname()[1]
    # if 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:
    #     print('Warning, some problem of parallel running in build_adjacent_map_of_polygons on curc, '
    #           'but ok in my laptop and uist, change process_num = 1')
    #     process_num = 1
    ############################################################
    ## build adjacent_matrix then merge for entire raster
    # adjacent_matrix = vector_gpd.build_adjacent_map_of_polygons(remain_polyons, process_num=process_num)
    # print(timeTools.get_now_time_str(), 'finish building adjacent_matrix')
    #
    # if adjacent_matrix is False:
    #     return None
    # merged_polygons = vector_features.merge_touched_polygons(remain_polyons, adjacent_matrix)

    ############################################################
    # ## build adjacent_matrix then merge, patch by patch (not too many improvements)
    # label_id_range_txt = os.path.splitext(in_shp)[0] + '_label_IDrange.txt'
    # merged_polygons = merge_polygons_patchBYpatch(label_id_range_txt, remain_polyons, DN_list, process_num=process_num)

    ############################################################
    ## merge polygons using rasterize
    label_raster = os.path.splitext(in_shp)[0] + '_label.tif'
    merged_polygons = merge_polygon_rasterize(label_raster, remain_polyons)

    print(timeTools.get_now_time_str(), 'finish merging touched polygons, get %d ones' % (len(merged_polygons)))

    # remove large ones
    remain_polyons = []
    rm_max_area_count = 0
    for poly in merged_polygons:
        if poly.area > max_area:
            rm_max_area_count += 1
            continue
        remain_polyons.append(poly)

    print('remove %d polygons based on max_area, remain %d' % (rm_max_area_count, len(remain_polyons)))

    polyons_noMulti = [vector_gpd.MultiPolygon_to_polygons(idx, poly) for idx, poly in enumerate(remain_polyons)]
    remain_polyons = []
    for polys in polyons_noMulti:
        polys = [poly for poly in polys if poly.area > min_area]  # remove tiny polygon before buffer
        remain_polyons.extend(polys)
    print('convert MultiPolygon (filter_merge_polygons) to polygons and remove small ones, remain %d' % (len(remain_polyons)))

    if len(remain_polyons) < 1:
        return None

    # calcualte attributes of remain ones: area, dem_diff: mean, std
    merged_pd = pd.DataFrame({'Polygon': remain_polyons})
    vector_gpd.save_polygons_to_files(merged_pd, 'Polygon', wkt, merged_shp)

    # based on the merged polygons, calculate the mean dem diff
    raster_statistic.zonal_stats_multiRasters(merged_shp, dem_diff_tif, tile_min_overlap=tile_min_overlap,
                                              stats=['mean', 'std', 'count'], prefix='demD',process_num=process_num)

    return merged_shp

def get_surrounding_polygons(remain_polyons,surrounding_shp,wkt, dem_diff_tif,buffer_surrounding,process_num):
    if os.path.isfile(surrounding_shp):
        # also check the file is complete
        surr_polys, surr_demD = vector_gpd.read_polygons_attributes_list(surrounding_shp,'demD_mean')
        if len(surr_polys) < len(remain_polyons) or surr_demD is None or len(surr_demD) < len(remain_polyons):
            basic.outputlogMessage('%s already exists, but not complete, will be overwritten'%surrounding_shp)
        else:
            basic.outputlogMessage('%s already exists, skip'%surrounding_shp)
            return surrounding_shp

    # based on the merged polygons, calculate the relative dem_diff
    surrounding_polygons = vector_gpd.get_surrounding_polygons(remain_polyons, buffer_surrounding)
    surr_pd = pd.DataFrame({'Polygon': surrounding_polygons})
    vector_gpd.save_polygons_to_files(surr_pd, 'Polygon', wkt, surrounding_shp)
    raster_statistic.zonal_stats_multiRasters(surrounding_shp, dem_diff_tif, tile_min_overlap=tile_min_overlap,
                                              stats=['mean', 'std', 'count'],prefix='demD', process_num=process_num)
    return surrounding_shp

def remove_polygons_based_relative_dem_diff(remain_polyons,merged_shp,surrounding_shp,wkt, save_shp, min_area, dem_diff_thread_m):
    if os.path.isfile(save_shp):
        basic.outputlogMessage('%s exists, skip'%save_shp)
        return save_shp

    # calculate the relative dem diff
    surr_dem_diff_list = vector_gpd.read_attribute_values_list(surrounding_shp,'demD_mean')
    merge_poly_dem_diff_list = vector_gpd.read_attribute_values_list(merged_shp,'demD_mean')
    # convert to float type (can change None to nan)
    surr_dem_diff_list = np.array(surr_dem_diff_list, dtype=float)
    merge_poly_dem_diff_list = np.array(merge_poly_dem_diff_list, dtype=float)

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

        # when convert MultiPolygon to Polygon, may create some small polygons (in function merge shp)
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

    return save_shp

def remove_polygons_based_shapeinfo(in_shp, output_shp,area_limit,circularit_limit, holes_count ):
    # remove polygons if they are larege (> area_limit) and narrow (<circularit_limit)
    # if too many holse, may consider remove them as well.

    # remove based on: INarea, INperimete,circularit,ratio_w_h, hole_count

    if os.path.isfile(output_shp):
        basic.outputlogMessage('%s exists, skip'%output_shp)
        return output_shp

    polygons = vector_gpd.read_polygons_gpd(in_shp)

    # add some shape info
    shape_info_list = [ vector_gpd.calculate_polygon_shape_info(item)  for item in polygons]
    shapeinfo_all_dict = vector_gpd.list_to_dict(shape_info_list)
    vector_gpd.add_attributes_to_shp(in_shp,shapeinfo_all_dict)

    shapefile = gpd.read_file(in_shp)

    # remove relative large but narrow ones.
    remove_count = 0
    for idx, row in shapefile.iterrows():
        shape_info = shape_info_list[idx]

        # remove quite large but narrow ones
        if shape_info['INarea'] > area_limit and shape_info['circularit'] < circularit_limit :
            shapefile.drop(idx, inplace=True)
            remove_count += 1
            continue

        # remove holes
        if shape_info['hole_count'] > holes_count:
            shapefile.drop(idx, inplace=True)
            remove_count += 1
            continue

    basic.outputlogMessage('remove %d polygons based on shapeinfo, remain %d ones saving to %s' %
                           (remove_count, len(shapefile.geometry.values), output_shp))
    # save results
    shapefile.to_file(output_shp, driver='ESRI Shapefile')
    return output_shp

def remove_based_slope(in_shp, output_shp,slope_files, max_slope,process_num):
    # not too many polygons, not sure if slope info has been calculate, so just let it calculate agian.
    # if os.path.isfile(output_shp):
    #     basic.outputlogMessage('%s exists, skip'%output_shp)
    #     return output_shp

    # calcuate slope info
    raster_statistic.zonal_stats_multiRasters(in_shp, slope_files, tile_min_overlap=tile_min_overlap,
                                              stats=['mean', 'std', 'count'], prefix='slope',process_num=process_num)

    # remove sloep greater than max_slope
    bsmaller = False
    if vector_gpd.remove_polygons(in_shp,'slope_mean',max_slope,bsmaller,output_shp) is False:
        return False
    return output_shp

def get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-0.5, min_area=40, max_area=100000000, process_num=1,
                                 b_rm_files=False):

    save_shp = io_function.get_name_by_adding_tail(in_shp, 'post')
    if os.path.isfile(save_shp):
        basic.outputlogMessage('%s exists, skip'%save_shp)
        return save_shp


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

    # merge polygons touch each others
    wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)
    merged_shp = io_function.get_name_by_adding_tail(in_shp, 'merged')
    if filter_merge_polygons(in_shp,merged_shp,wkt, min_area,max_area,dem_diff_tif,dem_diff_thread_m,process_num) is None:
        return None

    # in merge_polygons, it will remove some big polygons, convert MultiPolygon to Polygons, so neeed to update remain_polyons
    remain_polyons = vector_gpd.read_polygons_gpd(merged_shp)

    # check MultiPolygons again.
    polyons_noMulti = [vector_gpd.MultiPolygon_to_polygons(idx, poly) for idx, poly in enumerate(remain_polyons)]
    remain_polyons = []
    for polys in polyons_noMulti:
        polys = [poly for poly in polys if poly.area > min_area]  # remove tiny polygon before buffer
        remain_polyons.extend(polys)
    print('convert MultiPolygon to polygons and remove small ones, remain %d' % (len(remain_polyons)))

    # based on the merged polygons, surrounding polygons
    buffer_surrounding = 20  # meters
    surrounding_shp = io_function.get_name_by_adding_tail(in_shp, 'surrounding')
    get_surrounding_polygons(remain_polyons, surrounding_shp, wkt, dem_diff_tif, buffer_surrounding, process_num)

    rm_reldemD_shp = io_function.get_name_by_adding_tail(in_shp, 'rmreldemD')
    if remove_polygons_based_relative_dem_diff(remain_polyons, merged_shp, surrounding_shp, wkt, rm_reldemD_shp, min_area,dem_diff_thread_m) is None:
        return None

    rm_shapeinfo_shp = io_function.get_name_by_adding_tail(in_shp, 'rmshapeinfo')
    area_limit = 10000
    circularit_limit = 0.1
    holes_count = 20
    remove_polygons_based_shapeinfo(rm_reldemD_shp, rm_shapeinfo_shp, area_limit, circularit_limit, holes_count)

    # remove based on slope
    # use the slope derived from ArcitcDEM mosaic
    slope_tif_list = io_function.get_file_list_by_ext('.tif',dem_common.arcticDEM_tile_slope_dir,bsub_folder=False)
    basic.outputlogMessage('Find %d slope files in %s'%(len(slope_tif_list), dem_common.arcticDEM_tile_slope_dir))
    rm_slope_shp = io_function.get_name_by_adding_tail(in_shp, 'rmslope')
    max_slope = 20
    if remove_based_slope(rm_shapeinfo_shp, rm_slope_shp,slope_tif_list, max_slope,process_num) is False:
        return None

    # copy
    io_function.copy_shape_file(rm_slope_shp,save_shp)

    # add date difference if they are available
    date_diff_base = os.path.basename(dem_diff_tif).replace('DEM_diff','date_diff')
    date_diff_tif = os.path.join(os.path.dirname(dem_diff_tif) , date_diff_base)
    if os.path.isfile(date_diff_tif):
        raster_statistic.zonal_stats_multiRasters(save_shp, date_diff_tif,tile_min_overlap=tile_min_overlap,
                                                  stats=['mean', 'std'], prefix='dateD',process_num=process_num)

    return save_shp



def segment_subsidence_grey_image(dem_diff_grey_8bit, dem_diff, save_dir,process_num, subsidence_thr_m=-0.5, min_area=40, max_area=100000000,
                                  b_rm_files=False):
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

    # check if the final result exist
    final_shp_path = io_function.get_name_by_adding_tail(segment_shp_path, 'post')
    if os.path.isfile(final_shp_path):
        basic.outputlogMessage('Warning, Final results (%s) of subsidence shapefile exists, skip'%final_shp_path)
        return True


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
        raster_statistic.zonal_stats_multiRasters(segment_shp_path, dem_diff, tile_min_overlap=tile_min_overlap,
                                                  stats=['mean', 'std','count'], prefix='demD',process_num=process_num)

    # get DEM diff information for each polygon.
    dem_diff_shp = get_dem_subscidence_polygons(segment_shp_path, dem_diff, dem_diff_thread_m=subsidence_thr_m,
                                 min_area=min_area, max_area=max_area, process_num=process_num,b_rm_files=b_rm_files)

    if dem_diff_shp is None:
        id_str = re.findall('grid\d+', os.path.basename(dem_diff))[0][4:]
        if len(id_str) > 1:
            grid_id = int(id_str)
            save_id_grid_no_subsidence(grid_id)
    else:
        basic.outputlogMessage('obtain elevation reduction polygons: %s'%dem_diff_shp)

    ## remove files, only keep the final results.
    if b_rm_files:
        io_function.delete_file_or_dir(label_path)
        IDrange_txt = os.path.splitext(label_path)[0] + '_IDrange.txt'
        io_function.delete_file_or_dir(IDrange_txt)
        io_function.delete_shape_file(segment_shp_path)

        # other intermediate files
        other_shp_names = ['merged','surrounding','rmreldemD','rmshapeinfo','rmslope']
        for name in other_shp_names:
            rm_shp = io_function.get_name_by_adding_tail(segment_shp_path, name)
            io_function.delete_shape_file(rm_shp)

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


# def test_get_dem_subscidence_polygons():
#     # in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_parallel_sub/WR_dem_diff_DEM_diff_prj_8bit_sub.shp')
#     # dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/WR_dem_diff_DEM_diff_prj.tif')
#     # get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-1, min_area=40, max_area=100000000)
#
#     # in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_result_grid8868/ala_north_slo_extent_latlon_grid_ids_DEM_diff_grid8868_8bit.shp')
#     # dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/ala_north_slo_extent_latlon_grid_ids_DEM_diff_grid8868.tif')
#     # get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-2, min_area=40, max_area=100000000, process_num=10)
#
#     in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_result_grid9274/WR_extent_grid_ids_DEM_diff_grid9274_8bit.shp')
#     # in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_result_grid9274/WR_extent_grid_ids_DEM_diff_grid9274_8bit_sub.shp')
#     dem_diff_tif = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/WR_extent_grid_ids_DEM_diff_grid9274.tif')
#     get_dem_subscidence_polygons(in_shp, dem_diff_tif, dem_diff_thread_m=-2, min_area=40, max_area=100000000, process_num=10)

# def test_merge_touched_polygons():
#     in_shp = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_result_grid9274/'
#                                 'WR_extent_grid_ids_DEM_diff_grid9274_8bit_sub_rmdemDArea_2.shp')
#
#     remain_polyons = vector_gpd.read_polygons_gpd(in_shp)
#
#     print(timeTools.get_now_time_str(), 'start building adjacent_matrix')
#     # adjacent_matrix = vector_features.build_adjacent_map_of_polygons(remain_polyons)
#
#     process_num = 1
#     adjacent_matrix = vector_gpd.build_adjacent_map_of_polygons(remain_polyons, process_num=process_num)
#     print(timeTools.get_now_time_str(), 'finish building adjacent_matrix')
#
#     if adjacent_matrix is False:
#         return False
#     merged_polygons = vector_features.merge_touched_polygons(remain_polyons, adjacent_matrix)
#     print(timeTools.get_now_time_str(), 'finish merging touched polygons, get %d ones' % (len(merged_polygons)))
#
#     # save
#     wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)
#     save_path = io_function.get_name_by_adding_tail(in_shp,'merged')
#     save_pd = pd.DataFrame({'Polygon': merged_polygons})
#     vector_gpd.save_polygons_to_files(save_pd, 'Polygon', wkt, save_path)


def test_remove_based_slope():

    dir = os.path.expanduser('~/Data/dem_processing/segment_parallel_9274')
    in_shp = os.path.join(dir,'WR_extent_grid_ids_DEM_diff_grid9274_8bit.shp')
    rm_shapeinfo_shp = os.path.join(dir,'WR_extent_grid_ids_DEM_diff_grid9274_8bit_rmshapeinfo.shp')
    process_num = 1

    # remove based on slope
    # use the slope derived from ArcitcDEM mosaic
    slope_tif_list = io_function.get_file_list_by_ext('.tif',dem_common.arcticDEM_tile_slope_dir,bsub_folder=False)
    basic.outputlogMessage('Find %d slope files in %s'%(len(slope_tif_list), dem_common.arcticDEM_tile_slope_dir))
    rm_slope_shp = io_function.get_name_by_adding_tail(in_shp, 'rmslope')
    max_slope = 20
    remove_based_slope(rm_shapeinfo_shp, rm_slope_shp,slope_tif_list, max_slope,process_num)


def main(options, args):

    dem_diff_path = args[0]
    io_function.is_file_exist(dem_diff_path)
    save_dir = options.save_dir
    process_num = options.process_num
    # segment_subsidence_on_dem_diff(dem_diff_path,save_dir)
    b_rm_tmp_files = options.b_remove_tmp_files

    dem_diff_grey_8bit = options.dem_diff_8bit
    if dem_diff_grey_8bit is None:
        dem_diff_grey_8bit  = get_dem_diff_8bit(dem_diff_path)

    if save_dir is None:
        save_dir = get_save_dir(dem_diff_path)

    ele_diff_thr = options.ele_diff_thr
    min_area = options.min_area
    max_area = options.max_area

    segment_subsidence_grey_image(dem_diff_grey_8bit, dem_diff_path, save_dir, process_num, subsidence_thr_m=ele_diff_thr,
                                  min_area=min_area, max_area=max_area,b_rm_files=b_rm_tmp_files)

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

    parser.add_option("", "--b_remove_tmp_files",
                      action="store_true", dest="b_remove_tmp_files",default=False,
                      help="if set, will remove intermediate files")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

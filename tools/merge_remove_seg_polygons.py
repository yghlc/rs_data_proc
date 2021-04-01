#!/usr/bin/env python
# Filename: merge_remove_seg_polygons 
"""
introduction: merge or remove some polygons from image segmentation.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 31 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import vector_gpd
import vector_features
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.timeTools as timeTools

import pandas as pd


def remove_merge_polygon_in_one_shp(in_shp, org_raster, attribute_name, attribute_range, min_area, max_area, process_num=1):
    # attribute_range: [min, max],

    lower = attribute_range[0]
    upper = attribute_range[1]

    save_shp = io_function.get_name_by_adding_tail(in_shp, 'post')
    if os.path.isfile(save_shp):
        basic.outputlogMessage('%s exists, skip'%save_shp)
        return save_shp

    shp_pre = io_function.get_name_no_ext(in_shp)
    # read polygons and label from segment algorithm, note: some polygons may have the same label
    polygons, attr_value_list = vector_gpd.read_polygons_attributes_list(in_shp,attribute_name)
    print('Read %d polygons'%len(polygons))
    if attr_value_list is None:
        raise ValueError('%s not in %s, need to remove it and then re-create'%(attribute_name,in_shp))


    remain_polyons = []
    rm_min_area_count = 0
    rm_att_value_count = 0
    for poly, att_value in zip(polygons, attr_value_list):
        if poly.area < min_area:
            rm_min_area_count += 1
            continue
        if lower is None:
            if att_value >= upper:
                rm_att_value_count += 1
                continue
        elif upper is None:
            if att_value <= lower:
                rm_att_value_count += 1
                continue
        else:
            # out of range, rmeove
            if att_value < lower or att_value > upper:
                rm_att_value_count += 1
                continue
        remain_polyons.append(poly)

    print('remove %d polygons based on min_area, %d polygons based on attribute_range, remain %d ones'%(rm_min_area_count, rm_diff_thr_count,len(remain_polyons)))

    if len(remain_polyons) > 1:
        # we should only merge polygon with similar reduction, but we already remove polygons with mean reduction > threshhold
        # merge touch polygons
        print(timeTools.get_now_time_str(), 'start building adjacent_matrix')
        # adjacent_matrix = vector_features.build_adjacent_map_of_polygons(remain_polyons)
        machine_name = os.uname()[1]
        # if 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:
        #     print('Warning, some problem of parallel running in build_adjacent_map_of_polygons on curc, but ok in my laptop and uist, change process_num = 1')
        #     process_num = 1
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
    print('convert MultiPolygon to polygons, remove some small polygons, remain %d' % (len(remain_polyons)))


    # based on the merged polygons, calculate the mean dem diff, relative dem_diff
    buffer_surrounding = 20  # meters
    surrounding_polygons = vector_gpd.get_surrounding_polygons(remain_polyons,buffer_surrounding)
    surrounding_shp = io_function.get_name_by_adding_tail(in_shp, 'surrounding')
    surr_pd = pd.DataFrame({'Polygon': surrounding_polygons})
    vector_gpd.save_polygons_to_files(surr_pd, 'Polygon', wkt, surrounding_shp)
    raster_statistic.zonal_stats_multiRasters(surrounding_shp, org_raster, stats=['mean', 'std', 'count'], prefix='demD',process_num=process_num)


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



def main(options, args):

    b_merge_labels = options.b_merge_labels
    b_merge_touch_polys = options.b_merge_touch_polygons
    process_num = options.process_num
    org_raster = options.org_raster

    shp_path_list = [item for item in args if io_function.is_file_exist(item)]


    # get DEM diff information for each polygon.
    # dem_diff_shp = get_dem_subscidence_polygons(segment_shp_path, dem_diff, dem_diff_thread_m=subsidence_thr_m,
    #                              min_area=min_area, max_area=max_area, process_num=process_num)


    pass

if __name__ == '__main__':
    usage = "usage: %prog [options] shp_path1 shp_path2 ...  "
    parser = OptionParser(usage=usage, version="1.0 2021-3-31")
    parser.description = 'Introduction: merge and remove polygons from the image segmentation.'

    parser.add_option("-m", "--b_merge_labels",
                      action="store_true", dest="b_merge_labels", default=False,
                      help="for multiple labels, merge them after polygonize")

    parser.add_option("-t", "--b_merge_touch_polygons",
                      action="store_true", dest="b_merge_touch_polygons", default=False,
                      help="merge polygons if they touch each other and have similar attributes")

    parser.add_option("-o", "--org_raster",
                      action="store", dest="org_raster",
                      help="the original raster for calculating the attributes")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

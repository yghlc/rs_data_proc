#!/usr/bin/env python
# Filename: rm_large_small_polygons 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 April, 2021
"""

import os,sys
from optparse import OptionParser
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd

import geopandas as gpd

import basic_src.io_function as io_function
import basic_src.basic as basic

def list_to_dict(list_dict):
    out_dict = {}
    for dict_obj in list_dict:
        for key in dict_obj.keys():
            if key in out_dict.keys():
                out_dict[key].append(dict_obj[key])
            else:
                out_dict[key] = [dict_obj[key]]
    return out_dict

def refine_dem_reduction_polygons(input_shp,min_area,max_area):

    rm_min_max_area = os.path.basename(io_function.get_name_by_adding_tail(input_shp, 'MinMaxArea'))
    vector_gpd.remove_polygons_not_in_range(input_shp,'poly_area', min_area, max_area,rm_min_max_area)

    # rm_min_max_area = input_shp

    # add some shape info
    shapefile = gpd.read_file(rm_min_max_area)
    polygons = vector_gpd.read_polygons_gpd(rm_min_max_area)
    shape_info_list = [ vector_gpd.calculate_polygon_shape_info(item)  for item in polygons]

    shapeinfo_all_dict = list_to_dict(shape_info_list)
    vector_gpd.add_attributes_to_shp(rm_min_max_area,shapeinfo_all_dict)

    # INarea, INperimete, circularit, hole_count,   WIDTH & HEIGHT (minimum_rotated_rectangle), ratio_w_h

    # remove based on shape info
    output = os.path.basename(io_function.get_name_by_adding_tail(input_shp, 'final'))
    remove_count = 0
    for idx, row in shapefile.iterrows():
        shape_info = shape_info_list[idx]

        # remove quite large but narrow ones
        if shape_info['INarea'] > 10000 and shape_info['circularit'] < 0.3 :
            shapefile.drop(idx, inplace=True)
            remove_count += 1

    basic.outputlogMessage('remove %d polygons based on shapeinfo, remain %d ones saving to %s' %
                           (remove_count, len(shapefile.geometry.values), output))
    # save results
    shapefile.to_file(output, driver='ESRI Shapefile')

def main(options, args):

    input_shp = args[0]
    io_function.is_file_exist(input_shp)

    min_area = options.min_area
    max_area = options.max_area
    refine_dem_reduction_polygons(input_shp,min_area,max_area)


def test_main():
    shp=os.path.expanduser('~/ala_north_slo_extent_latlon_DEM_diff_grid11473_8bit_post_MinMaxArea.shp')
    refine_dem_reduction_polygons(shp, 120, 1000*1000)



if __name__ == "__main__":
    usage = "usage: %prog [options] input_shp "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment subsidence based on DEM difference  '

    parser.add_option("-l", "--min_area",
                      action="store", dest="min_area", type=float, default=40,
                      help="the minimum area for each polygon")

    parser.add_option("-u", "--max_area",
                      action="store", dest="max_area", type=float, default=100000000,  # 10 km by 10 km
                      help="the maximum area for each polygon")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

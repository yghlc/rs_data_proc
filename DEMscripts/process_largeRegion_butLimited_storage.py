#!/usr/bin/env python
# Filename: process_largeRegion_butLimited_storage.py
"""
introduction: to process large areas of ArcticDEM but with limited storage,
we will process grid by grid.

We will first download ArcticDEM strip for needed grids, then remove the strip file if the grids it covered have been processed

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 11 October, 2021
"""

import os,sys
from optparse import OptionParser
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import vector_gpd
import slurm_utility
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.timeTools as timeTools

import re
from datetime import datetime
import time

import multiprocessing
from multiprocessing import Pool


machine_name = os.uname()[1]

# some files
from dem_common import grid_20_shp,grid_20_id_raster,dem_strip_shp,dem_tile_shp

# log or txt files
from dem_common import process_log_dir, grid_complete_list, strip_dem_cover_grids, tile_dem_cover_grids
if os.path.isdir(process_log_dir) is False:
    io_function.mkdir(process_log_dir)

from produce_DEM_diff_ArcticDEM import get_grid_20


def download_arctic_dem():
    # download tarball
    pass

def get_overlap_grids_for_one_extent(all_ids,all_grid_polys, dem_poly, dem_name, idx, dem_poly_count):
    print(timeTools.get_now_time_str(), idx, dem_poly_count)
    index = vector_gpd.get_poly_index_within_extent(all_grid_polys, dem_poly)
    gird_ids = [all_ids[idx] for idx in index]
    return dem_name,gird_ids

def build_dict_of_dem_cover_grid_ids(dem_info_shp,grid_20_shp,save_dict_txt):
    # this will take time, but only need to run once at the beginning
    if os.path.isfile(save_dict_txt):
        print('warning, %s exists, skip build_dict_of_dem_cover_grid_ids'%save_dict_txt)
        return True

    # extent polygons and projection (proj4)
    dem_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(dem_info_shp)
    if dem_shp_prj == '':
        raise ValueError('get proj4 of %s failed' % dem_info_shp)
    grid_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)
    if grid_shp_prj == '':
        raise ValueError('get proj4 of %s failed' % grid_20_shp)

    if dem_shp_prj != grid_shp_prj:
        raise ValueError('%s and %s do not have the same projection'% (dem_info_shp, grid_20_shp))

    # read DEM info
    dem_polygons, dem_names = vector_gpd.read_polygons_attributes_list(dem_info_shp, 'name',b_fix_invalid_polygon=False)
    # dem_name: eg. SETSM_GE01_20090818_1050410001E0CF00_1050410001D80200_seg1_2m_v3.0  or 11_27_2_1_2m_v3.0
    dem_poly_count = len(dem_polygons)
    # check if there is duplicate dem names
    if len(dem_names) != len(set(dem_names)):
        raise ValueError('some duplicate dem name in %s'%dem_info_shp)

    # read grid polygons and ids
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')

    dem_cover_grids = {}
    # this will take time.
    # for idx, (dem_poly,dem_name) in enumerate(zip(dem_polygons, dem_names)):
    #     print(timeTools.get_now_time_str(), idx, dem_poly_count)
    #     index = vector_gpd.get_poly_index_within_extent(all_grid_polys, dem_poly)
    #     gird_ids = [ all_ids[idx] for idx in index ]
    #     # if dem_name in dem_cover_grids.keys():
    #     #     basic.outputlogMessage('\n Warning, %s already in dict \n'%dem_name)
    #     dem_cover_grids[dem_name] = gird_ids

    ### parallel version
    theadPool = Pool(multiprocessing.cpu_count())  # multi processes
    parameters_list = [(all_ids,all_grid_polys, dem_poly,dem_name, idx, dem_poly_count) for idx, (dem_poly,dem_name) in enumerate(zip(dem_polygons, dem_names))]
    results = theadPool.starmap(get_overlap_grids_for_one_extent, parameters_list)  # need python3
    for res in results:
        dem_name, gird_ids = res
        dem_cover_grids[dem_name] = gird_ids

    # save to dict
    io_function.save_dict_to_txt_json(save_dict_txt,dem_cover_grids)
    return True

def test_build_dict_of_dem_cover_grid_ids():
    # for test, only have three polygons
    dem_info_shp = os.path.expanduser('~/Data/tmp_data/test_ArcticDEM/ArcticDEM_Tile_Index_Rel7.shp')
    save_dict_txt = tile_dem_cover_grids
    build_dict_of_dem_cover_grid_ids(dem_info_shp, grid_20_shp, save_dict_txt)

def main(options, args):
    extent_shp = args[0]

    # build map dem cover grid (take time, but only need to run once at the beginning)
    build_dict_of_dem_cover_grid_ids(dem_strip_shp, grid_20_shp, strip_dem_cover_grids)
    build_dict_of_dem_cover_grid_ids(dem_tile_shp, grid_20_shp, tile_dem_cover_grids)

    # get grids for processing

    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    grid_polys, grid_ids = get_grid_20(extent_shp,all_grid_polys, all_ids)


    if machine_name == 'ubuntu' or machine_name == 'uist':
        # download and unpack ArcticDEM, do registration, send to curc
        download_arctic_dem()
    elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:  # curc
        # process ArcticDEM using the computing resource on CURC
        pass
    else:
        print('unknown machine name')




if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp "
    parser = OptionParser(usage=usage, version="1.0 2021-10-11")
    parser.description = 'Introduction: process ArcticDEM with limited storage  '

    parser.add_option("-n", "--max_grids",
                      action="store", dest="max_grids", type=int, default=100,
                      help="the maximum number of grids for process in parallel, large storage allows large number")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

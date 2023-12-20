#!/usr/bin/env python
# Filename: select_grids_data_overlap_vectors.py 
"""
introduction: select grid data (DEM difference raster, polygons, and others) if they overlap input polygons or points

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 December, 2023
"""

import os, sys
from optparse import OptionParser
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import pandas as pd
import geopandas as gpd
import vector_gpd

arcticDEM_res_dir ="/Volumes/Seagate8T/ArcticDEM_results"

def find_grids_overlap_vector_shp(grid_indexes_shp,vector_shp_list):

    if isinstance(vector_shp_list,list) is False:
        vector_shp_list = [vector_shp_list]

    overlap_touch_list = []
    for vector_shp in vector_shp_list:
        overlap_touch = vector_gpd.geometries_overlap_another_group(grid_indexes_shp,vector_shp)
        overlap_touch_list.append(overlap_touch)

    overlap_touch_all = pd.concat(overlap_touch_list)

    # remove duplicated geometries in overlap_touch
    overlap_touch_all = overlap_touch_all.drop_duplicates(subset=['geometry'])    # only check geometry
    overlap_touch_all.to_file('overlap_touch_all.shp')

    return overlap_touch_all

def save_grids_ids_to_txt(sel_grids_gpd, save_txt_path):
    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    grid_id_list = [str(item) for item in grid_id_list]
    io_function.save_list_to_txt(save_txt_path,grid_id_list)

def save_grids_fileurl_to_txt(sel_grids_gpd, save_txt_path):
    grid_id_list = sel_grids_gpd['fileurl'].to_list()
    io_function.save_list_to_txt(save_txt_path, grid_id_list)

    ext_num_list = [ item.split('/')[1] for item in grid_id_list]
    ext_num_list = list(set(ext_num_list))
    io_function.save_list_to_txt('ext_num_list.txt', ext_num_list)



def find_dem_difference_raster(sel_grids_gpd, save_txt_path):
    grid_id_list = sel_grids_gpd['cell_id'].to_list()
    fileurl_list = sel_grids_gpd['fileurl'].to_list()
    ext_num_list = [item.split('/')[1] for item in fileurl_list]

    dem_difference_file_list = []
    for grid_id, ext_num in zip(grid_id_list,ext_num_list):
        pattern = '%s*/grid_dem_diffs/*grid%d*.*'%(ext_num,grid_id)
        file_list = io_function.get_file_list_by_pattern(arcticDEM_res_dir,pattern)
        pattern2 = '%s*/grid_dem_diffs/*grid%d*_*' % (ext_num, grid_id)
        file_list2 = io_function.get_file_list_by_pattern(arcticDEM_res_dir, pattern2)
        file_list.extend(file_list2)
        if len(file_list) > 0:
            dem_difference_file_list.extend(file_list)
        else:
            print('warning: there is no files for grid %d in ext: %s'%(grid_id, ext_num))
            print(arcticDEM_res_dir, pattern)

    io_function.save_list_to_txt(save_txt_path, dem_difference_file_list)


def test_find_grids_overlap_vector_shp():
    data_dir = os.path.expanduser('~/Data/dem_processing/products_derived_from_ArcticDEM/Index_shp')
    grid_indexes_shp = os.path.join(data_dir,'elevation-differences_Index/elevation-differences_Index.shp')
    # vector_shp = os.path.expanduser('~/codes/PycharmProjects/rts-site-data-review/rts_research_sites_point.shp')

    vector_shp = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_point.shp')
    # vector_shp2 = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_edit.shp')
    vector_shp2 = os.path.expanduser('~/Data/slump_demdiff_classify/rts_research_sites/rts_research_sites_arctic_small.shp')
    find_grids_overlap_vector_shp(grid_indexes_shp,[vector_shp,vector_shp2])
    # find_grids_overlap_vector_shp(grid_indexes_shp,vector_shp2)


def main(options, args):
    grid_indexes_shp= args[0]
    vector_shp_list = args[1:]
    if len(vector_shp_list) < 1:
        raise ValueError("No input vector files")

    #   cell_id (Integer64) = 58601
    #   tarball (String) = dem_diffs_2m_grid58601.tar.gz
    #   fileurl (String) = elevation-differences/ext11/dem_diffs_2m_grid58601.tar.gz
    #   POLYGON ((3320000 -120000,3340000 -120000,3340000 -140000,3320000 -140000,3320000 -120000))
    # read the grid index file
    # grid_polys, grid_attributes = vector_gpd.read_polygons_attributes_list(grid_indexes_shp,['cell_id','fileurl'],
    #                                                                        b_fix_invalid_polygon=False)

    sel_grids_gpd = find_grids_overlap_vector_shp(grid_indexes_shp, vector_shp_list)
    save_grids_ids_to_txt(sel_grids_gpd,'select_grids_ids.txt')
    save_grids_fileurl_to_txt(sel_grids_gpd,'select_grids_fileurls.txt')
    find_dem_difference_raster(sel_grids_gpd, 'select_dem_diff_raster_list.txt')



if __name__ == '__main__':
    usage = "usage: %prog [options] grid_indexes_shp vector_shp1 vector_shp2 ... "
    parser = OptionParser(usage=usage, version="1.0 2023-12-18")
    parser.description = 'Introduction: select grid data (DEM diff, polygons etc) if they overlap input polygons or points '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='select_grid_data',
                      help="the folder to save DEMs")

    # test_find_grids_overlap_vector_shp()
    # sys.exit(0)

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

#!/usr/bin/env python
# Filename: process_panArctic.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 19 January, 2022
"""


import os,sys
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic
import vector_gpd

from dem_common import grid_ids_txt_dir,get_extent_grid_id_txt_done_files
from dem_common import grid_20_shp
from produce_DEM_diff_ArcticDEM import get_grid_20

ext_shp_dir = os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/ArcticDEM_subsets')
if os.path.isdir(grid_ids_txt_dir) is False:
    io_function.mkdir(grid_ids_txt_dir)

def read_grid_ids_from_other_extent():
    grid_txt_list = io_function.get_file_list_by_ext('.txt', grid_ids_txt_dir, bsub_folder=False)
    other_grid_ids = []
    for txt in grid_txt_list:
        id_list = io_function.read_list_from_txt(txt)
        other_grid_ids.extend(id_list)

    other_grid_ids = [ int(item) for item in other_grid_ids ]
    return other_grid_ids

def produce_corresponding_grid_ids_txt(extent_shp,local_grid_id_txt, log_grid_ids_txt):

    # if it in the logdir, not the current dir, then copy it
    if os.path.isfile(log_grid_ids_txt) and os.path.isfile(local_grid_id_txt) is False:
        io_function.copy_file_to_dst(log_grid_ids_txt,local_grid_id_txt,overwrite=False)
        return True

    # if not in the local dir, then generate it
    if os.path.isfile(local_grid_id_txt) is False:
        # read grids and ids
        time0 = time.time()
        all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
        print('time cost of read polygons and attributes', time.time() - time0)

        # this will create local_grid_id_txt
        grid_polys, grid_ids = get_grid_20(extent_shp, all_grid_polys, all_ids)

        # modify local_grid_id_txt by excluding grid_id already in adjacent extent
        other_grid_ids = read_grid_ids_from_other_extent()
        grid_ids = [id for id in grid_ids if id not in other_grid_ids ]

        # over write local_grid_id_txt file
        grid_ids_str = [str(item) for item in grid_ids]
        io_function.copy_file_to_dst(local_grid_id_txt,io_function.get_name_by_adding_tail(local_grid_id_txt, 'noRMadj')) # save a copy
        io_function.save_list_to_txt(local_grid_id_txt, grid_ids_str)

        # copy to log dir
        io_function.copy_file_to_dst(local_grid_id_txt,log_grid_ids_txt)

    return True


def process_one_extent(extent_shp):

    # local_grid_id_txt is in the current dir
    # log_grid_ids_txt, log_grid_ids_txt_done is in grid_ids_txt_dir
    local_grid_id_txt, log_grid_ids_txt, log_grid_ids_txt_done = get_extent_grid_id_txt_done_files(extent_shp)

    # check if it has been complete
    if os.path.isfile(log_grid_ids_txt_done):
        basic.outputlogMessage('Tasks for extent %s have been completed'%extent_shp)
        return True

    produce_corresponding_grid_ids_txt(extent_shp,local_grid_id_txt, log_grid_ids_txt)

    res = os.system('./process_largeRegion.sh %s')
    if res != 0:
        sys.exit(res)


def main():

    ext_shps = io_function.get_file_list_by_ext('.shp',ext_shp_dir,bsub_folder=False)
    basic.outputlogMessage('%d extent shapefiles in %s'%(len(ext_shps), ext_shp_dir))
    for shp in ext_shps:
        process_one_extent(shp)


if __name__ == '__main__':
    main()
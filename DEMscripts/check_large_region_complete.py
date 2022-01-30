#!/usr/bin/env python
# Filename: check_large_region_complete.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 28 January, 2022
"""
import os,sys
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic

from datetime import datetime

from process_panArctic import get_ext_shps, ext_shp_dir
from dem_common import get_extent_grid_id_txt_done_files
from process_largeRegion_butLimited_storage import get_complete_ignore_grid_ids,save_grid_ids_need_to_process, make_note_all_task_done

from process_largeRegion_butLimited_storage import update_complete_grid_list

def check_one_extent(extent_shp):
    print('start to check %s' % extent_shp)

    # local_grid_id_txt is in the current dir
    # log_grid_ids_txt, log_grid_ids_txt_done is in grid_ids_txt_dir
    local_grid_id_txt, log_grid_ids_txt, log_grid_ids_txt_done = get_extent_grid_id_txt_done_files(extent_shp)
    if os.path.isfile(local_grid_id_txt) is False:
        print('the _grid_ids.txt for %s does not exist, maybe it has started'%extent_shp)
        return False

    # check if it has been complete
    if os.path.isfile(log_grid_ids_txt_done):
        basic.outputlogMessage('Tasks for extent %s have been completed' % extent_shp)
        return True

    grid_ids_to_process_txt = io_function.get_name_no_ext(extent_shp) + '_' + 'grid_ids_to_process.txt'

    # read from txt file directly
    grid_ids = [int(item) for item in io_function.read_list_from_txt(local_grid_id_txt)]

    update_complete_grid_list(grid_ids, task_list)

    # check complete files, to see if it's done
    # remove grids that has been complete or ignored
    ignore_ids = get_complete_ignore_grid_ids()
    num_grid_ids = save_grid_ids_need_to_process(grid_ids, ignore_ids=ignore_ids, save_path=grid_ids_to_process_txt)
    if num_grid_ids < 1:
        print(datetime.now(),' %s is marked as completed'%extent_shp)
        make_note_all_task_done(extent_shp)
    else:
        print(datetime.now(), ' %s has not completed, %d grids to process, total: %d' % (extent_shp, num_grid_ids, len(grid_ids)))

    return True

def main():
    ext_shps = get_ext_shps()

    basic.outputlogMessage('%d extent shapefiles in %s'%(len(ext_shps), ext_shp_dir))
    for shp in ext_shps:
        check_one_extent(shp)

if __name__ == '__main__':
    task_list = ['dem_diff','dem_headwall_grid', 'hillshade_headwall_line','segment']
    main()

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
from process_largeRegion_butLimited_storage import get_subset_info_txt_list
from process_largeRegion_butLimited_storage import get_subset_info
from process_largeRegion_butLimited_storage import subset_shp_dir
from process_largeRegion_butLimited_storage import msg_file_pre
from dem_common import subset_message_dir


from dem_segment_subsidence_jobs import curc_node

process_node = '$curc_host' #if options.process_node is None else options.process_node
r_working_dir = '/scratch/summit/lihu9680/Arctic/dem_processing'  #if options.remote_working_dir is None else options.remote_working_dir

def find_fail_grids_in_complete_subsets(extent_shp, grid_ids_to_process_txt):


    msg_file_pre_name = io_function.get_name_no_ext(extent_shp) + '_' + msg_file_pre
    print('msg_file_pre_name:', msg_file_pre_name)

    # subset shp dir
    subset_shp_txt_dir = subset_shp_dir + '_' + io_function.get_name_no_ext(extent_shp)
    subset_shp_txt_dir = os.path.join(subset_message_dir, subset_shp_txt_dir)
    print('subset_shp_txt_dir:',subset_shp_txt_dir)

    # subset_txt_list = get_subset_info_txt_list('proc_status', ['done'], local_folder='./')
    remote_sub_done_list = get_subset_info_txt_list('proc_status', ['done'], remote_node=process_node,
                                              remote_folder=r_working_dir, msg_pre=msg_file_pre_name)
    to_process_grids = io_function.read_list_from_txt(grid_ids_to_process_txt)
    all_fail_grids = []
    for done_subset in remote_sub_done_list:

        fail_grids = []

        subset_info_txt_path = os.path.join(subset_message_dir,done_subset)
        print('subset_info_txt_path:',subset_info_txt_path)
        subset_info = get_subset_info(subset_info_txt_path,dir='./')
        sub_ext_shp = subset_info['shp']
        print('sub_ext_shp:',sub_ext_shp)
        subset_grid_txt = os.path.join(subset_shp_txt_dir, os.path.splitext(sub_ext_shp)[0] + '_grid_ids.txt')
        print('subset_grid_txt:',subset_grid_txt)
        subset_grid_ids = io_function.read_list_from_txt(subset_grid_txt)
        for id in subset_grid_ids:
            if id in to_process_grids:
                print('fail grid %s of %s '%(id,sub_ext_shp))
                fail_grids.append(id)

        if len(fail_grids) > 0:
            subset_fail_grid_txt = io_function.get_name_by_adding_tail(os.path.basename(subset_grid_txt),'fail')
            io_function.save_list_to_txt(subset_fail_grid_txt,fail_grids)

        all_fail_grids.extend(fail_grids)

    if len(all_fail_grids) > 0:
        fail_grid_txt = io_function.get_name_by_adding_tail(os.path.basename(grid_ids_to_process_txt), 'fail')
        io_function.save_list_to_txt(fail_grid_txt, all_fail_grids)

    return True

def check_one_extent(extent_shp):
    print('start to check %s' % extent_shp)

    # local_grid_id_txt is in the current dir
    # log_grid_ids_txt, log_grid_ids_txt_done is in grid_ids_txt_dir
    local_grid_id_txt, log_grid_ids_txt, log_grid_ids_txt_done = get_extent_grid_id_txt_done_files(extent_shp)
    if os.path.isfile(local_grid_id_txt) is False and os.path.isfile(log_grid_ids_txt):
        io_function.copy_file_to_dst(log_grid_ids_txt,local_grid_id_txt)
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
        make_note_all_task_done(extent_shp,curc_node)
    else:
        print(datetime.now(), ' %s has not completed, %d grids to process, total: %d' % (extent_shp, num_grid_ids, len(grid_ids)))
        find_fail_grids_in_complete_subsets(extent_shp,grid_ids_to_process_txt)

    return True

def main():
    ext_shps = get_ext_shps()

    basic.outputlogMessage('%d extent shapefiles in %s'%(len(ext_shps), ext_shp_dir))
    for shp in ext_shps:
        check_one_extent(shp)

if __name__ == '__main__':
    task_list = ['dem_diff','dem_headwall_grid', 'hillshade_headwall_line','segment']
    main()

#!/usr/bin/env python
# Filename: dem_segment_subsidence_jobs.py 
"""
introduction: because CURC are too busy, I decide to run segment jobs in local workstation (tesia, uist, and donostia)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 January, 2022
"""

import os,sys
import re
import time
from datetime import datetime

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import scp_communicate

import dem_common
from dem_common import get_grid_id_from_path

from dem_common import grid_dem_diffs_segment_dir
curc_node = '$curc_host'
r_seg_res_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir/grid_dem_diffs_segment_results'

machine_name = os.uname()[1]

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, code_dir)
from tools.move_old_files_folders import check_file_or_dir_is_old

def get_dem_diff_list_to_seg():

    dem_diff_list = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_dir, '*DEM_diff_grid*.tif')
    dem_diff_list_copy = dem_diff_list.copy()
    if len(dem_diff_list) < 1:
        print('No *DEM_diff*.tif in %s' % dem_common.grid_dem_diffs_dir)
        return []

    dem_diff_ids = [int(re.findall('grid\d+', os.path.basename(item))[0][4:]) for item in dem_diff_list]

    dem_subsidence_shps = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_segment_dir,
                                                               'segment_result*/*_post.shp')
    subsidence_ids = [int(re.findall('grid\d+', os.path.basename(item))[0][4:]) for item in dem_subsidence_shps]
    no_subsidence_ids = []
    if os.path.isfile(dem_common.grid_no_subscidence_poly_txt):
        no_subsidence_ids = [int(item) for item in
                             io_function.read_list_from_txt(dem_common.grid_no_subscidence_poly_txt)]
    subsidence_ids.extend(no_subsidence_ids)

    # remove dem diff already been segmented or no subsidence
    for id, dem_diff in zip(dem_diff_ids, dem_diff_list_copy):
        if id in subsidence_ids:
            dem_diff_list.remove(dem_diff)

    return dem_diff_list

def get_dem_diff_old_enough(dem_diff_list):
    dem_diff_old = []
    # if a dem diff is 24 hours old,
    for demdiff in dem_diff_list:
        if check_file_or_dir_is_old(demdiff,24):
            dem_diff_old.append(demdiff)

    return dem_diff_old

def read_dem_diff_assigned_to_other_machine(job_list_pre):

    dem_list_txts = io_function.get_file_list_by_pattern(dem_common.process_log_dir, job_list_pre + '*.txt')
    assign_dem_diff = []
    for txt in dem_list_txts:
        assign_dem_diff.extend(io_function.read_list_from_txt(txt))

    return assign_dem_diff

def copy_segment_result_to_curc(grid_ids):
    '''after complete, copy the results to curc, avoid submit a new jobs again'''
    for idx, g_id in enumerate(grid_ids):
        print('(%d / %d) Transfer DEM diff segment results to CURC '%(idx + 1,len(grid_ids)))
        dem_diff_seg_res_folder = os.path.join(grid_dem_diffs_segment_dir,'segment_result_grid%d'%g_id)
        scp_communicate.copy_file_folder_to_remote_machine(curc_node,r_seg_res_dir,dem_diff_seg_res_folder)

def produce_products_dem_subsidence(b_remove_job_folder=True):
    # run segment jobs in local workstations.

    task = 'segment'
    max_list_count = 20
    if 'donostia' in machine_name:
        max_list_count = 8      # donostia is really slow, assigined less task to it
    job_list_pre = 'job_seg_dem_diff_list_'

    if os.path.isdir(dem_common.process_log_dir) is False:
        io_function.mkdir(dem_common.process_log_dir)

    dem_list_txt = os.path.join(dem_common.process_log_dir, job_list_pre + machine_name + '.txt')

    # when submit segment of dem_diff, no need ext_shp
    ext_shp = "monitor_fail_segment_jobs"

    while True:
        dem_diff_list = get_dem_diff_list_to_seg()

        # only handle file are old enough
        dem_diff_list = get_dem_diff_old_enough(dem_diff_list)

        dem_diff_ids = [get_grid_id_from_path(item) for item in dem_diff_list]
        print('dem_diff_ids')
        print(dem_diff_ids)

        # remove dem_diff already assigined for other machine
        if os.path.isfile(dem_list_txt):
            io_function.delete_file_or_dir(dem_list_txt)
        dem_diff_assigned = read_dem_diff_assigned_to_other_machine(job_list_pre)
        assigned_ids = [get_grid_id_from_path(item) for item in dem_diff_assigned]
        print('assigned_ids')
        print(assigned_ids)
        keep_idx = [idx for idx,id in enumerate(dem_diff_ids) if id not in assigned_ids ]
        dem_diff_list = [ dem_diff_list[item] for item in keep_idx]

        if len(dem_diff_list) < 1:
            print(datetime.now(), 'there is no DEM_diff for %s to seg, wait 10 minutes'%machine_name)
            time.sleep(600)     # wait 10 min
            continue

        # save some of them to txt, for "parallel_processing_curc.py"
        dem_diff_list = dem_diff_list[:max_list_count]
        save_ids = [get_grid_id_from_path(item) for item in dem_diff_list]
        print('save_ids')
        print(save_ids)

        io_function.save_list_to_txt(dem_list_txt,dem_diff_list)

        res = os.system('./run.sh %s %s' % (ext_shp, task))
        if res != 0:
            sys.exit(1)

        copy_segment_result_to_curc(save_ids)

        if b_remove_job_folder:
            os.system('rm -r seg_dem_diff_*')
            io_function.delete_file_or_dir(dem_list_txt)


def test_copy_segment_result_to_curc():
    ids = [5095,5275]
    copy_segment_result_to_curc(ids)

def main():
    tasks =  ['segment']
    produce_products_dem_subsidence()
    pass

if __name__ == '__main__':
    # test_copy_segment_result_to_curc()
    main()
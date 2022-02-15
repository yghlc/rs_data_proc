#!/usr/bin/env python
# Filename: monitor_fail_jobs.py
"""
introduction: monitor and processing failed jobs.

# run on CURC compile node.

Currently, all task includes: dem_diff dem_headwall_grid  hillshade_headwall_line segment

        segment is handle by "dem_segment_subsidence_jobs.py" on UIST.
        This is script will run on CURC compile node, handling dem_diff and dem_headwall_grid

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 15 February, 2022
"""


import os,sys
from datetime import datetime
import time

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, code_dir)
from dem_common import process_log_dir
from tools.move_old_files_folders import check_file_or_dir_is_old

machine_name = os.uname()[1]

def get_failed_grid_ids(task):
    if task== 'dem_diff':
        fail_log_dir = os.path.join(process_log_dir,'get_dem_diff')
    elif task == 'dem_headwall_grid':
        fail_log_dir = os.path.join(process_log_dir, 'extract_headwall_from_slope_grid')
    else:
        raise ValueError('Unknow task: %s'%str(task))

    fail_id_txt_list = io_function.get_file_list_by_ext('.txt',fail_log_dir,bsub_folder=False)
    old_id_txts = []

    # check they are old enough (24 hours)
    for txt in fail_id_txt_list:
        if check_file_or_dir_is_old(txt,24,print_time=True):
            old_id_txts.append(txt)
    return old_id_txts

def merge_grid_ids_txt(task, fail_id_txt_list):
    id_list = []
    for txt in fail_id_txt_list:
        id_list.extend(io_function.read_list_from_txt(txt))
    id_list = list(set(id_list)) # remove redudant ones
    save_path = '%s_fail_grid_ids.txt'%task
    io_function.save_list_to_txt(save_path,id_list)
    return save_path

def monitor_process_failed_grids():
    tasks = ['dem_diff','dem_headwall_grid']

    while True:
        all_fail_list = []
        for task in tasks:
            all_fail_list.extend(get_failed_grid_ids(task))
        if len(all_fail_list) < 1:
            print(datetime.now(), 'there is no failed jobs for %s to process, wait 10 minutes' % machine_name)
            time.sleep(600)
            continue

        for task in tasks:
            fail_id_txt_list = get_failed_grid_ids(task)
            if len(fail_id_txt_list) < 1:
                continue
            print('\n',task,'fail_id_txt_list:',fail_id_txt_list,'\n')

            fail_ids_txt = merge_grid_ids_txt(task,fail_id_txt_list)
            # start processing
            res = os.system('./run.sh %s %s' % (fail_ids_txt, task))
            if res != 0:
                sys.exit(1)

            # wait all job done? Yes, in "parallel_processing_curc.py"

            # remove fail txt and working dir
            for txt in fail_id_txt_list:
                io_function.delete_file_or_dir(txt)
            if task == 'dem_diff':
                os.system('rm -r dem_diff_*')
            if task == 'dem_headwall_grid':
                os.system('rm -r extract_headwall_grid_*')



def main():
    monitor_process_failed_grids()

if __name__ == '__main__':
    # test_copy_segment_result_to_curc()
    main()


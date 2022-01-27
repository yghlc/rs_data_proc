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

import dem_common

machine_name = os.uname()[1]

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


def read_dem_diff_assigned_to_other_machine(job_list_pre):

    dem_list_txts = io_function.get_file_list_by_pattern(dem_common.process_log_dir, job_list_pre + '*.txt')
    assign_dem_diff = []
    for txt in dem_list_txts:
        assign_dem_diff.extend(io_function.read_list_from_txt(txt))

    return assign_dem_diff


def produce_products_dem_subsidence(b_remove_job_folder=True):
    # run segment jobs in local workstations.

    task = 'segment'
    max_list_count = 20
    job_list_pre = 'job_seg_dem_diff_list_'

    if os.path.isdir(dem_common.process_log_dir) is False:
        io_function.mkdir(dem_common.process_log_dir)

    dem_list_txt = os.path.join(dem_common.process_log_dir, job_list_pre + machine_name + '.txt')

    # when submit segment of dem_diff, no need ext_shp
    ext_shp = "None"

    while True:
        dem_diff_list = get_dem_diff_list_to_seg()
        dem_diff_name_list = [os.path.basename(item) for item in dem_diff_list]

        # remove dem_diff already assigined for other machine
        dem_diff_assigned = read_dem_diff_assigned_to_other_machine(job_list_pre)
        dem_diff_assigned = [os.path.basename(item) for item in dem_diff_assigned]
        for name, full_name in zip(dem_diff_name_list,dem_diff_list):
            if name in dem_diff_assigned:
                dem_diff_list.remove(full_name)

        if len(dem_diff_list) < 1:
            print(datetime.now(), 'there is no DEM_diff for %s to seg, wait 10 minutes'%machine_name)
            time.sleep(600)     # wait 10 min

        # save some of them to txt, for "parallel_processing_curc.py"
        dem_diff_list = dem_diff_list[:max_list_count]
        io_function.save_list_to_txt(dem_list_txt,dem_diff_list)

        res = os.system('./run.sh %s %s' % (ext_shp, task))
        if res != 0:
            sys.exit(1)

        if b_remove_job_folder:
            os.system('rm -r seg_dem_diff_*')
            io_function.delete_file_or_dir(dem_list_txt)



def main():
    tasks =  ['segment']
    produce_products_dem_subsidence()
    pass

if __name__ == '__main__':
    main()
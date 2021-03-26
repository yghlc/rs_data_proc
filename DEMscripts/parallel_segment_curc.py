#!/usr/bin/env python
# Filename: parallel_segment_curc.py 
"""
introduction: run segmentation for DEM difference on curc

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 March, 2021
"""

import os,sys
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import slurm_utility
import basic_src.io_function as io_function
import basic_src.basic as basic

import re
from datetime import datetime
import time

# some path
import dem_common
from dem_segment_subsidence import segment_subsidence_grey_image

root_dir = os.path.expanduser('/scratch/summit/lihu9680/Arctic/dem_processing')    # run in this folder
jobsh_dir = os.path.expanduser('/scratch/summit/lihu9680/Arctic/dem_processing/dem_proc_jobs')

machine_name = os.uname()[1]
curc_username = 'lihu9680'


def segment_a_dem_diff(dem_diff_path, process_num, ele_diff_thr, min_area, max_area, job_id=0):

    basic.setlogfile('log_segment_dem_diff_job_%d.txt'%job_id)

    # find 8bit one
    tif_8bit = io_function.get_name_by_adding_tail(dem_diff_path, '8bit')
    demD_8bit= os.path.join(dem_common.grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
    if os.path.isfile(demD_8bit):
        basic.outputlogMessage('error, 8bit DEM diff not exists: %s '%demD_8bit)
        return False

    grid_id = int(re.findall('grid\d+',os.path.basename(dem_diff_path))[0][4:])

    save_dir = os.path.join(dem_common.grid_dem_diffs_segment_dir, 'segment_result_grid%d'%grid_id)
    return segment_subsidence_grey_image(demD_8bit, dem_diff_path, save_dir, process_num,
                                  subsidence_thr_m=ele_diff_thr,
                                  min_area=min_area, max_area=max_area)

def seg_working_dir_string(trial_id, root=None):
    if root is not None:
        return os.path.join(root, 'seg_dem_diff_' + str(trial_id).zfill(5))
    else:
        return 'seg_dem_diff_' + str(trial_id).zfill(5)


def copy_curc_seg_job_files(sh_dir, work_dir):
    sh_list = ['seg_dem_diff.sh','job.sh','run_INsingularity_curc_GPU_tf.sh']
    for sh in sh_list:
        io_function.copy_file_to_dst(os.path.join(sh_dir, sh), os.path.join(work_dir, sh)) #, overwrite=True

def submit_segment_dem_diff_job(dem_diff_list, idx):

    while True:
        job_count = slurm_utility.get_submit_job_count(curc_username, job_name_substr='seg')
        if job_count >= 50:
            print(machine_name, datetime.now(),'You have submitted 50 or more jobs, wait ')
            time.sleep(60) #
            continue
        break


    job_name = 'seg%d'%idx
    work_dir = seg_working_dir_string(idx, root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        io_function.save_list_to_txt('dem_diff_list.txt',dem_diff_list)

        # run segmentation
        copy_curc_seg_job_files(jobsh_dir, work_dir)
        slurm_utility.modify_slurm_job_sh('job.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job'%work_dir)
            return


    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    res = os.system('sbatch job.sh' )
    if res != 0:
        sys.exit(1)

    os.chdir(curr_dir_before_start)

    return


def main():

    # write these parameters in the exe.sh file
    # process_num = 8 #  process
    # ele_diff_thr = -2 # meters
    # min_area = 40  # m^2
    # max_area = 10*10*100000 # 10 km by 10 km

    dem_diff_list = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_dir, '*DEM_diff_grid*.tif')
    if len(dem_diff_list) < 1:
        raise ValueError('No *DEM_diff*.tif in %s'%dem_common.grid_dem_diffs_dir)

    print('total %d DEM differnce tifs'%len(dem_diff_list))

    n_tif_per_jobs = 50  # each job, have how many tif to segment
    total_tif_count = len(dem_diff_list)
    # divide to many sub jobs
    tif_groups = [dem_diff_list[i:i + n_tif_per_jobs] for i in range(0, total_tif_count, n_tif_per_jobs)]

    for idx, tif_group in enumerate(tif_groups):

        print(datetime.now(), 'processing %d group of DEM diff, total %d ones'%(idx, len(tif_groups)))
        submit_segment_dem_diff_job(dem_diff_list, idx)


    pass


if __name__ == '__main__':
    curr_dir_before_start = os.getcwd()
    print('\n\ncurrent folder before start: ', curr_dir_before_start, '\n\n')
    main()

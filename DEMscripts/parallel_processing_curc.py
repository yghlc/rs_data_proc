#!/usr/bin/env python
# Filename: parallel_segment_curc.py 
"""
introduction: run segmentation for DEM difference on curc

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 March, 2021
"""

import os,sys
from optparse import OptionParser
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import vector_gpd
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
jobsh_dir = os.path.expanduser('/projects/lihu9680/Data/Arctic/dem_proc_jobs')

machine_name = os.uname()[1]
curc_username = 'lihu9680'

def wait_if_reach_max_jobs(max_job_count,job_name_substr):
    while True:
        job_count = slurm_utility.get_submit_job_count(curc_username, job_name_substr=job_name_substr)
        if job_count >= max_job_count:
            print(machine_name, datetime.now(),'You have submitted %d or more jobs, wait '%max_job_count)
            time.sleep(60) #
            continue
        break

def segment_a_dem_diff(dem_diff_path, process_num, ele_diff_thr, min_area, max_area, job_id=0):

    basic.setlogfile('log_segment_dem_diff_job_%d.txt'%job_id)

    # find 8bit one
    tif_8bit = io_function.get_name_by_adding_tail(dem_diff_path, '8bit')
    demD_8bit= os.path.join(dem_common.grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
    if os.path.isfile(demD_8bit):
        print('error, 8bit DEM diff not exists: %s '%demD_8bit)
        return False

    grid_id = int(re.findall('grid\d+',os.path.basename(dem_diff_path))[0][4:])

    save_dir = os.path.join(dem_common.grid_dem_diffs_segment_dir, 'segment_result_grid%d'%grid_id)
    return segment_subsidence_grey_image(demD_8bit, dem_diff_path, save_dir, process_num,
                                  subsidence_thr_m=ele_diff_thr,
                                  min_area=min_area, max_area=max_area)

def working_dir_string(trial_id, pre_name, root=None):
    if root is not None:
        return os.path.join(root, pre_name + str(trial_id).zfill(5))
    else:
        return pre_name + str(trial_id).zfill(5)


def copy_curc_job_files(sh_dir, work_dir, sh_list):
    for sh in sh_list:
        io_function.copy_file_to_dst(os.path.join(sh_dir, sh), os.path.join(work_dir, sh)) #, overwrite=True

def submit_segment_dem_diff_job(dem_diff_list, idx,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'seg')

    job_name = 'seg%d'%idx
    work_dir = working_dir_string(idx, 'seg_dem_diff_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        io_function.save_list_to_txt('dem_diff_list.txt',dem_diff_list)

        # run segmentation
        sh_list = ['seg_dem_diff.sh', 'job_segment.sh']
        copy_curc_job_files(jobsh_dir, work_dir,sh_list)
        slurm_utility.modify_slurm_job_sh('job_segment.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job'%work_dir)
            return

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done'%work_dir)
            return


    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    res = os.system('sbatch job_segment.sh' )
    if res != 0:
        sys.exit(1)

    os.chdir(curr_dir_before_start)

    return

def run_segment_jobs(max_job_count,n_tif_per_jobs):

    dem_diff_list = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_dir, '*DEM_diff_grid*.tif')
    if len(dem_diff_list) < 1:
        raise ValueError('No *DEM_diff*.tif in %s'%dem_common.grid_dem_diffs_dir)

    print('total %d DEM differnce tifs'%len(dem_diff_list))


    total_tif_count = len(dem_diff_list)
    # divide to many sub jobs
    tif_groups = [dem_diff_list[i:i + n_tif_per_jobs] for i in range(0, total_tif_count, n_tif_per_jobs)]

    for idx, tif_group in enumerate(tif_groups):

        print(datetime.now(), 'processing %d group of segmenting DEM diff, total %d ones'%(idx, len(tif_groups)))
        submit_segment_dem_diff_job(tif_group, idx,max_job_count)

def submit_produce_dem_diff_job(ids_list, idx,grid_base_name,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'demD')

    job_name = 'demD%d'%idx
    work_dir = working_dir_string(idx, 'dem_diff_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        ids_list = [str(item) for item in ids_list]
        io_function.save_list_to_txt(grid_base_name+'.txt',ids_list)

        # run segmentation
        sh_list = ['produce_dem_diff.sh', 'job_dem_diff.sh']
        copy_curc_job_files(jobsh_dir, work_dir,sh_list)
        slurm_utility.modify_slurm_job_sh('job_dem_diff.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job'%work_dir)
            return

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done' % work_dir)
            return


    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    res = os.system('sbatch job_dem_diff.sh' )
    if res != 0:
        sys.exit(1)

    os.chdir(curr_dir_before_start)

    return

def run_grid_jobs(max_job_count,n_tif_per_jobs,task_name,extent_shp):

    from dem_common import grid_20_shp, grid_dem_diffs_dir
    from produce_DEM_diff_ArcticDEM import get_grid_20

    if os.path.isdir(grid_dem_diffs_dir) is  False:
        io_function.mkdir(grid_dem_diffs_dir)

    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    # grid_base_name = os.path.splitext(os.path.basename(extent_shp))[0]
    grid_base_name = 'grid_ids'
    grid_polys, grid_ids = get_grid_20(extent_shp,all_grid_polys, all_ids)

    # divide grid_ids to many groups
    grid_ids_count = len(grid_ids)
    grid_ids_groups = [grid_ids[i:i + n_tif_per_jobs] for i in range(0, grid_ids_count, n_tif_per_jobs)]

    for idx, ids_group in enumerate(grid_ids_groups):
        if task_name=='dem_diff':
            print(datetime.now(), 'processing %d group for DEM diff, total %d ones'%(idx, len(grid_ids_groups)))
            submit_produce_dem_diff_job(ids_group, idx,grid_base_name, max_job_count)
        elif task_name=='dem_headwall_grid':
            print(datetime.now(), 'processing %d group for headwall extraction (per grid), total %d ones'%(idx, len(grid_ids_groups)))
            submit_extract_headwall_grid_job(ids_group, idx,grid_base_name, max_job_count)
        else:
            raise ValueError('unknow task name for grid: %s'%task_name)


def submit_extract_headwall_grid_job(ids_list, idx, grid_base_name,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'gHW')

    job_name = 'gHW%d'%idx
    work_dir = working_dir_string(idx, 'extract_headwall_grid_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        ids_list = [str(item) for item in ids_list]
        io_function.save_list_to_txt(grid_base_name+'.txt',ids_list)

        # run segmentation
        sh_list = ['extract_headwall_from_slope_grid.sh', 'job_healwall_grid.sh']
        copy_curc_job_files(jobsh_dir, work_dir,sh_list)
        slurm_utility.modify_slurm_job_sh('job_healwall_grid.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job'%work_dir)
            return

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done' % work_dir)
            return

    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    res = os.system('sbatch job_healwall_grid.sh' )
    if res != 0:
        sys.exit(1)

    os.chdir(curr_dir_before_start)


def submit_extract_headwall_job(slope_tifs, idx, max_job_count):

    wait_if_reach_max_jobs(max_job_count,'HW')

    job_name = 'HW%d'%idx
    work_dir = working_dir_string(idx, 'extract_headwall_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        io_function.save_list_to_txt('slope_tif_list.txt',slope_tifs)

        # run segmentation
        sh_list = ['job_healwall.sh', 'extract_headwall_from_slope.sh']
        copy_curc_job_files(jobsh_dir, work_dir,sh_list)
        slurm_utility.modify_slurm_job_sh('job_healwall.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done' % work_dir)
            return

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job'%work_dir)
            return

    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    res = os.system('sbatch job_healwall.sh' )
    if res != 0:
        sys.exit(1)

    os.chdir(curr_dir_before_start)

    return


def run_extract_headwall_jobs(max_job_count, n_tif_per_jobs):

    from dem_common import dem_headwall_shp_dir,dem_slope_dir

    if os.path.isdir(dem_headwall_shp_dir) is  False:
        io_function.mkdir(dem_headwall_shp_dir)

    # get slope file list
    slope_tifs = io_function.get_file_list_by_ext('.tif',dem_slope_dir,bsub_folder=False)
    print('Found %d tif in %s'%(len(slope_tifs), dem_slope_dir))

    # divide grid_ids to many groups
    slope_tif_count = len(slope_tifs)
    slope_tif_groups = [slope_tifs[i:i + n_tif_per_jobs] for i in range(0, slope_tif_count, n_tif_per_jobs)]

    for idx, slope_tifs_group in enumerate(slope_tif_groups):

        print(datetime.now(), 'processing %d group for extracting headwall, total %d ones'%(idx, len(slope_tif_groups)))
        submit_extract_headwall_job(slope_tifs_group, idx, max_job_count)



def main(options, args):

    task_name = args[0]
    max_job_count = options.max_job_count
    print('max_job_count', max_job_count)
    n_tif_per_jobs = options.n_tif_per_job  # each job, have how many tif to segment
    extent_shp = options.extent_shp

    if task_name == 'segment':
        run_segment_jobs(max_job_count, n_tif_per_jobs)
    elif task_name == 'dem_diff':
        run_grid_jobs(max_job_count, n_tif_per_jobs,'dem_diff',extent_shp)
    elif task_name == 'dem_headwall_grid':
        run_grid_jobs(max_job_count, n_tif_per_jobs, 'dem_headwall_grid',extent_shp)
    elif task_name == 'dem_headwall':
        run_extract_headwall_jobs(max_job_count, n_tif_per_jobs)
    else:
        print('unknow task name: %s'%task_name)
        pass



if __name__ == '__main__':
    usage = "usage: %prog [options] task_name (dem_diff, segment, dem_headwall) "
    parser = OptionParser(usage=usage, version="1.0 2021-4-25")
    parser.description = 'Introduction:  parallel processing DEM on CURC '

    parser.add_option("-j", "--max_job_count",
                      action="store", dest="max_job_count", type=int, default=50,
                      help="number of jobs to submit at the same time")

    parser.add_option("-n", "--n_tif_per_job",
                      action="store", dest="n_tif_per_job", type=int, default=10,
                      help="number of input files (task) per job")

    parser.add_option("-e", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the extent shapefile")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    curr_dir_before_start = os.getcwd()
    print('\ncurrent folder before start: ', curr_dir_before_start, '\n')

    main(options, args)

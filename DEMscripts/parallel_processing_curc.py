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
from multiprocessing import Process
local_tasks = []
b_run_job_local = False

root_dir = os.path.expanduser('/scratch/summit/lihu9680/Arctic/dem_processing')    # run in this folder
jobsh_dir = os.path.expanduser('/projects/lihu9680/Data/Arctic/dem_proc_jobs')
data_dir = None

machine_name = os.uname()[1]
curc_username = 'lihu9680'

def check_length_jobname(job_name):
    if len(job_name) > 8:
        raise ValueError('the length job name exceed 8 letters, will be cut off to 8, leading to troubles')

def wait_if_reach_max_jobs(max_job_count,job_name_substr):
    if b_run_job_local:
        # in the local machine

        basic.check_exitcode_of_process(local_tasks) # if there is one former job failed, then quit

        while True:
            job_count = basic.alive_process_count(local_tasks)
            if job_count >= max_job_count:
                print(machine_name, datetime.now(),'You are running %d or more tasks in parallel, wait '%max_job_count)
                time.sleep(60) #
                continue
            break
        basic.close_remove_completed_process(local_tasks)

    else:
        while True:
            job_count = slurm_utility.get_submit_job_count(curc_username, job_name_substr=job_name_substr)
            if job_count >= max_job_count:
                print(machine_name, datetime.now(),'You have submitted %d or more jobs, wait '%max_job_count)
                time.sleep(60) #
                continue
            break

def run_a_script(proc_sh, place_holder=None):
    # simulate a SLURM_JOB_ID, which will be used in script for log
    os.environ['SLURM_JOB_ID'] = str(os.getpid()) + '-' + str(int(time.time()))
    res = os.system('./%s'%proc_sh)
    if res != 0:
        sys.exit(1)

def submit_job_curc_or_run_script_local(job_sh, proc_sh):
    """
    submit a job to curc or start running the script locally
    :param job_sh: job to submit
    :param proc_sh:  processing script.
    :return:
    """
    # 'shas0136' 'shas0137' are compile node on CURC
    # if machine_name == 'ubuntu' or machine_name == 'uist-int-colorado-edu' or 'colorado.edu' in machine_name or 'MacBook' in machine_name:
    #     b_run_job_local = True

    if b_run_job_local:
        print(datetime.now(),'will run the job on local machine, not to submit a slurm job')
        # run the job in local computer
        # command_str = './%s'%proc_sh
        # res = os.system(command_str) # this will wait until the job exist

        sub_process = Process(target=run_a_script, args=(proc_sh,None)) # start a process, don't wait
        sub_process.start()
        local_tasks.append(sub_process)
    else:
        # submit the job
        res = os.system('sbatch %s'%job_sh)
        if res != 0:
            sys.exit(1)


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
    check_length_jobname(job_name)
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
    submit_job_curc_or_run_script_local('job_segment.sh','seg_dem_diff.sh')

    os.chdir(curr_dir_before_start)

    return

def run_segment_jobs(max_job_count,n_tif_per_jobs,extent_or_id_txt=None):

    dem_list_txt = os.path.join(dem_common.process_log_dir, 'job_seg_dem_diff_list_' + machine_name + '.txt')
    if os.path.isfile(dem_list_txt) and extent_or_id_txt=='monitor_fail_segment_jobs':
        dem_diff_list = io_function.read_list_from_txt(dem_list_txt)
        print(datetime.now(), 'total %d DEM differnce tifs to seg ' % (len(dem_diff_list)))
    else:
        dem_diff_list = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_dir, '*DEM_diff_grid*.tif')
        dem_diff_list_copy = dem_diff_list.copy()
        if len(dem_diff_list) < 1:
            raise ValueError('No *DEM_diff*.tif in %s'%dem_common.grid_dem_diffs_dir)

        dem_diff_ids = [int(re.findall('grid\d+', os.path.basename(item))[0][4:]) for item in dem_diff_list]

        dem_subsidence_shps = io_function.get_file_list_by_pattern(dem_common.grid_dem_diffs_segment_dir,'segment_result*/*_post.shp')
        subsidence_ids = [int(re.findall('grid\d+', os.path.basename(item))[0][4:]) for item in dem_subsidence_shps]
        no_subsidence_ids = []
        if os.path.isfile(dem_common.grid_no_subscidence_poly_txt):
            no_subsidence_ids = [int(item) for item in io_function.read_list_from_txt(dem_common.grid_no_subscidence_poly_txt) ]
        subsidence_ids.extend(no_subsidence_ids)

        # remove dem diff already been segmented or no subsidence
        new_dem_diff_ids, new_dem_diff_list = [], []
        for id, dem_diff in zip(dem_diff_ids,dem_diff_list_copy):
            if id not in subsidence_ids:
                new_dem_diff_ids.append(id)
                new_dem_diff_list.append(dem_diff)
        dem_diff_ids = new_dem_diff_ids
        dem_diff_list = new_dem_diff_list

        # only keep the ids within in extent
        if extent_or_id_txt is not None:
            new_dem_diff_ids, new_dem_diff_list = [], []
            grid_ids = get_grid_ids_extent(extent_or_id_txt)
            for id, dem_diff in zip(dem_diff_ids, dem_diff_list):
                if id in grid_ids:
                    new_dem_diff_ids.append(id)
                    new_dem_diff_list.append(dem_diff)
            dem_diff_ids = new_dem_diff_ids
            dem_diff_list = new_dem_diff_list


        print('total %d DEM differnce tifs, %d of them need to segment'%(len(dem_diff_list_copy), len(dem_diff_list)))


    total_tif_count = len(dem_diff_list)
    # divide to many sub jobs
    tif_groups = [dem_diff_list[i:i + n_tif_per_jobs] for i in range(0, total_tif_count, n_tif_per_jobs)]

    for idx, tif_group in enumerate(tif_groups):

        print(datetime.now(), 'processing %d group of segmenting DEM diff, total %d ones'%(idx, len(tif_groups)))
        submit_segment_dem_diff_job(tif_group, idx,max_job_count)

def submit_produce_dem_diff_job(ids_list, idx,grid_base_name,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'demD')

    job_name = 'demD%d'%idx
    check_length_jobname(job_name)
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
    submit_job_curc_or_run_script_local('job_dem_diff.sh','produce_dem_diff.sh')

    os.chdir(curr_dir_before_start)

    return

def submit_hillshade_newest_headwall_line_grid_job(ids_list, idx, grid_base_name,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'dLi')   # draw Line on hillshade

    job_name = 'dLi%d'%idx
    check_length_jobname(job_name)
    work_dir = working_dir_string(idx, 'hillshade_newest_headwall_line_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        ids_list = [str(item) for item in ids_list]
        io_function.save_list_to_txt(grid_base_name+'.txt',ids_list)

        # prepare job
        sh_list = ['hillshade_headwall_line_grid.sh', 'job_hillshade_headwall_line_grid.sh']
        copy_curc_job_files(jobsh_dir, work_dir,sh_list)
        slurm_utility.modify_slurm_job_sh('job_hillshade_headwall_line_grid.sh', 'job-name', job_name)
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
    submit_job_curc_or_run_script_local('job_hillshade_headwall_line_grid.sh','hillshade_headwall_line_grid.sh')

    os.chdir(curr_dir_before_start)

from produce_DEM_diff_ArcticDEM import get_grid_20
from dem_common import grid_20_shp

def get_grid_ids_extent(extent_shp):
    if 'ArcticDEM_grid_20km' in os.path.basename(extent_shp):
        print('input %s like a grid files, read grid polygons and ids from it directly'%extent_shp)
        grid_polys, grid_ids = vector_gpd.read_polygons_attributes_list(extent_shp, 'grid_id')
        file_name_base = os.path.splitext(os.path.basename(extent_shp))[0]
        shp_corresponding_grid_ids_txt = file_name_base + '_grid_ids.txt'
        io_function.save_list_to_txt(shp_corresponding_grid_ids_txt,[str(item) for item in grid_ids])
    else:
        # read grids and ids
        time0 = time.time()
        all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')   # in this file, it's "id", not "grid_id"
        print('time cost of read polygons and attributes', time.time() - time0)
        grid_polys, grid_ids = get_grid_20(extent_shp,all_grid_polys, all_ids)

    return grid_ids

def run_grid_jobs(max_job_count,n_tif_per_jobs,task_name,extent_shp):

    from dem_common import grid_dem_diffs_dir

    if os.path.isdir(grid_dem_diffs_dir) is  False:
        io_function.mkdir(grid_dem_diffs_dir)

    # get grid ids based on input extent
    # grid_base_name = os.path.splitext(os.path.basename(extent_shp))[0]
    grid_base_name = 'grid_ids'
    grid_ids = get_grid_ids_extent(extent_shp)

    # divide grid_ids to many groups
    grid_ids_count = len(grid_ids)
    grid_ids_groups = [grid_ids[i:i + n_tif_per_jobs] for i in range(0, grid_ids_count, n_tif_per_jobs)]

    for idx, ids_group in enumerate(grid_ids_groups):
        if task_name=='dem_diff':
            print(datetime.now(), 'processing %d group for DEM diff, total %d ones'%(idx, len(grid_ids_groups)))
            submit_produce_dem_diff_job(ids_group, idx,grid_base_name, max_job_count)
        elif task_name=='hillshade_headwall_line':
            print(datetime.now(),'processing %d group for hillshade newest and headwall Line (per grid), total %d ones' % (idx, len(grid_ids_groups)))
            submit_hillshade_newest_headwall_line_grid_job(ids_group, idx,grid_base_name, max_job_count)
        elif task_name=='dem_headwall_grid':
            print(datetime.now(), 'processing %d group for headwall extraction (per grid), total %d ones'%(idx, len(grid_ids_groups)))
            submit_extract_headwall_grid_job(ids_group, idx,grid_base_name, max_job_count)
        else:
            raise ValueError('unknow task name for grid: %s'%task_name)


def submit_extract_headwall_grid_job(ids_list, idx, grid_base_name,max_job_count):

    wait_if_reach_max_jobs(max_job_count,'gHW')

    job_name = 'gHW%d'%idx
    check_length_jobname(job_name)
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
    submit_job_curc_or_run_script_local('job_healwall_grid.sh','extract_headwall_from_slope_grid.sh')

    os.chdir(curr_dir_before_start)


def submit_extract_headwall_job(slope_tifs, idx, max_job_count):

    wait_if_reach_max_jobs(max_job_count,'HW')

    job_name = 'HW%d'%idx
    check_length_jobname(job_name)
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
    submit_job_curc_or_run_script_local('job_healwall.sh','extract_headwall_from_slope.sh')

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

def submit_headwall_ripple_job(headwall_shp_list, idx,max_job_count):

    wait_if_reach_max_jobs(max_job_count, 'hwR')

    job_name = 'hwR%d' % idx
    check_length_jobname(job_name)
    work_dir = working_dir_string(idx, 'headwall_ripple_', root=root_dir)
    if os.path.isdir(work_dir) is False:
        io_function.mkdir(work_dir)
        os.chdir(work_dir)

        io_function.save_list_to_txt('headwall_shp_list.txt', headwall_shp_list)

        # run segmentation
        sh_list = ['headwall_ripple_stastics.sh', 'job_headwall_ripple.sh','setting.ini']
        copy_curc_job_files(jobsh_dir, work_dir, sh_list)
        slurm_utility.modify_slurm_job_sh('job_headwall_ripple.sh', 'job-name', job_name)

    else:
        os.chdir(work_dir)

        submit_job_names = slurm_utility.get_submited_job_names(curc_username)
        if job_name in submit_job_names:
            print('The folder: %s already exist and the job has been submitted, skip submitting a new job' % work_dir)
            return

        # job is completed
        if os.path.isfile('done.txt'):
            print('The job in the folder: %s is Done' % work_dir)
            return

    # submit the job
    # sometime, when submit a job, end with: singularity: command not found,and exist, wired, then try run submit a job in scomplie note
    submit_job_curc_or_run_script_local('job_headwall_ripple.sh', 'headwall_ripple_stastics.sh')

    os.chdir(curr_dir_before_start)

    return

def run_headwall_ripple_jobs(max_job_count,n_tif_per_jobs,extent_or_id_txt=None):

    shp_list_txt = os.path.join(dem_common.process_log_dir, 'job_headwall_ripple_list_' + machine_name + '.txt')
    if os.path.isfile(shp_list_txt) and extent_or_id_txt == 'monitor_fail_headwall_ripple_jobs':
        headwall_shps = io_function.read_list_from_txt(shp_list_txt)
        print(datetime.now(), 'total %d headwall ripple shps to run ' % (len(headwall_shps)))
    else:
        grid_dem_headwall_shp_dir = data_dir
        headwall_shps = io_function.get_file_list_by_pattern(grid_dem_headwall_shp_dir, 'headwall_shps_grid*/*.shp')
        # remove rippleSel
        headwall_shps = [ item for item in headwall_shps if 'rippleSel' not in os.path.basename(item)]

        # headwall_shps_copy = headwall_shps.copy()
        if len(headwall_shps) < 1:
            raise ValueError('No headwall_shps_grid*/*.shp in %s' % grid_dem_headwall_shp_dir)

    total_count = len(headwall_shps)
    # divide to many sub jobs
    tif_groups = [headwall_shps[i:i + n_tif_per_jobs] for i in range(0, total_count, n_tif_per_jobs)]

    for idx, tif_group in enumerate(tif_groups):
        print(datetime.now(), 'processing %d group of headwall ripple, total %d ones' % (idx, len(tif_groups)))
        submit_headwall_ripple_job(tif_group, idx, max_job_count)



def main(options, args):

    task_name = args[0]
    max_job_count = options.max_job_count
    print('max_job_count', max_job_count)
    n_tif_per_jobs = options.n_tif_per_job  # each job, have how many tif to segment
    extent_shp_or_id_txt = options.extent_shp
    global b_run_job_local
    if options.b_run_job_local:
        b_run_job_local = True

    if options.user_name is not None:
        global curc_username
        curc_username = options.user_name
    if options.working_dir is not None:
        global root_dir
        root_dir = options.working_dir
    if options.script_dir is not None:
        global jobsh_dir
        jobsh_dir = options.script_dir
    if options.data_dir is not None:
        global data_dir
        data_dir = options.data_dir

    if task_name == 'segment':
        run_segment_jobs(max_job_count, n_tif_per_jobs,extent_or_id_txt=extent_shp_or_id_txt)
    elif task_name == 'dem_diff':
        run_grid_jobs(max_job_count, n_tif_per_jobs,'dem_diff',extent_shp_or_id_txt)
    elif task_name == 'hillshade_headwall_line':
        run_grid_jobs(max_job_count, n_tif_per_jobs, 'hillshade_headwall_line', extent_shp_or_id_txt)
    elif task_name == 'dem_headwall_grid':
        run_grid_jobs(max_job_count, n_tif_per_jobs, 'dem_headwall_grid',extent_shp_or_id_txt)
    elif task_name == 'dem_headwall':
        run_extract_headwall_jobs(max_job_count, n_tif_per_jobs)
    elif task_name == 'headwall_ripple':
        run_headwall_ripple_jobs(max_job_count, n_tif_per_jobs,extent_or_id_txt=extent_shp_or_id_txt)
    else:
        print('unknow task name: %s'%task_name)
        pass

    time.sleep(10)

    # wait all local task finished
    while basic.b_all_process_finish(local_tasks) is False:
        print(datetime.now(),'wait 5 minutes to let all local tasks to complete')
        time.sleep(60*5)


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
                      help="the extent shapefile or grid_id_list in txt")

    parser.add_option("-u", "--user_name",
                      action="store", dest="user_name",
                      help="the username of the server")

    parser.add_option("", "--b_run_job_local",
                      action="store_true", dest="b_run_job_local",default=False,
                      help="if set (True), will run the job on the machine instead of submitting a slurm job")

    parser.add_option("", "--working_dir",
                      action="store", dest="working_dir",
                      help="the working directory")
    parser.add_option("", "--script_dir",
                      action="store", dest="script_dir",
                      help="the directory that template scripts ")
    parser.add_option("", "--data_dir",
                      action="store", dest="data_dir",
                      help="the directory for data ")


    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    curr_dir_before_start = os.getcwd()
    print('\ncurrent folder before start: ', curr_dir_before_start, '\n')

    main(options, args)

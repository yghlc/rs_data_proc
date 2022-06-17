#!/usr/bin/env python
# Filename: process_largeRegion_butLimited_storage.py
"""
introduction: to process large areas of ArcticDEM but with limited storage,
we will process grid by grid.

We will first download ArcticDEM strip for needed grids, then remove the strip file if the grids it covered have been processed

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 11 October, 2021
"""

import os,sys
from optparse import OptionParser
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import vector_gpd
import slurm_utility
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import raster_io
import scp_communicate

import re
from datetime import datetime
import time
import numpy as np
import pandas as pd

import multiprocessing
from multiprocessing import Pool


machine_name = os.uname()[1]

# some files
from dem_common import grid_20_shp,grid_20_id_raster,dem_strip_shp,dem_tile_shp

# log or txt files
from dem_common import process_log_dir, grid_complete_list_txt, grid_excluded_list_txt,strip_dem_cover_grids_txt, tile_dem_cover_grids_txt
from dem_common import grid_no_dem_txt,grid_no_valid_dem_txt,grid_dem_diff_less2dem_txt,grid_no_headwall_txt
if os.path.isdir(process_log_dir) is False:
    io_function.mkdir(process_log_dir)

from dem_common import check_create_lock, release_lock

from dem_common import arcticDEM_tile_tarball_dir,arcticDEM_tile_reg_tif_dir,tarball_dir,arcticDEM_reg_tif_dir,ArcticDEM_tmp_dir

# results dir
from dem_common import grid_hillshade_newest_HDLine_dir, grid_dem_diffs_dir,grid_dem_headwall_shp_dir

from dem_common import grid_dem_diffs_segment_dir, grid_no_subscidence_poly_txt, arcticDEM_tile_slope_dir

from dem_common import get_corresponding_grid_ids_txt,grid_ids_txt_dir,get_extent_grid_id_txt_done_files

from produce_DEM_diff_ArcticDEM import get_grid_20

from parallel_processing_curc import curc_username

# scripts
dem_download_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/download_arcticDEM.py')
dem_unpack_reg_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/ArcticDEM_unpack_registration.py')
subset_shp_dir = 'subset_grids_shp'
msg_file_pre = 'subset'     # prename of the message file
download_ahead_proc = 4
from dem_common import subset_message_dir
if os.path.isdir(subset_message_dir) is False:
    io_function.mkdir(subset_message_dir)

def update_subset_info(txt_path, key_list=None, info_list=None):
    # maintain a info of subset for processing, dict
    # id: subset_id
    # shp: the shapefile contain all grids in this subset
    # "pre_status": the status of downloading and registration of ArcticDEM, has values: 'notYet', 'working', 'done'
    # 'proc_status': the status of processing ArcticDEM, has values of 'notYet', 'working', 'done'

    info_dict = {}
    if os.path.isfile(txt_path):
        info_dict = io_function.read_dict_from_txt_json(txt_path)
    if isinstance(key_list,str):
        key_list = [key_list]
    if isinstance(info_list,str):
        info_list = [info_list]
    if len(key_list) != len(info_list):
        raise ValueError('the lengths of keys and info are different')
    for key, info in zip(key_list, info_list):
        info_dict[key] = info
    io_function.save_dict_to_txt_json(txt_path,info_dict)

def get_subset_info(txt_path):
    info_dict = io_function.read_dict_from_txt_json(txt_path)
    # change file name to absolute path depending on machine
    info_dict['shp'] = os.path.join(subset_shp_dir,info_dict['shp'])
    return info_dict

def get_subset_info_txt_list(key,values,remote_node=None, remote_folder=None, local_folder='./'):
    # get subset info with specific key and values
    # 'subset%d.txt' % subset_id
    # sync 'subset%d.txt' files from remote node
    if remote_node is not None:
        # this will overwrite 'subset%d.txt' in local
        scp_communicate.copy_file_folder_from_remote_machine(remote_node,os.path.join(remote_folder,msg_file_pre+'*.txt'),
                                                             os.path.join(local_folder,'.'))

    select_txt_list = []
    subinfo_txt_list = io_function.get_file_list_by_pattern(local_folder,msg_file_pre+'*.txt')
    for txt in subinfo_txt_list:
        info_dict = get_subset_info(txt)
        if info_dict[key] in values:
            select_txt_list.append(txt)
    return sorted(select_txt_list)


def download_process_send_arctic_dem(subset_info_txt, r_working_dir, remote_node,task_list, b_send_data=True):
    # this function run on tesia

    subset_info = get_subset_info(subset_info_txt)

    if subset_info['pre_status'] == 'done':
        print(datetime.now(),'pre_status for %s is done, skip'%subset_info_txt)
        # copy to remote machine
        if b_send_data:
            if scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir, subset_info_txt):
                scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir, subset_shp_dir)
            return True

    # if subset_id for download is far more ahead than processing (curc), then wait, in case occupy too much storage
    while True and b_send_data:
        remote_sub_txt = get_subset_info_txt_list('proc_status',['notYet', 'working'],remote_node=remote_node,remote_folder=r_working_dir)
        if len(remote_sub_txt) > download_ahead_proc:
            print(datetime.now(),'there is %d subset have not complete,'
                                 ' wait 300 seconds, avoid downloading too many data'%len(remote_sub_txt))
            time.sleep(300)
        else:
            break

    update_subset_info(subset_info_txt, key_list=['pre_status','pre_node'], info_list=['working',machine_name])

    # download strip tarball, upack, registration
    res = os.system(dem_download_py + ' ' + subset_info['shp'] + ' ' + dem_strip_shp )
    if res != 0:
        sys.exit(1)
    # res = os.system(dem_unpack_reg_py + ' ' + tarball_dir  + ' -d ' + arcticDEM_reg_tif_dir + ' -r ' )  # -r: --remove_inter_data
    # if res != 0:
    #     sys.exit(1)

    # download tile tarball and unpack
    res = os.system(dem_download_py + ' ' + subset_info['shp'] + ' ' + dem_tile_shp )
    if res != 0:
        sys.exit(1)
    # res = os.system(dem_unpack_reg_py + ' ' +arcticDEM_tile_tarball_dir + ' -d ' + arcticDEM_tile_reg_tif_dir + ' -r ' )
    # if res != 0:
    #     sys.exit(1)
    # parallel unpack and registration for Strip and Tile DEM
    res = os.system('./parallel_unpack.sh')
    if res != 0:
        sys.exit(1)

    # ArcticDEM (mosaic) to slope, needed for segmentation of dem diff
    if 'segment' in task_list:
        res = os.system('./dem_to_slope.sh')
        if res != 0:
            sys.exit(1)



    # send to remote machine
    rsync_sh = os.path.join(ArcticDEM_tmp_dir,'rsync_to_curc.sh')
    if b_send_data:
        # create a lock file (make sure only one workstation is sending data to curc)
        rsync_to_curc_lock = os.path.join(ArcticDEM_tmp_dir, 'rsync_to_curc_lock.txt')
        check_create_lock(rsync_to_curc_lock, 'because other program is sending data (rsync_to_curc.sh)')

        res = os.system(rsync_sh)
        if res != 0:
            sys.exit(1)
        release_lock(rsync_to_curc_lock)

    update_subset_info(subset_info_txt,key_list=['pre_status','pre_done_time'],info_list=['done',str(datetime.now())])
    # copy to remote machine
    if b_send_data:
        if scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir, subset_info_txt):
            scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir,subset_shp_dir)


def copy_results_from_remote_node():
    # rsync results from remote machine (it's ok if fil)
    rsync_sh = os.path.join(ArcticDEM_tmp_dir,'rsync_from_curc.sh')
    res = os.system(rsync_sh)
    if res != 0:
        sys.exit(1)




def get_overlap_grids_for_one_extent(all_ids,all_grid_polys, dem_poly, dem_name, idx, dem_poly_count):
    print(timeTools.get_now_time_str(), idx, dem_poly_count)
    index = vector_gpd.get_poly_index_within_extent(all_grid_polys, dem_poly)
    gird_ids = [all_ids[idx] for idx in index]
    return dem_name,gird_ids

def build_dict_of_dem_cover_grid_ids(dem_info_shp,grid_20_shp,save_dict_txt):
    # this will take time, but only need to run once at the beginning
    if os.path.isfile(save_dict_txt):
        print('warning, %s exists, skip build_dict_of_dem_cover_grid_ids'%save_dict_txt)
        return True

    # extent polygons and projection (proj4)
    dem_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(dem_info_shp)
    if dem_shp_prj == '':
        raise ValueError('get proj4 of %s failed' % dem_info_shp)
    grid_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)
    if grid_shp_prj == '':
        raise ValueError('get proj4 of %s failed' % grid_20_shp)

    if dem_shp_prj != grid_shp_prj:
        raise ValueError('%s and %s do not have the same projection'% (dem_info_shp, grid_20_shp))

    # read DEM info
    dem_polygons, dem_names = vector_gpd.read_polygons_attributes_list(dem_info_shp, 'name',b_fix_invalid_polygon=False)
    # dem_name: eg. SETSM_GE01_20090818_1050410001E0CF00_1050410001D80200_seg1_2m_v3.0  or 11_27_2_1_2m_v3.0
    dem_poly_count = len(dem_polygons)
    # check if there is duplicate dem names
    if len(dem_names) != len(set(dem_names)):
        raise ValueError('some duplicate dem name in %s'%dem_info_shp)

    # read grid polygons and ids
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')

    dem_cover_grids = {}
    # this will take time.
    # for idx, (dem_poly,dem_name) in enumerate(zip(dem_polygons, dem_names)):
    #     print(timeTools.get_now_time_str(), idx, dem_poly_count)
    #     index = vector_gpd.get_poly_index_within_extent(all_grid_polys, dem_poly)
    #     gird_ids = [ all_ids[idx] for idx in index ]
    #     # if dem_name in dem_cover_grids.keys():
    #     #     basic.outputlogMessage('\n Warning, %s already in dict \n'%dem_name)
    #     dem_cover_grids[dem_name] = gird_ids

    ### parallel version
    theadPool = Pool(multiprocessing.cpu_count())  # multi processes
    parameters_list = [(all_ids,all_grid_polys, dem_poly,dem_name, idx, dem_poly_count) for idx, (dem_poly,dem_name) in enumerate(zip(dem_polygons, dem_names))]
    results = theadPool.starmap(get_overlap_grids_for_one_extent, parameters_list)  # need python3
    for res in results:
        dem_name, gird_ids = res
        dem_cover_grids[dem_name] = gird_ids

    # save to dict
    io_function.save_dict_to_txt_json(save_dict_txt,dem_cover_grids)
    theadPool.close()
    return True

def test_build_dict_of_dem_cover_grid_ids():
    # for test, only have three polygons
    dem_info_shp = os.path.expanduser('~/Data/tmp_data/test_ArcticDEM/ArcticDEM_Tile_Index_Rel7.shp')
    save_dict_txt = tile_dem_cover_grids_txt
    build_dict_of_dem_cover_grid_ids(dem_info_shp, grid_20_shp, save_dict_txt)

def get_not_completed_grids(grid_polys, grid_ids):
    completed_id_list = []
    if os.path.isfile(grid_complete_list_txt):
        completed_id_list =  [int(item) for item in io_function.read_list_from_txt(grid_complete_list_txt)]

    if len(completed_id_list) < 1:
        return grid_polys, grid_ids

    no_complete_polys =[]
    no_complete_ids = []
    for poly, id in zip(grid_polys, grid_ids):
        if id in completed_id_list:
            pass
        else:
            no_complete_polys.append(poly)
            no_complete_ids.append(id)
    return no_complete_polys, no_complete_ids

def get_complete_ignore_grid_ids():
    id_list = []
    if os.path.isfile(grid_complete_list_txt):
        completed_id_list =  [int(item) for item in io_function.read_list_from_txt(grid_complete_list_txt)]
        id_list.extend(completed_id_list)
    # get manual excluded ids
    if os.path.isfile(grid_excluded_list_txt):
        exclude_id_list = [int(item) for item in io_function.read_list_from_txt(grid_excluded_list_txt)]
        id_list.extend(exclude_id_list)
    # get ids that don't have DEM
    if os.path.isfile(grid_no_dem_txt):
        nodem_id_list = [int(item) for item in io_function.read_list_from_txt(grid_no_dem_txt)]
        id_list.extend(nodem_id_list)
    # get ids that don't have DEM
    if os.path.isfile(grid_no_valid_dem_txt):
        no_valid_dem_id_list = [int(item) for item in io_function.read_list_from_txt(grid_no_valid_dem_txt)]
        id_list.extend(no_valid_dem_id_list)

    return id_list

def save_grid_ids_need_to_process(grid_ids,ignore_ids=None, save_path='grid_ids_to_process.txt'):
    ''' save a list to txt, contain grid ids need to process, return the number of grids to process'''
    if ignore_ids is None:
        id_list = get_complete_ignore_grid_ids()
    else:
        id_list = ignore_ids
    ids_need_to_proc = [str(id) for id in grid_ids if id not in id_list]
    io_function.save_list_to_txt(save_path,ids_need_to_proc)
    return len(ids_need_to_proc)


def b_exist_grid_headwall_shp(id):

    # if a grid don't have headwall (due to really flat or other reasons), then we think it's complete
    if os.path.isfile(grid_no_headwall_txt):
        grid_ids_no_headwall = [int(item) for item in io_function.read_list_from_txt(grid_no_headwall_txt)]
        if id in grid_ids_no_headwall:
            return True

    headwall_shps_dir = io_function.get_file_list_by_pattern(grid_dem_headwall_shp_dir, '*_grid%d' % id)
    if len(headwall_shps_dir) == 1:
        return True
    elif len(headwall_shps_dir) > 1:
        basic.outputlogMessage('warning, There are multiple headwall shps for grid: %d' % id)
        for item in headwall_shps_dir:
            basic.outputlogMessage(item)
        return True
    else:
        return False

def b_exist_dem_hillshade_newest_HWLine_grid(id):

    # if a grid don't have headwall (due to really flat or other reasons), then we think it's complete
    if os.path.isfile(grid_no_headwall_txt):
        grid_ids_no_headwall = [int(item) for item in io_function.read_list_from_txt(grid_no_headwall_txt)]
        if id in grid_ids_no_headwall:
            return True

    hillshade_newest_HDLine_tifs = io_function.get_file_list_by_pattern(grid_hillshade_newest_HDLine_dir, '*_grid%d.tif' % id)
    if len(hillshade_newest_HDLine_tifs) == 1:
        return True
    elif len(hillshade_newest_HDLine_tifs) > 1:
        basic.outputlogMessage('warning, There are multiple hillshade (newest) HDLine tif for grid: %d' % id)
        for item in hillshade_newest_HDLine_tifs:
            basic.outputlogMessage(item)
        return True
    else:
        return False

def b_exist_gid_dem_diff(id):
    # if an id don't have enough dem for dem diff, then we think it's complete
    if os.path.isfile(grid_dem_diff_less2dem_txt):
        grid_ids_less_2dem = [int(item) for item in io_function.read_list_from_txt(grid_dem_diff_less2dem_txt)]
        if id in grid_ids_less_2dem:
            return True

    dem_diff_files = io_function.get_file_list_by_pattern(grid_dem_diffs_dir, '*_DEM_diff_grid%d.tif' % id)

    if len(dem_diff_files) == 1:
        return True
    elif len(dem_diff_files) > 1:
        basic.outputlogMessage('warning, There are multiple DEM difference for grid: %d' % id)
        for item in dem_diff_files:
            basic.outputlogMessage(item)
            return True
    else:
        return False

def b_exist_grid_dem_subsidence(id):
    # if an grid don't have dem subsidence, then we think it's complete
    if os.path.isfile(grid_no_subscidence_poly_txt):
        grid_ids_no_subsidence = [int(item) for item in io_function.read_list_from_txt(grid_no_subscidence_poly_txt)]
        if id in grid_ids_no_subsidence:
            return True

    # if an id don't have enough dem for dem diff, then we think it's complete
    if os.path.isfile(grid_dem_diff_less2dem_txt):
        grid_ids_less_2dem = [int(item) for item in io_function.read_list_from_txt(grid_dem_diff_less2dem_txt)]
        if id in grid_ids_less_2dem:
            return True

    dem_subsidence_shps = io_function.get_file_list_by_pattern(grid_dem_diffs_segment_dir, '*_grid%d/*_grid%d_8bit_post.shp' % (id,id))

    if len(dem_subsidence_shps) == 1:
        return True
    elif len(dem_subsidence_shps) > 1:
        basic.outputlogMessage('warning, There are multiple DEM subsidence for grid: %d' % id)
        for item in dem_subsidence_shps:
            basic.outputlogMessage(item)
            return True
    else:
        return False


def update_complete_grid_list(grid_ids, task_list):
    # based on some criteria, to check if results exist, then update grid_complete_list_txt
    completed_id_list = []
    if os.path.isfile(grid_complete_list_txt):
        completed_id_list = [int(item) for item in io_function.read_list_from_txt(grid_complete_list_txt)]
    n_task = len(task_list)
    if n_task < 1:
        raise ValueError('No task in %s'%str(task_list))

    for g_id in grid_ids:
        if g_id in completed_id_list:
            continue
        # check if it has been completed based on multiple criteria
        complete_count = 0
        if 'dem_diff' in task_list and b_exist_gid_dem_diff(g_id):
            complete_count += 1
        if 'hillshade_headwall_line' in task_list and b_exist_dem_hillshade_newest_HWLine_grid(g_id):
            complete_count += 1
        if 'dem_headwall_grid' in task_list and b_exist_grid_headwall_shp(g_id):
            complete_count += 1
        if 'segment' in task_list and b_exist_grid_dem_subsidence(g_id):
            complete_count += 1
        # we may check more task results: segment, dem_headwall

        if complete_count == n_task:
            completed_id_list.append(g_id)

    # save the txt
    completed_id_list = [str(item) for item in completed_id_list]
    io_function.save_list_to_txt(grid_complete_list_txt,completed_id_list)


def find_neighbours_2d(grid_ids_2d,visited,seed,connect):
    '''
    find neighbourhood voxels
    :param grid_ids_2d: 2D data
    :param visited: indicates pixels has been checked
    :param seed: a seed
    :param connect: pixel connectivity
    :return: a list containing voxels
    '''

    height, width = grid_ids_2d.shape
    y,x = seed[0],seed[1]
    visited[y, x] = 1

    neigh_range = [-1,0,1]
    neighbours = [[i,j] for i in neigh_range for j in neigh_range  ]
    neighbours.remove([0,0])

    # distance within 1
    if connect==4:
        connectivity =  [ [y+dy, x+dx] for (dy,dx) in neighbours if (dy*dy + dx*dx) <= 1 ]
    # distance within sqrt(2)
    elif connect==8:
        connectivity = [[y + dy, x + dx] for (dy, dx) in neighbours if (dy * dy + dx * dx) <= 2]
    else:
        raise ValueError('Only accept connectivity of 4 or 8')

    new_seeds = []
    for [y,x] in connectivity:
        # out extent
        if y<0 or x<0 or y >= height or x >= width:
            continue
        # already visited
        if visited[y,x]:
            continue
        new_seeds.append([y,x])

        # masked as visited
        visited[y,x] = 1

    return new_seeds


def get_grids_for_download_process(grid_polys, grid_ids, ignore_ids,max_grid_count, grid_ids_2d, visit_np, save_path, proj=None):
    '''
    get grids for donwload ArcticDEM and ids
    :param grid_polys:
    :param grid_ids:
    :param ignore_ids: ids need to be ignored (complete or excluded)
    :param max_grid_count: max grid count, it will return number close to this value
    :param grid_ids_2d: gird_ids in 2d array
    :param visit_np: visit_np like a mask, to indicate which pixels has been checked, 0 no visit previously, 1 visited
    :param proj: projection information, for save grid_poygons
    :return:
    '''

    # find a connected region with for donwload and processing, and save to files
    seed_loc = np.where(visit_np == 0)
    if len(seed_loc[0]) < 1:
        print('warning, all pixels have been visited')
        return None, None
    # seed_loc = np.where(grid_ids_2d == grid_ids[0])
    y, x = seed_loc[0][0], seed_loc[1][0]
    selected_gird_id_list = [grid_ids_2d[y, x]]
    seed_list = [ [y, x]]
    while len(selected_gird_id_list) < max_grid_count and len(seed_list) > 0:
        # find neighbours
        new_seeds = find_neighbours_2d(grid_ids_2d,visit_np,seed_list[0],8)
        del seed_list[0]
        seed_list.extend(new_seeds)
        # find new ids
        for seed in new_seeds:
            row, col = seed
            selected_gird_id_list.append( grid_ids_2d[row, col])

    # remove some ids
    # selected_gird_id_list = [id for id in selected_gird_id_list if id not in ignore_ids]
    ignore_ids_in_selected = list(set(ignore_ids).intersection(selected_gird_id_list))  # intersection is faster
    _ = [selected_gird_id_list.remove(rm_id) for rm_id in ignore_ids_in_selected ]
    if len(selected_gird_id_list) < 1:
        return [], []

    select_grid_polys = [ grid_polys[grid_ids.index(item) ] for item in selected_gird_id_list ]

    save_selected_girds_and_ids(selected_gird_id_list,select_grid_polys,proj,save_path)


    return select_grid_polys, selected_gird_id_list


def save_selected_girds_and_ids(selected_gird_id_list,select_grid_polys,proj,save_path):
    # save to shapefile to download and processing
    # change numpy.uint16 to int, avoid become negative when saving to shapefile
    selected_gird_id_list = [int(item) for item in selected_gird_id_list]
    save_pd = pd.DataFrame({'grid_id':selected_gird_id_list, 'Polygon':select_grid_polys})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',proj,save_path)
    basic.outputlogMessage('saved %d grids to %s'%(len(select_grid_polys), save_path))
    # save the ids to txt
    save_id_txt = os.path.splitext(save_path)[0] + '_grid_ids.txt'
    selected_grid_ids_str = [str(item) for item in selected_gird_id_list]
    io_function.save_list_to_txt(save_id_txt, selected_grid_ids_str)


def save_list_no_need_dem_files(file_name,file_list):
    if len(file_list) < 1:
        return True
    # update the file list
    save_list = []
    if os.path.isfile(file_name):
        save_list = io_function.read_list_from_txt(file_name)
    for item in file_list:
        if item in save_list:
            continue
        save_list.append(item)
    return io_function.save_list_to_txt(file_name,save_list)

def remove_no_need_dem_files(b_remove=True):
    # if os.path.isfile(grid_complete_list_txt):
    #     completed_id_list =  [int(item) for item in io_function.read_list_from_txt(grid_complete_list_txt)]
    # else:
    #     print(datetime.now(), 'no complete grids')
    #     return True
    #
    # if os.path.isfile(grid_excluded_list_txt):
    #     exclude_id_list = [int(item) for item in io_function.read_list_from_txt(grid_excluded_list_txt)]
    #     completed_id_list.extend(exclude_id_list)

    completed_id_list = get_complete_ignore_grid_ids()
    if len(completed_id_list) < 1:
        print(datetime.now(), 'no complete grids')
        return True

    if len(completed_id_list) < 1:
        return True

    completed_id_set = set(completed_id_list)

    # check four folders: arcticDEM_tile_tarball_dir,arcticDEM_tile_reg_tif_dir,tarball_dir,arcticDEM_reg_tif_dir
    strip_dem_cover_grids = io_function.read_dict_from_txt_json(strip_dem_cover_grids_txt)

    strip_no_need_list = [strip for strip in strip_dem_cover_grids.keys()
                          if set(strip_dem_cover_grids[strip]).issubset(completed_id_set) ]

    tile_dem_cover_grids = io_function.read_dict_from_txt_json(tile_dem_cover_grids_txt)
    tile_no_need_list = [tile for tile in tile_dem_cover_grids.keys() if
                         set(tile_dem_cover_grids[tile]).issubset(completed_id_set)]

    if b_remove is False:
        save_list_no_need_dem_files('no_need_ArcticDEM_strip_names.txt',strip_no_need_list)
        save_list_no_need_dem_files('no_need_ArcticDEM_mosaic_names.txt',tile_no_need_list)
    else:
        # remove
        basic.outputlogMessage('there are %d no need strip DEM, downloaded files will be or have been removed'%len(strip_no_need_list))
        for strip in strip_no_need_list:
            file_list = io_function.get_file_list_by_pattern(tarball_dir,strip+'*')
            file_list_2 = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,strip+'*')
            file_list.extend(file_list_2)
            if len(file_list) > 0:
                for path in file_list:
                    basic.outputlogMessage('removing %s' % path)
                    io_function.delete_file_or_dir(path)


        basic.outputlogMessage('there are %d no need tile DEM, downloaded files will be or have been removed'%len(tile_no_need_list))
        for tile in tile_no_need_list:
            file_list = io_function.get_file_list_by_pattern(arcticDEM_tile_tarball_dir,tile+'*')
            file_list_2 = io_function.get_file_list_by_pattern(arcticDEM_tile_reg_tif_dir,tile+'*')
            file_list.extend(file_list_2)
            # remove slope file derived ArcticDEM (mosaic)
            file_list_3 = io_function.get_file_list_by_pattern(arcticDEM_tile_slope_dir, tile + '*')
            file_list.extend(file_list_3)
            if len(file_list) > 0:
                for path in file_list:
                    basic.outputlogMessage('removing %s' % path)
                    io_function.delete_file_or_dir(path)

no_subset_to_proc = 0
def produce_dem_products(tasks,b_remove_job_folder=True,b_remove_dem=True,no_slurm=False,message_dir='./'):
    # this function run on process node, such as curc
    global no_subset_to_proc

    subset_txt_list = get_subset_info_txt_list('proc_status',['notYet', 'working'],local_folder=message_dir)
    if len(subset_txt_list) < 1:
        print(datetime.now(), 'checking times: %d: No subset (%s) for processing, wait 300 seconds'%(no_subset_to_proc, msg_file_pre+'*.txt'))
        time.sleep(300)
        no_subset_to_proc += 1
        if no_subset_to_proc > 60:   # if has continued waiting from 6o times (10 hours), then return Flase, will exit the while loop
            return False
        return True

    no_subset_to_proc = 0  # reset the count if it found a job to process
    # task_job_names are from parallel_processing_curc.py
    task_job_name = {'dem_diff':'demD',
                   'dem_headwall_grid':'gHW',
                   'hillshade_headwall_line':'dLi',
                   'segment':'seg'}
    task_depend = {'dem_diff':[],
                   'dem_headwall_grid':[],
                   'hillshade_headwall_line':['dem_headwall_grid'],
                   'segment':['dem_diff']}


    subset_txt_list = sorted(subset_txt_list)
    for sub_txt in subset_txt_list:
        update_subset_info(sub_txt,key_list=['proc_status'],info_list=['working'])
        subset_info = get_subset_info(sub_txt)
        ext_shp = subset_info['shp']

        # submit tasks with dependencies
        tasks_no_depend = [item for item in tasks if len(task_depend[item]) < 1 ]
        for task in tasks_no_depend:
            if no_slurm:
                res = os.system('./run_local.sh %s %s' % (ext_shp, task))
            else:
                res = os.system('./run.sh %s %s'%(ext_shp,task))
            if res !=0:
                sys.exit(1)
        time.sleep(5)   # wait

        # submit tasks with dependencies
        tasks_with_depend = [item for item in tasks if len(task_depend[item]) > 0 ]
        while len(tasks_with_depend) > 0:
            for task in tasks_with_depend:
                depend_tasks = task_depend[task]
                job_count_list = [ slurm_utility.get_submit_job_count(curc_username, job_name_substr=task_job_name[item])
                                   for item in depend_tasks ]
                if sum(job_count_list) > 0:
                    print(machine_name, datetime.now(),
                          'task: %s need results of task:%s whose jobs are not completed, need to wait'%(task,str(depend_tasks)))
                else:
                    if no_slurm:
                        res = os.system('./run_local.sh %s %s' % (ext_shp, task))
                    else:
                        res = os.system('./run.sh %s %s' % (ext_shp, task))
                    if res != 0:
                        sys.exit(1)
                    tasks_with_depend.remove(task)  # if submit, remove it

            time.sleep(60)

        # wait until all jobs finished
        while True:
            if slurm_utility.get_submit_job_count(curc_username, job_name_substr=None) > 0:
                print(machine_name, datetime.now(),'wait 300 seconds until all submitted jobs are completed ')
                time.sleep(300)
                continue

            if b_remove_job_folder:
                # remove temporal folders
                if 'dem_diff' in tasks:
                    os.system('rm -r dem_diff_*')
                if 'segment' in tasks:
                    os.system('rm -r seg_dem_diff_*')
                if 'dem_headwall_grid' in tasks:
                    os.system('rm -r extract_headwall_grid_*')
                if 'hillshade_headwall_line' in tasks:
                    os.system('rm -r hillshade_newest_headwall_line_*')
            break

        # if allow grid has been submit, then marked as done, we don't check results for each grids here
        update_subset_info(sub_txt, key_list=['proc_status','proc_done_time'], info_list=['done',str(datetime.now())])
        # remove no need dem files
        remove_no_need_dem_files(b_remove=b_remove_dem)

    return True


def sync_log_files(process_node,r_log_dir,process_log_dir):
    # copy complete id list, dem info, grid_no_dem_ids.txt to remote machine
    files_to_processNode = ['strip_dem_cover_grids.txt','tile_dem_cover_grids.txt','grid_complete_ids.txt','grid_no_dem_ids.txt']
    for file in files_to_processNode:
        scp_communicate.copy_file_folder_to_remote_machine(process_node, os.path.join(r_log_dir,file),os.path.join(process_log_dir, file))

    files_from_processNode = ['grid_dem_diff_less2dem_ids.txt','grid_no_valid_dem_ids.txt','grid_no_headwall_ids.txt',
                              'grid_no_subscidence_poly_ids.txt','grid_no_watermask_ids.txt']

    remote_name = process_node[1:].replace('_host', '')  # change $curc_host to curc
    for file in files_from_processNode:
        # copy the file, do not overwrite the local file
        remote_file = os.path.join(process_log_dir, io_function.get_name_by_adding_tail(file,remote_name))
        scp_communicate.copy_file_folder_from_remote_machine(process_node, os.path.join(r_log_dir,file),remote_file)
        # if they are new ids, then merged to "file"
        local_file = os.path.join(process_log_dir, file)
        remote_ids = io_function.read_list_from_txt(remote_file) if os.path.isfile(remote_file) else []    # no need, to int
        local_ids = io_function.read_list_from_txt(local_file) if os.path.isfile(local_file) else []
        new_ids = [id for id in remote_ids if id not in local_ids ]
        if len(new_ids) < 1:
            continue
        else:
            local_ids.extend(new_ids)
            io_function.save_list_to_txt(local_file,local_ids)


def make_note_all_task_done(extent_shp, reomte_node):
    if os.path.isdir(grid_ids_txt_dir) is False:
        io_function.mkdir(grid_ids_txt_dir)

    shp_grid_id_txt, log_grid_ids_txt, log_grid_ids_txt_done= get_extent_grid_id_txt_done_files(extent_shp)

     # shp_grid_id_txt should be in the current folder
    if os.path.isfile(log_grid_ids_txt) is False:
        io_function.copy_file_to_dst(shp_grid_id_txt,log_grid_ids_txt)

    if os.path.isfile(log_grid_ids_txt_done) is False:
        io_function.save_list_to_txt(log_grid_ids_txt_done,['Done'])
        # copy the curc
        r_grid_ids_txt_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir/grid_ids_txt'
        scp_communicate.copy_file_folder_to_remote_machine(reomte_node,r_grid_ids_txt_dir,log_grid_ids_txt_done)

def main(options, args):
    extent_shp = args[0]
    task_list = [args[item] for item in range(1, len(args)) ]
    # task_name = args[1]
    if len(task_list) < 1:
        raise ValueError('There is no task: %s'%str(task_list))

    # local_grid_id_txt is in the current dir
    # log_grid_ids_txt, log_grid_ids_txt_done is in grid_ids_txt_dir
    local_grid_id_txt, log_grid_ids_txt, log_grid_ids_txt_done = get_extent_grid_id_txt_done_files(extent_shp)
    # check if it has been complete
    if os.path.isfile(log_grid_ids_txt_done):
        basic.outputlogMessage('Tasks for extent %s have been completed'%extent_shp)
        return True

    r_working_dir = '/scratch/summit/lihu9680/Arctic/dem_processing' if options.remote_working_dir is None else options.remote_working_dir
    r_log_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir/log_dir' if options.remote_log_dir is None else options.remote_log_dir
    process_node = '$curc_host' if options.process_node is None else options.process_node
    download_node = '$curc_host' if options.download_node is None else options.download_node

    max_grid_count = options.max_grids
    b_remove_tmp_folders = options.b_remove_tmp_folders
    b_dont_remove_DEM_files = options.b_dont_remove_DEM_files
    b_no_slurm = options.b_no_slurm
    b_divide_to_subsets = True
    b_main_preProc = options.b_main_preProc

    # modify the folder name of subsets
    global subset_shp_dir
    subset_shp_dir = subset_shp_dir + '_' +io_function.get_name_no_ext(extent_shp)
    subset_shp_dir = os.path.join(subset_message_dir,subset_shp_dir)
    global msg_file_pre
    msg_file_pre = io_function.get_name_no_ext(extent_shp) + '_' + msg_file_pre

    grid_ids_to_process_txt = os.path.join(subset_message_dir,io_function.get_name_no_ext(extent_shp) +'_' + 'grid_ids_to_process.txt')

    # build map dem cover grid (take time, but only need to run once at the beginning)
    build_dict_of_dem_cover_grid_ids(dem_strip_shp, grid_20_shp, strip_dem_cover_grids_txt)
    build_dict_of_dem_cover_grid_ids(dem_tile_shp, grid_20_shp, tile_dem_cover_grids_txt)

    # get grids for processing
    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    gird_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)

    # get grid ids based on input extent
    grid_polys, grid_ids = get_grid_20(extent_shp,all_grid_polys, all_ids)

    # based on extent shape, subset grid_20_id_raster
    # # using gdalwarp to crop the mask, also have 0.5 pixel offset, so dont use it
    # grid_20_id_raster_sub = io_function.get_name_by_adding_tail(os.path.basename(grid_20_id_raster),'sub')
    # if RSImageProcess.subset_image_by_shapefile(grid_20_id_raster,extent_shp,save_path=grid_20_id_raster_sub) is False:
    #     return False

    # read grid_ids_2d, then mask it
    grid_ids_2d, grid_nodata = raster_io.read_raster_one_band_np(grid_20_id_raster)  # 2d array of gird ids
    # rasterize grid_polys, will served as mask.
    grid_ids_2d_mask = raster_io.burn_polygons_to_a_raster(grid_20_id_raster,grid_polys,1,None)
    # raster_io.save_numpy_array_to_rasterfile(grid_ids_2d_mask,'grid_ids_2d_mask.tif',grid_20_id_raster,nodata=255)  # save to disk for checking
    loc_masked_out = np.where(grid_ids_2d_mask != 1)
    # grid_ids_2d[ loc_masked_out ] = grid_nodata
    visit_np = np.zeros_like(grid_ids_2d, dtype=np.uint8)
    visit_np[loc_masked_out] = 1    # 1 indicate already visited
    visit_np[np.where(grid_ids_2d == grid_nodata)] = 1    # 1 indicate already visited

    subset_id = -1
    # on tesia, uist, vpn-connected laptop
    if b_main_preProc:
        io_function.mkdir(subset_shp_dir)
        sync_log_files(process_node, r_log_dir, process_log_dir)
        update_complete_grid_list(grid_ids, task_list)

    # divide a region into many subsets
    if b_main_preProc:
        # remove grids that has been complete or ignored
        ignore_ids = get_complete_ignore_grid_ids()
        num_grid_ids = save_grid_ids_need_to_process(grid_ids, ignore_ids=ignore_ids,
                                                     save_path=grid_ids_to_process_txt)
        if num_grid_ids < 1:
            make_note_all_task_done(extent_shp, process_node)
        else:
            # divide a region into many subsets
            while True:
                subset_id += 1
                # if the input is not a shapefile, then don't divide it to many subsets
                if extent_shp.endswith('.txt'):
                    select_grid_polys, selected_gird_ids = grid_polys, grid_ids
                    if len(selected_gird_ids) > 2000:
                        raise ValueError('There are too many grid to process once')
                    b_divide_to_subsets = False
                    subset_id = 999999
                    select_grids_shp = os.path.join(subset_shp_dir, io_function.get_name_no_ext(extent_shp) + '_sub%d' % subset_id + '.shp')
                    save_selected_girds_and_ids(selected_gird_ids, select_grid_polys, gird_prj, select_grids_shp)
                    break

                else:
                    select_grids_shp = os.path.join(subset_shp_dir, io_function.get_name_no_ext(extent_shp) + '_sub%d' % subset_id + '.shp')
                    # when re-run this, each subset will be the same or some grids in the subset would be removed if they has been completed (or ignored)
                    select_grid_polys, selected_gird_ids = get_grids_for_download_process(grid_polys, grid_ids, ignore_ids, max_grid_count,
                                                                                          grid_ids_2d, visit_np,select_grids_shp,proj=gird_prj)
                if selected_gird_ids is None:
                    break  # no more grids
                if len(selected_gird_ids) < 1:
                    continue

                subset_info_txt = os.path.join(subset_message_dir, msg_file_pre + str(subset_id).zfill(6)+'.txt')
                if os.path.isfile(subset_info_txt) is False:
                    # init the file
                    update_subset_info(subset_info_txt,
                                       key_list=['id', 'createTime', 'shp', 'pre_status', 'pre_node','proc_status'],
                                       info_list=[subset_id, str(datetime.now()), os.path.basename(select_grids_shp),
                                                  'notYet', 'unknown' ,'notYet'])

    b_preProc_complete = False
    while True:
        # on tesia, uist, vpn-connected laptop
        if machine_name == 'ubuntu' or machine_name == 'uist-int-colorado-edu' or 'colorado.edu' in machine_name or 'MacBook' in machine_name:

            # if subset_id for download is far more ahead than processing (curc), then wait, in case occupy too much storage
            while True and b_no_slurm is False:
                remote_sub_txt_list = get_subset_info_txt_list('proc_status', ['notYet', 'working'], remote_node=process_node,remote_folder=r_working_dir)
                if len(remote_sub_txt_list) > download_ahead_proc:
                    print(datetime.now(), 'there is %d subset have not complete,'
                                          ' wait 300 seconds, avoid downloading too many data' % len(remote_sub_txt_list))
                    time.sleep(300)
                else:
                    break

            subset_txt_list = get_subset_info_txt_list('pre_status', ['notYet', 'working'], local_folder=subset_message_dir)
            if len(subset_txt_list) > 0:
                subset_info_txt = None
                for txt in subset_txt_list:
                    info_dict = get_subset_info(txt)
                    if info_dict['pre_status'] == 'working' and info_dict['pre_node'] == machine_name:
                        subset_info_txt = txt
                        break
                    if info_dict['pre_status'] == 'notYet':
                        subset_info_txt = txt
                        break

                if subset_info_txt is not None:
                    # download and unpack ArcticDEM, do registration, send to curc
                    if download_process_send_arctic_dem(subset_info_txt, r_working_dir,process_node,task_list,
                                                        b_send_data = b_no_slurm==False) is True:
                        continue
                else:
                    b_preProc_complete = True
            else:
                subset_txt_all_list = get_subset_info_txt_list('pre_status', ['notYet', 'working','done'],local_folder=subset_message_dir)
                if len(subset_txt_all_list)==0:
                    # wait 10 minutes, until main computer create new subset info text files
                    print(datetime.now(), 'No subset info txt, wait 10 minutes')
                    time.sleep(600)
                else:
                    b_preProc_complete = True

            # copy file from remote machine
            if b_no_slurm is False and b_main_preProc:
                copy_results_from_remote_node()

                sync_log_files(process_node, r_log_dir, process_log_dir)

                # update complete id list
                update_complete_grid_list(grid_ids, task_list)

            # save this to disk, to check progress, if there are not too many grids (<100),
            # we can use this one to process withtou divide grids to many subsets
            if b_main_preProc:
                num_grid_ids = save_grid_ids_need_to_process(grid_ids,save_path=grid_ids_to_process_txt)
                if num_grid_ids < 1:
                    make_note_all_task_done(extent_shp,process_node)

            if b_no_slurm:
                # process ArcticDEM using local computing resource
                if produce_dem_products(task_list, b_remove_job_folder=b_remove_tmp_folders,b_remove_dem=b_dont_remove_DEM_files,no_slurm=b_no_slurm, message_dir=subset_message_dir) is False:
                    break

            if b_divide_to_subsets is False or b_preProc_complete is True:
                break

        elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:  # curc
            # process ArcticDEM using the computing resource on CURC
            if produce_dem_products(task_list,b_remove_job_folder=b_remove_tmp_folders,b_remove_dem=b_dont_remove_DEM_files, message_dir=subset_message_dir) is False:
                break

        else:
            print('unknown machine : %s '%machine_name)
            break

        if b_main_preProc:
            # remove no need dem files
            remove_no_need_dem_files(b_remove=b_dont_remove_DEM_files)

    # monitor results in remote computer
    check_time = 200
    while check_time > 0 and b_no_slurm==False and b_main_preProc is True:
        # on tesia, uist, vpn-connected laptop
        if machine_name == 'ubuntu' or machine_name == 'uist-int-colorado-edu' or 'colorado.edu' in machine_name or 'MacBook' in machine_name:
            print(datetime.now(),'wait 10 min for results in computing nodes')
            time.sleep(600)
            # copy file from remote machine
            copy_results_from_remote_node()
            # sync complete id list, dem info, no dem grids etcs.
            sync_log_files(process_node, r_log_dir, process_log_dir)
            # update complete id list
            update_complete_grid_list(grid_ids, task_list)
            # remove no need dem files
            remove_no_need_dem_files(b_remove=b_dont_remove_DEM_files)
            remote_sub_txt = get_subset_info_txt_list('proc_status', ['notYet', 'working'], remote_node=process_node,
                                                      remote_folder=r_working_dir)
            if len(remote_sub_txt) < 1 and check_time != 1:
                check_time = 1  # set to 1, then will only check one more time
            else:
                check_time -= 1
        else:
            break





if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp task_name ... (dem_diff, segment, dem_headwall)"
    parser = OptionParser(usage=usage, version="1.0 2021-10-11")
    parser.description = 'Introduction: process ArcticDEM with limited storage  '

    parser.add_option("-n", "--max_grids",
                      action="store", dest="max_grids", type=int, default=200,
                      help="the maximum number of grids for process in parallel, large storage allows large number")

    parser.add_option("-w", "--remote_working_dir",
                      action="store", dest="remote_working_dir",
                      help="the working directory in remote machine.")

    parser.add_option("-l", "--remote_log_dir",
                      action="store", dest="remote_log_dir",
                      help="the log directory in remote machine.")

    parser.add_option("-p", "--process_node",
                      action="store", dest="process_node",
                      help="the username and machine for processing ")

    parser.add_option("-d", "--download_node",
                      action="store", dest="download_node",
                      help="the username and machine for download ")

    parser.add_option("", "--b_dont_remove_tmp_folders",
                      action="store_false", dest="b_remove_tmp_folders",default=True,
                      help="if set, then dont remove processing folders of each job")

    parser.add_option("", "--b_dont_remove_DEM_files",
                      action="store_false", dest="b_dont_remove_DEM_files",default=True,
                      help="if set, then dont ArcticDEM (strip and mosaic) that have been processed")

    parser.add_option("", "--b_no_slurm",
                      action="store_true", dest="b_no_slurm",default=False,
                      help="if set, dont submit a slurm job, run job using local machine ")

    parser.add_option("", "--b_main_preProc",
                      action="store_true", dest="b_main_preProc",default=False,
                      help="if set, this computer is the main computer for pre-processing, "
                           "only the main computer would divide a region into subsets, check completeness,"
                           "the main computer must be run")



    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

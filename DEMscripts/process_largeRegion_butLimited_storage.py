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
from dem_common import process_log_dir, grid_complete_list_txt, strip_dem_cover_grids_txt, tile_dem_cover_grids_txt
if os.path.isdir(process_log_dir) is False:
    io_function.mkdir(process_log_dir)

from dem_common import arcticDEM_tile_tarball_dir,arcticDEM_tile_reg_tif_dir,tarball_dir,arcticDEM_reg_tif_dir,ArcticDEM_tmp_dir

# results dir
from dem_common import grid_hillshade_newest_HDLine_dir, grid_dem_diffs_dir

from produce_DEM_diff_ArcticDEM import get_grid_20

# scripts
dem_download_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/download_arcticDEM.py')
dem_unpack_reg_py = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/DEMscripts/ArcticDEM_unpack_registration.py')
subset_shp_dir = 'subset_grids_shp'



def update_subset_info(txt_path, key_list=None, info_list=None):
    # maintain a info of subset for processing, dict
    # id: subset_id
    # shp: the shapefile contain all grids in this subset
    # "pre_status": the status of downloading and registration of ArcticDEM, has values: 'notYet', 'working', 'done'
    # 'proc_status': the status of processing ArcticDEM, has values of 'notYet', 'working', 'done'

    info_dict = {}
    if os.path.isfile(txt_path):
        info_dict = io_function.read_dict_from_txt_json(txt_path)
    for key, info in zip(key_list, info_list):
        info_dict[key] = info
    io_function.save_dict_to_txt_json(txt_path,info_dict)

def get_subset_info(txt_path):
    return io_function.read_dict_from_txt_json(txt_path)

def get_subset_info_txt_list(key,values,remote_node=None, remote_folder=None, local_folder='./'):
    # get subset info with specific key and values
    # 'subset%d.txt' % subset_id
    # sync 'subset%d.txt' files from remote node
    if remote_node is not None:
        # this will overwrite 'subset%d.txt' in local
        scp_communicate.copy_file_folder_from_remote_machine(remote_node,os.path.join(remote_folder,'subset*.txt'),
                                                             os.path.join(local_folder,'.'))

    select_txt_list = []
    subinfo_txt_list = io_function.get_file_list_by_pattern(local_folder,'subset*.txt')
    for txt in subinfo_txt_list:
        info_dict = get_subset_info(txt)
        if info_dict[key] in values:
            select_txt_list.append(txt)
    return select_txt_list


def download_process_send_arctic_dem(subset_info_txt, r_working_dir, remote_node):
    # this function run on tesia

    subset_info = get_subset_info(subset_info_txt)

    # if subset_id for download is far more ahead than processing (curc), then wait, in case occupy too much storage
    while True:
        remote_sub_txt = get_subset_info_txt_list('proc_status',['notYet', 'working'],remote_node=remote_node,remote_folder=r_working_dir)
        if len(remote_sub_txt) > 5:
            print(datetime.now(),'there is %d subset in remote machine have not complete,'
                                 ' wait 300 seconds to download too many data'%len(remote_sub_txt))
            time.sleep(300)
        else:
            break

    update_subset_info(subset_info_txt, key_list=['pre_status'], info_list=['working'])

    # download strip tarball, upack, registration
    res = os.system(dem_download_py + ' ' + subset_info['shp'] + ' ' + dem_strip_shp )
    if res != 0:
        sys.exit(1)
    res = os.system(dem_unpack_reg_py + ' ' + tarball_dir  + ' -d ' + arcticDEM_reg_tif_dir + ' -r ' )  # -r: --remove_inter_data
    if res != 0:
        sys.exit(1)

    # download tile tarball and unpack
    res = os.system(dem_download_py + ' ' + subset_info['shp'] + ' ' + dem_tile_shp )
    if res != 0:
        sys.exit(1)
    res = os.system(dem_unpack_reg_py + ' ' +arcticDEM_tile_tarball_dir + ' -d ' + arcticDEM_tile_reg_tif_dir + ' -r ' )
    if res != 0:
        sys.exit(1)

    # send to remote machine
    rsync_sh = os.path.join(ArcticDEM_tmp_dir,'rsync_to_curc.sh')
    res = os.system(rsync_sh)
    if res != 0:
        sys.exit(1)

    update_subset_info(subset_info,key_list=['pre_status'],info_list=['done'])
    # copy to remote machine
    if scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir, subset_info_txt):
        scp_communicate.copy_file_folder_to_remote_machine(remote_node, r_working_dir,subset_shp_dir)


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


def b_exist_dem_hillshade_newest_HWLine_grid(id):
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
    dem_diff_files = io_function.get_file_list_by_pattern(grid_dem_diffs_dir, '*_DEM_diff_grid%d.tif' % id)
    if len(dem_diff_files) == 1:
        return True
    elif len(dem_diff_files) > 1:
        basic.outputlogMessage('warning, There are multiple hillshade (newest) HDLine tif for grid: %d' % id)
        for item in dem_diff_files:
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
        # we may check more task results: segment, dem_headwall_grid, dem_headwall

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


def get_grids_for_download_process(grid_polys, grid_ids, max_grid_count, grid_ids_2d, visit_np, save_path, proj=None):
    '''
    get grids for donwload ArcticDEM and ids
    :param grid_polys:
    :param grid_ids:
    :param max_grid_count: max grid count, it will return number close to this value
    :param grid_ids_2d: gird_ids in 2d array
    :param visit_np: visit_np like a mask, to indicate which pixels has been checked, 0 no visit previously, 1 visited
    :param proj: projection information, for save grid_poygons
    :return:
    '''

    # find a connected region with for donwload and processing, and save to files
    seed_loc = np.where(grid_ids_2d == grid_ids[0])
    selected_gird_id_list = [grid_ids[0]]
    seed_list = [ [seed_loc[0][0], seed_loc[1][0]]]
    while len(selected_gird_id_list) < max_grid_count and len(seed_list) > 0:
        # find neighbours
        new_seeds = find_neighbours_2d(grid_ids_2d,visit_np,seed_list[0],8)
        del seed_list[0]
        seed_list.extend(new_seeds)
        # find new ids
        for seed in new_seeds:
            row, col = seed
            selected_gird_id_list.append( grid_ids_2d[row, col])

    select_grid_polys = [ grid_polys[grid_ids.index(item) ] for item in selected_gird_id_list ]

    # save to shapefile to download and processing
    save_pd = pd.DataFrame({'id':selected_gird_id_list, 'Polygon':select_grid_polys})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',proj,save_path)

    return select_grid_polys, selected_gird_id_list


def remove_no_need_dem_files():
    pass


def main(options, args):
    extent_shp = args[0]
    task_list = [args[item] for item in range(1, len(args)) ]
    # task_name = args[1]
    if len(task_list) < 1:
        raise ValueError('There is no task: %s'%str(task_list))

    r_working_dir = '/home/lihu9680/Data/dem_processing' if options.remote_working_dir is None else options.remote_working_dir
    process_node = '$curc_host' if options.process_node is None else options.process_node
    download_node = '$curc_host' if options.download_node is None else options.download_node

    max_grid_count = options.max_grids


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

    subset_id = 0
    io_function.mkdir(subset_shp_dir)

    # grid is being process but not complete yet
    running_grid_ids = []

    while True:
        # on tesia, uist, vpn-connected laptop
        if machine_name == 'ubuntu' or machine_name == 'uist' or 'colorado.edu' in machine_name or 'MacBook' in machine_name:

            no_complete_polys, no_complete_ids = get_not_completed_grids(grid_polys, grid_ids)
            if len(no_complete_ids) < 1:
                break
            # update with running girds
            to_work_ids = [item for item in no_complete_ids if item not in running_grid_ids]
            if len(to_work_ids) < 1:
                print('No new grids need to download or process')
                break
            to_polys = [no_complete_polys[no_complete_ids.index(item)] for item in to_work_ids]

            select_grids_shp = os.path.join(subset_shp_dir,io_function.get_name_by_adding_tail(os.path.basename(grid_20_shp),'sub%d' % subset_id))
            select_grid_polys, selected_gird_ids = get_grids_for_download_process(to_polys, to_work_ids, max_grid_count,
                                                                                  grid_ids_2d, visit_np,
                                                                                  select_grids_shp, proj=gird_prj)
            # print(len(select_grid_polys),len(selected_gird_ids),selected_gird_ids)
            running_grid_ids.extend(selected_gird_ids)
            subset_info_txt = 'subset%d.txt'%subset_id
            update_subset_info(subset_info_txt,key_list=['id','shp','pre_status','proc_status'],
                               info_list=[subset_id,select_grids_shp, 'notYet','notYet'])

            # download and unpack ArcticDEM, do registration, send to curc
            download_process_send_arctic_dem(subset_info_txt, r_working_dir,process_node)

        elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:  # curc
            # process ArcticDEM using the computing resource on CURC
            pass
        else:
            print('unknown machine name')
            break

        # update complete id list
        update_complete_grid_list(grid_ids, task_list)

        # remove no need dem files

        subset_id += 1





if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp task_name ... (dem_diff, segment, dem_headwall)"
    parser = OptionParser(usage=usage, version="1.0 2021-10-11")
    parser.description = 'Introduction: process ArcticDEM with limited storage  '

    parser.add_option("-n", "--max_grids",
                      action="store", dest="max_grids", type=int, default=100,
                      help="the maximum number of grids for process in parallel, large storage allows large number")

    parser.add_option("-w", "--remote_working_dir",
                      action="store", dest="remote_working_dir",
                      help="the working directory in remote machine.")

    parser.add_option("-p", "--process_node",
                      action="store", dest="process_node",
                      help="the username and machine for processing ")

    parser.add_option("-d", "--download_node",
                      action="store", dest="download_node",
                      help="the username and machine for download ")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

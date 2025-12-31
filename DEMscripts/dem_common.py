#!/usr/bin/env python
# Filename: dem_common.py 
"""
introduction: put some variables

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 March, 2021
"""

import os,sys
import re

machine_name = os.uname()[1]
import time
from datetime import datetime
import numpy as np

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import parameters

# some folder paths
if machine_name == 'uist-int-colorado-edu':
    ArcticDEM_tmp_dir = '/home/lhuang/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
    data_dir = '/home/lhuang/Data'

elif machine_name == 'donostia.int.colorado.edu':  # donostia
    ArcticDEM_tmp_dir = '/home/lhuang/Bhaltos2/ArcticDEM_tmp_dir'
    data_dir = '/home/lhuang/Data'

elif machine_name == 'ubuntu':  # tesia
    ArcticDEM_tmp_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
    data_dir = '/home/lihu9680/Data'

elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:   # curc
    ArcticDEM_tmp_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir'
    data_dir = '/projects/lihu9680/Data'
else:
    ArcticDEM_tmp_dir = './'
    data_dir = './'

download_ahead_proc = 3

setting_ini = 'setting.ini'
if os.path.isfile(setting_ini):
    ArcticDEM_tmp_dir = parameters.get_directory(setting_ini,'ArcticDEM_tmp_dir')
    data_dir = parameters.get_directory(setting_ini,'data_dir')
    tmp = parameters.get_digit_parameters_None_if_absence(setting_ini,'download_ahead_proc', 'int')
    if tmp is not None:
        download_ahead_proc = tmp


# tarball_dir = os.path.join(ArcticDEM_tmp_dir,'tarballs')    # strip version of ArcticDEM
tarball_dir = 'tarballs'    # strip version of ArcticDEM

arcticDEM_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir,'registration_tifs')
relative_dem_dir = os.path.join(ArcticDEM_tmp_dir,'dem_relative_8bit')

grid_dem_diffs_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs')
grid_dem_diffs_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_8bit')
grid_dem_diffs_color_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_color')
grid_dem_diffs_segment_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_segment_results')

grid_matchtag_sum_dir = os.path.join(ArcticDEM_tmp_dir,'grid_matchtag_sum_tifs')

dem_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'dem_slope_8bit')
dem_slope_dir = os.path.join(ArcticDEM_tmp_dir,'dem_slope')
dem_hillshade_dir = os.path.join(ArcticDEM_tmp_dir,'dem_hillshade')
dem_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_tpi_8bit')

grd_hillshade_newest_on_top_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_newest_top_grid')

dem_headwall_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_headwall_shp')
grid_dem_headwall_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_headwall_shp_grid')
grid_hillshade_newest_HDLine_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_newest_HWLine_grid')


dem_hillshade_subImages_headwall = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_subImages_headwall')


grid_dem_subsidence_select = os.path.join(ArcticDEM_tmp_dir,'grid_dem_subsidence_select')

# the mosaic version of AricticDEM
# arcticDEM_tile_tarball_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_tarballs')
arcticDEM_tile_tarball_dir = 'arcticdem_mosaic_tarballs'
arcticDEM_tile_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_reg_tifs')
arcticDEM_tile_hillshade_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_hillshade')
arcticDEM_tile_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope_8bit')
arcticDEM_tile_slope_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope')
arcticDEM_tile_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_tpi_8bit')
# dem_pattern = '*reg_dem.tif'


# surface water mask
mask_water_dir = os.path.join(data_dir, 'global_surface_water' , 'extent_epsg3413')

grid_20_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp')
grid_20_id_raster = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km_id.tif')
dem_strip_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7/ArcticDEM_Strip_Index_Rel7.shp')
dem_tile_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Tile_Index_Rel7/ArcticDEM_Tile_Index_Rel7.shp')

if os.path.isfile(setting_ini):
    dem_strip_shp = parameters.get_file_path_parameters(setting_ini,'dem_strip_shp')
    dem_tile_shp = parameters.get_file_path_parameters(setting_ini,'dem_tile_shp')


# some log and information files
subset_message_dir = os.path.join(ArcticDEM_tmp_dir, 'subset_message_dir')
if 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:   # curc
    subset_message_dir = './'
process_log_dir = os.path.join(ArcticDEM_tmp_dir, 'log_dir')
grid_complete_list_txt = os.path.join(process_log_dir,'grid_complete_ids.txt')  # store ids of grids that have completed
# manually exclude some grids that dont have enough data
grid_excluded_list_txt = os.path.join(process_log_dir,'grid_exclude_ids.txt')   # store ids of grids that manually exclude

grid_dem_diff_less2dem_txt = os.path.join(process_log_dir,'grid_dem_diff_less2dem_ids.txt')  # store ids of grids that has less than 2 DEM (cannot calculate DEM differnce)
grid_no_dem_txt = os.path.join(process_log_dir,'grid_no_dem_ids.txt')  # store ids of grids that don't have DEM (strip version) overlap

# store ids of grids that have overlap of DEM  (strip version), but the coverage is too smaller or all overlap DEM are invalid.
grid_no_valid_dem_txt = os.path.join(process_log_dir,'grid_no_valid_dem_ids.txt')

# some grid in very high latitude dont overlap with the global surface water, record these ids
grid_no_water_mask_txt = os.path.join(process_log_dir,'grid_no_watermask_ids.txt')

# some place that is really flat, cannot detect headwall based on slope from it
grid_no_headwall_txt = os.path.join(process_log_dir,'grid_no_headwall_ids.txt')

# based on the criteria, there is no subsidence polygons
grid_no_subscidence_poly_txt = os.path.join(process_log_dir,'grid_no_subscidence_poly_ids.txt')

#based on the criteria, there is no headwall lines after filtering
grid_no_rippleSel_headwall_line_txt = os.path.join(process_log_dir,'grid_no_rippleSel_headwall_line_ids.txt')

strip_dem_cover_grids_txt = os.path.join(process_log_dir,'strip_dem_cover_grids.txt') # each strip cover how many grids (ids), dict
tile_dem_cover_grids_txt = os.path.join(process_log_dir,'tile_dem_cover_grids.txt') # each tile cover how many grids (ids), dict

# store grid id txt files for each extent, and indicator of them have been completed.
grid_ids_txt_dir = os.path.join(ArcticDEM_tmp_dir, 'grid_ids_txt')

# rts results
grid_rts_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'grid_rts_shp')



#########################
## some common function
def get_corresponding_grid_ids_txt(extent_shp):
    # if it's in the same folder of the extent_shp
    grid_ids_txt = os.path.splitext(extent_shp)[0] + '_grid_ids.txt'
    if os.path.isfile(grid_ids_txt):
        print('grid_ids.txt for %s is in the same directory'%extent_shp)
        return grid_ids_txt

    # otherwise, just return a file name in current folder (don't have to be existed)
    grid_ids_txt = os.path.splitext(os.path.basename(extent_shp))[0] +'_grid_ids.txt'
    return grid_ids_txt

def get_extent_grid_id_txt_done_files(extent_shp):
    grid_ids_txt = get_corresponding_grid_ids_txt(extent_shp)
    log_grid_ids_txt = os.path.join(grid_ids_txt_dir,os.path.basename(grid_ids_txt))
    log_grid_ids_txt_done = log_grid_ids_txt + '_done'
    return grid_ids_txt, log_grid_ids_txt, log_grid_ids_txt_done

def get_grid_id_from_path(item):
    return int(re.findall('grid\d+', os.path.basename(item))[0][4:])

def save_id_grid_no_result(grid_id,file_path):
    if os.path.isdir(process_log_dir) is False:
        io_function.mkdir(process_log_dir)

    id_list = []
    if os.path.isfile(file_path):
        id_list = io_function.read_list_from_txt(file_path)  # no need covert to int
    id_str = str(grid_id)
    if id_str in id_list:
        return True
    else:
        # save by adding one line
        with open(file_path, 'a') as f_obj:
            f_obj.writelines(str(grid_id) + '\n')
        return True


def find_neighbours_grids(grid_ids_2d,grid_id,connect):
    '''
    find neighbourhood grids, a similar funciotn is "find_neighbours_2d"
    :param grid_ids_2d: 2D data, read from "ArcticDEM_grid_20km_id.tif"
    :param seed: a seed
    :param connect: pixel connectivity, 4 or 8
    :return: a list containing neighbours grids,
    '''

    height, width = grid_ids_2d.shape
    rows, cols  = np.where(grid_ids_2d == grid_id)
    # print(rows, cols)
    if len(rows) != 1:
        raise ValueError(f'There are {len(rows)} grids with the ID: {grid_id} in the numpy array, should be only one grid')

    y,x = rows[0],cols[0]

    neigh_range = [-1,0,1]
    neighbours = [[i,j] for i in neigh_range for j in neigh_range  ]
    neighbours.remove([0,0])

    # for grid ids, +- 1 or 191 may also got the id for their neighbour

    # distance within 1
    if connect==4:
        connectivity = [ [y+dy, x+dx] for (dy,dx) in neighbours if (dy*dy + dx*dx) <= 1 ]
    # distance within sqrt(2)
    elif connect==8:
        connectivity = [[y + dy, x + dx] for (dy, dx) in neighbours if (dy * dy + dx * dx) <= 2]
    else:
        raise ValueError('Only accept connectivity of 4 or 8')

    # new_seeds = []
    neighour_grid_ids = []
    for [y,x] in connectivity:
        # out extent
        if y<0 or x<0 or y >= height or x >= width:
            continue

        # new_seeds.append([y,x])
        neighour_grid_ids.append(grid_ids_2d[y,x])

    return neighour_grid_ids



def check_create_lock(lock_path, message):
    # lock_path: the absolute path for a lock file
    # message: information for why it's locked
    check_lock_time = 0
    while os.path.isfile(lock_path):
        print(datetime.now(),'checked %d times: wait 300 seconds'%check_lock_time,message)
        time.sleep(300)
        check_lock_time += 1
    # create a lock
    with open(lock_path, 'w') as f_obj:
        f_obj.writelines('locked at ' + str(datetime.now()) +' by %s:%d'%(machine_name, os.getpid()) + '\n')

def release_lock(lock_path):
    if os.path.exists(lock_path):
        os.remove(lock_path)
    else:
        print("Error: The lock file: %s does not exist"%lock_path)

if __name__ == '__main__':
    pass
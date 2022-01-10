#!/usr/bin/env python
# Filename: ArcticDEM_proc_grid.py
"""
introduction: Convert DEM difference to 8bit.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 20 March, 2021
"""

import os,sys
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic

# some folder paths

from dem_common import grid_dem_diffs_dir, grid_dem_diffs_8bit_dir

py8bit= os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py')

def dem_tif_to_8bit(input,output):
    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    dst_nodat = 255
    hist_max_percent=0.98
    hist_min_percent=0.02
    MIN_MAX_VALUE = '-500 500'     # five meters

    command_str = py8bit + ' ' + input + ' ' + output
    command_str += ' -n ' + str(dst_nodat)
    command_str += ' -u ' + str(hist_max_percent) + ' -l ' + str(hist_min_percent)
    command_str += ' -m ' + MIN_MAX_VALUE

    print(command_str)
    res = os.system(command_str)
    # basic.os_system_exit_code(command_str)
    if res == 0:
        return True
    else:
        return False

def test_dem_tif_to_8bit():
    dem_diff_list = io_function.get_file_list_by_pattern('./','*.tif')
    count = len(dem_diff_list)
    for idx, tif in enumerate(dem_diff_list):
        print('%d/%d convert %s to 8 bit'%(idx+1, count, tif))
        tif_8bit = io_function.get_name_by_adding_tail(tif, '8bit')
        output = os.path.join(grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
        dem_tif_to_8bit(tif,output)

def one_dem_diff_to_8bit(demDiff_tif):
    if os.path.isdir(grid_dem_diffs_8bit_dir) is False:
        io_function.mkdir(grid_dem_diffs_8bit_dir)
    tif_8bit = io_function.get_name_by_adding_tail(demDiff_tif, '8bit')
    output = os.path.join(grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
    if dem_tif_to_8bit(demDiff_tif, output) is False:
        basic.outputlogMessage('failed to generate 8bit grey image from DEM differnce')
        return False
    return True

def main():
    basic.setlogfile('log_convert_dem_diff_to8bit.txt')
    if os.path.isdir(grid_dem_diffs_8bit_dir) is False:
        io_function.mkdir(grid_dem_diffs_8bit_dir)

    dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir,'*DEM_diff_grid*.tif')
    count = len(dem_diff_list)
    failed_tifs = []
    for idx, tif in enumerate(dem_diff_list):
        print('%d/%d convert %s to 8 bit'%(idx+1, count, tif))
        tif_8bit = io_function.get_name_by_adding_tail(tif, '8bit')
        output = os.path.join(grid_dem_diffs_8bit_dir, os.path.basename(tif_8bit))
        if dem_tif_to_8bit(tif,output) is False:
            failed_tifs.append(tif)

    if len(failed_tifs)>0:
        io_function.save_list_to_txt('failed_dem_diff_to8bit.txt',failed_tifs)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# Filename: ArcticDEM_proc_grid.py
"""
introduction: Convert DEM difference to color relief

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 March, 2024
"""

import os,sys
from optparse import OptionParser
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic

# some folder paths

from dem_common import grid_dem_diffs_dir, grid_dem_diffs_color_dir

def dem_tif_to_colorReleif(input,output):
    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    color_text_file = 'dem_diff_color_5to5m.txt'

    command_str = 'gdaldem color-relief -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer ' \
                  + input + ' ' + ' %s '%color_text_file  + output

    print(command_str)
    res = os.system(command_str)
    # basic.os_system_exit_code(command_str)
    if res == 0:
        return True
    else:
        return False

# def test_dem_tif_to_colorReleif():
#     dem_diff_list = io_function.get_file_list_by_pattern('./','*.tif')
#     count = len(dem_diff_list)
#     for idx, tif in enumerate(dem_diff_list):
#         print('%d/%d convert %s to color relief '%(idx+1, count, tif))
#         tif_color = io_function.get_name_by_adding_tail(tif, 'color')
#         output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
#         dem_tif_to_colorReleif(tif,output)

def test_dem_tif_to_colorReleif_one():

    tif = os.path.expanduser('~/Data/dem_processing/grid_dem_diffs/grid_ids_DEM_diff_grid13965.tif')
    tif_color = io_function.get_name_by_adding_tail(tif, 'color')
    # output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
    dem_tif_to_colorReleif(tif,tif_color)

def one_dem_diff_to_colorRelief(demDiff_tif):
    if os.path.isdir(grid_dem_diffs_color_dir) is False:
        io_function.mkdir(grid_dem_diffs_color_dir)
    tif_color = io_function.get_name_by_adding_tail(demDiff_tif, 'color')
    output = os.path.join(grid_dem_diffs_color_dir, os.path.basename(tif_color))
    if dem_tif_to_colorReleif(demDiff_tif, output) is False:
        basic.outputlogMessage('failed to generate color relief from DEM difference')
        return False
    return True

def main(options, args):
    basic.setlogfile('log_convert_dem_diff_to_colorRelief.txt')

    if len(args) < 1:
        dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir, '*DEM_diff_grid*.tif')
        if os.path.isdir(grid_dem_diffs_color_dir) is False:
            io_function.mkdir(grid_dem_diffs_color_dir)
    else:
        dem_diff_file_or_dir = args[0]
        if os.path.isfile(dem_diff_file_or_dir):
            dem_diff_list = [dem_diff_file_or_dir]
        else:
            dem_diff_list = io_function.get_file_list_by_pattern(dem_diff_file_or_dir,'*DEM_diff_grid*.tif')
            if len(dem_diff_list) < 1:
                basic.outputlogMessage(f'No DEM diff files in {dem_diff_file_or_dir} for colorRelief')

    out_dir = grid_dem_diffs_color_dir if os.path.isdir(grid_dem_diffs_color_dir) else './'

    count = len(dem_diff_list)
    failed_tifs = []
    for idx, tif in enumerate(dem_diff_list):
        print('%d/%d convert %s to Color Relief'%(idx+1, count, tif))
        tif_color = io_function.get_name_by_adding_tail(tif, 'color')
        output = os.path.join(out_dir, os.path.basename(tif_color))
        if dem_tif_to_colorReleif(tif,output) is False:
            failed_tifs.append(tif)

    if len(failed_tifs)>0:
        io_function.save_list_to_txt('failed_dem_diff_to_color.txt',failed_tifs)


if __name__ == '__main__':
    usage = "usage: %prog [options] dem_diff or dem_diff_dir "
    parser = OptionParser(usage=usage, version="1.0 2024-3-21")
    parser.description = 'Introduction: producing DEM color relief '

    (options, args) = parser.parse_args()
    main(options, args)
    # test_dem_tif_to_colorReleif_one()

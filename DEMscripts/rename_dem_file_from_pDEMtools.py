#!/usr/bin/env python
# Filename: rename_dem_file_from_pDEMtools.py 
"""
introduction:

rename the DEM files that were saved by "download_arcticDEM_pDEMtools.py" before Feb 3, 2026,
making the filename be consistent with previous workflow


authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 04 February, 2026
"""

import os
import sys
from optparse import OptionParser

deeplabforRS = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.basic as basic
import basic_src.io_function as io_function

def rename_dem_files_in_one_grid_dir(grid_folder):
    if not os.path.isdir(grid_folder):
        basic.outputlogMessage(f'Warning, {grid_folder} is not a directory')
        return
    tif_list = io_function.get_file_list_by_ext('.tif', grid_folder, bsub_folder=False)
    json_list = io_function.get_file_list_by_ext('.json', grid_folder, bsub_folder=False)
    if len(tif_list) != len(json_list):
        raise ValueError('The count of tif and json files is different')

    # rename tif files
    tif_count = 0
    for tif in tif_list:
        if tif.endswith('_dem.tif'):
            continue
        tif_new = io_function.get_name_by_adding_tail(tif, 'dem')
        io_function.move_file_to_dst(tif, tif_new, overwrite=False, b_verbose=False)
        tif_count += 1

    # rename json files
    json_count = 0
    for j_file in json_list:
        if j_file.endswith('_dem.json'):
            continue
        json_new = io_function.get_name_by_adding_tail(j_file, 'dem')
        io_function.move_file_to_dst(j_file, json_new, overwrite=False, b_verbose=False)
        json_count += 1

    basic.outputlogMessage(f'Rename {tif_count} tif and {json_count} json files in {grid_folder}')

def main(options, args):
    if len(args) > 0:
        data_dir = args[0]
    else:
        data_dir = os.getcwd()
    basic.outputlogMessage(f'TO Re-name DEM file name under: {data_dir}')

    # e.g., dem_grid0005400357
    dem_grid_dir_list = io_function.get_file_list_by_pattern(data_dir, 'dem_grid??????????')
    basic.outputlogMessage(f'Found {len(dem_grid_dir_list)} grid folders')
    for idx, g_dir in enumerate(dem_grid_dir_list):
        rename_dem_files_in_one_grid_dir(g_dir)

    basic.outputlogMessage('Rename, Done')

if __name__ == '__main__':
    usage = "usage: %prog [options] data_dir "
    parser = OptionParser(usage=usage, version="1.0 2026-02-04")
    parser.description = 'Introduction: rename DEM file names'
    (options, args) = parser.parse_args()

    # Optionally, enforce required argument:
    # if len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)



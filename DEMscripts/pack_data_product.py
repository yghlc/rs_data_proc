#!/usr/bin/env python
# Filename: pack_data_product.py 
"""
introduction: organize data products derived from ArcticDEM, run on tesia: /tiampostorage/results_Lingcao/products_derived_from_ArcticDEM

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 24 October, 2022
"""

import os,sys
from datetime import datetime
import time

import tarfile

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic

from dem_common import get_grid_id_from_path

dir1="/BhaltosMount/Bhaltos/lingcaoHuang/ArcticDEM_results"
dir2="/Miavaig/Work/lingcaoHuang/ArcticDEM_results"


def readme_elevation_diff(grid_files):
    # Pixel-wise elevation differences with a spatial resolution of 2m, derived from ArctciDEM by:
    #
    # "the most recent elevation - the oldest elevation".
    #
    # grid_ids_DEM_diff_grid10741.tif: the elevation differences. The pixel value/ 100.0 is the difference (unit:  meters).
    #
    # grid_ids_date_diff_grid10741.txt : the indices and ArcticDEM files for producing the elevation differences, we can tell the acquisition dates from the file name.
    #
    # grid_ids_date_diff_grid10741_newIndex.tif: the index of the most recent ArcticDEM file at each pixel (grid_ids_date_diff_grid10741.txt)
    # grid_ids_date_diff_grid10741_oldIndex.tif: the index of the oldest ArcticDEM file at each pixel (grid_ids_date_diff_grid10741.txt)

    file_names = [ os.path.basename(item) for item in grid_files ]
    dem_diff_file = [item for item in file_names if 'DEM_diff' in item][0]
    file_names.remove(dem_diff_file)
    date_txt = [item for item in file_names if 'txt' in item][0]
    file_names.remove(date_txt)
    newIndex_file = [item for item in file_names if 'newIndex' in item][0]
    file_names.remove(newIndex_file)
    oldIndex_file = [item for item in file_names if 'oldIndex' in item][0]
    file_names.remove(oldIndex_file)
    date_diff = file_names[0]

    save_txt = os.path.abspath('readme.txt')

    with open(save_txt, 'w') as f_obj:
        f_obj.writelines('Pixel-wise elevation differences with a spatial resolution of 2m, derived from ArctciDEM by:\n')
        f_obj.writelines('"the most recent elevation - the oldest elevation".\n\n')

        f_obj.writelines('%s: the elevation differences. The pixel value/ 100.0 is the difference (unit:  meters).\n'%dem_diff_file)
        f_obj.writelines('%s: pixel-wise differences of acquisition dates between the most recent elevation and the oldest elevation.\n'%date_diff)
        f_obj.writelines('%s: the indices and ArcticDEM files for producing the elevation differences, we can tell the acquisition dates from the file names.\n'%date_txt)
        f_obj.writelines('%s: the index of the most recent ArcticDEM file at each pixel (%s).\n'%(newIndex_file,date_txt))
        f_obj.writelines('%s: the index of the oldest ArcticDEM file at each pixel (%s).\n'%(oldIndex_file,date_txt))

    return save_txt



def copy_pack_elevation_diff(ext_dir,ext_name):
    diff_dir = os.path.join(ext_dir,'grid_dem_diffs')
    diff_list = io_function.get_file_list_by_pattern(diff_dir,'*DEM_diff*.tif')
    basic.outputlogMessage('count of elevation difference: %d'%len(diff_list))

    save_dir = os.path.join('elevation-differences',ext_name)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    grid_ids = [get_grid_id_from_path(item) for item in diff_list ]
    for id in grid_ids:
        basic.outputlogMessage('packing data for grid %d'%id)
        save_tar = os.path.join(save_dir, 'dem_diffs_2m_grid%d.tar.gz' % id)
        if os.path.isfile(save_tar):
            basic.outputlogMessage('%s already exists, skip')
            continue

        grid_files = io_function.get_file_list_by_pattern(diff_dir,'*grid%d*'%id)
        # create a readme file
        readme_txt = readme_elevation_diff(grid_files)
        grid_files.append(readme_txt)

        # command_str = 'tar -czvf %s '%save_tar
        # for file in grid_files:
        #     command_str += ' %s'%file
        # basic.os_system_exit_code(command_str)

        with tarfile.open(save_tar, 'x:gz') as tar:
            for file in grid_files:
                tar.add(file, arcname=os.path.basename(file))



def main():

    for dir in [dir1, dir2]:
        basic.outputlogMessage('working on ArcticDEM_results, dir: %s'%dir)
        ext_folders = io_function.get_file_list_by_pattern(dir,'ext*')
        for ext_dir in ext_folders:
            basic.outputlogMessage('Working on ext folder, dir: %s' % ext_dir)
            ext_name =  os.path.basename(ext_dir).split('_')[0]
            basic.outputlogMessage('ext_name: %s' % ext_name)

            copy_pack_elevation_diff(ext_dir,ext_name)



if __name__ == '__main__':
    main()
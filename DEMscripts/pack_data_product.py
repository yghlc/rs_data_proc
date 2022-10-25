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
        basic.outputlogMessage('elevation-differences, packing data for grid %d'%id)
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

def readme_composited_image(grid_files):

    file_names = [os.path.basename(item) for item in grid_files]
    pixel_count_file = [item for item in file_names if '_count.tif' in item][0]
    file_names.remove(pixel_count_file)
    image_file = file_names[0]
    save_txt = os.path.abspath('readme.txt')
    with open(save_txt,'w') as f_obj:
        f_obj.writelines('composited imagery derived from ArcticDEM\n\n')
        f_obj.writelines('%s: a composited image, with lines of narrow-steep slopes on top and hillshades derived from the '
                         'most recent ArcticDEM in the background. '
                         'For the legend of lines in different color, please refer to '
                         'Fig. 2 at https://yghlc.github.io/validate-thaw-slump\n'%image_file)
        f_obj.writelines('%s: the number of lines of narrow-steep slopes at each pixel\n'%pixel_count_file)
    return save_txt

def copy_pack_composited_image(ext_dir,ext_name):
    hillshade_HWLine_dir = os.path.join(ext_dir,'dem_hillshade_newest_HWLine_grid')
    image_list = io_function.get_file_list_by_pattern(hillshade_HWLine_dir,'*.tif')
    image_list = [item for item in image_list if '_count.tif' not in os.path.basename(item)]    # remove *_count.tif
    basic.outputlogMessage('count of hillshade + HWLine image: %d' % len(image_list))

    save_dir = os.path.join('composited-images',ext_name)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    grid_ids = [get_grid_id_from_path(item) for item in image_list]
    for id in grid_ids:
        basic.outputlogMessage('composited-images, packing data for grid %d'%id)
        save_tar = os.path.join(save_dir, 'composited_image_2m_grid%d.tar.gz' % id)
        if os.path.isfile(save_tar):
            basic.outputlogMessage('%s already exists, skip')
            continue

        grid_files = io_function.get_file_list_by_pattern(hillshade_HWLine_dir,'*grid%d*'%id)
        # create a readme file
        readme_txt = readme_composited_image(grid_files)
        grid_files.append(readme_txt)

        with tarfile.open(save_tar, 'x:gz') as tar:
            for file in grid_files:
                tar.add(file, arcname=os.path.basename(file))

def readme_lines_slope_headwall(grid_files):
    file_names = [os.path.basename(item) for item in grid_files if item.endswith('.shp')]
    headwall_shp_file = [item for item in file_names if '_rippleSel.shp' in item][0]
    file_names.remove(headwall_shp_file)
    slope_shp_file = file_names[0]
    save_txt = os.path.abspath('readme.txt')
    with open(save_txt,'w') as f_obj:
        f_obj.writelines('Lines of narrow-steep slopes and potential headwalls of retrogressive thaw slumps\n\n')
        f_obj.writelines('%s: lines representing narrow-steep slopes\n'%slope_shp_file)
        f_obj.writelines('%s: lines reprenting potential headwalls of retrogressive thaw slumps\n'%headwall_shp_file)
    return save_txt

def copy_pack_lines_of_narrow_steep_slope(ext_dir,ext_name):
    lines_dir = os.path.join(ext_dir, 'dem_headwall_shp_grid')
    lines_shpDir_list = io_function.get_file_list_by_pattern(lines_dir, '*')
    basic.outputlogMessage('count of grid with lines_of_narrow_steep_slope_and_headwall: %d' % len(lines_shpDir_list))

    save_dir = os.path.join('lines-of-narrow-steep-slopes-and-headwall', ext_name)
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    grid_ids = [get_grid_id_from_path(item) for item in lines_shpDir_list]
    for id in grid_ids:
        basic.outputlogMessage('lines-of-narrow-steep-slopes-and-headwall, packing data for grid %d' % id)
        save_tar = os.path.join(save_dir, 'lines_of_narrow_steep_slope_and_headwall_grid%d.tar.gz' % id)
        if os.path.isfile(save_tar):
            basic.outputlogMessage('%s already exists, skip')
            continue

        grid_files = io_function.get_file_list_by_pattern(lines_dir, '*grid%d*/*' % id)
        # create a readme file
        readme_txt = readme_lines_slope_headwall(grid_files)
        grid_files.append(readme_txt)

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

            copy_pack_composited_image(ext_dir,ext_name)

            copy_pack_lines_of_narrow_steep_slope(ext_dir,ext_name)



if __name__ == '__main__':
    main()
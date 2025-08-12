#!/usr/bin/env python
# Filename: check_grid_id_in_area_ini.py 
"""
introduction: due to a bug:

        # this is a potential bug, for example: grid10, is in grid101, grid102, .....
        # if grid_str in basename:
        #     return shp

the grid shapefile, may not corresponding to the DEM diff grid in the area setting ini,
find these and remove them, allowing to re-run

run this in the working folder, such as: "ext09_for_ArcticDEM_proc"

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 12 August, 2025
"""
import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic
import vector_gpd

import parameters
from dem_common import get_grid_id_from_path

def check_one_area_ini_file(area_ini):
    all_polygons_labels = parameters.read_Parameters_file(area_ini,'all_polygons_labels')
    area_name = parameters.get_string_parameters(area_ini,'area_name')

    area_name_grid_id = get_grid_id_from_path(area_name)
    poly_grid_id = get_grid_id_from_path(all_polygons_labels)
    if area_name_grid_id == poly_grid_id:
        return True
    else:
        print(f'in {area_ini}, shp setting is wrong')
        return False


def remove_wrong_file_folder(area_ini):
    area_name = parameters.get_string_parameters(area_ini,'area_name')
    area_time = parameters.get_string_parameters(area_ini,'area_time')
    area_remark = parameters.get_string_parameters(area_ini,'area_remark')
    # grid50784_DEMdiffColor_2008-17
    area_folder = f'{area_name}_{area_remark}_{area_time}'

    # remove data in image_patches
    rm_folder_1 = os.path.join('image_patches',area_folder)
    if os.path.isdir(rm_folder_1):
        io_function.delete_file_or_dir(rm_folder_1)
        print(f'removed {rm_folder_1}')
    else:
        print(f'{rm_folder_1} does not exists')

    # remove results in classified_results
    rm_folder_2 = os.path.join('classified_results','exp3',area_folder)
    if os.path.isdir(rm_folder_2):
        io_function.delete_file_or_dir(rm_folder_2)
        print(f'removed {rm_folder_2}')
    else:
        print(f'{rm_folder_2} does not exist')

    # remove area ini
    io_function.delete_file_or_dir(area_ini)


def main():

    area_ini_list = io_function.read_list_from_txt('area_grid_ini_list.txt')

    wrong_area_ini_list = []

    for area_ini in area_ini_list:
        if check_one_area_ini_file(area_ini) is False:
            wrong_area_ini_list.append(area_ini)
            remove_wrong_file_folder(area_ini)

    print(f'Checked {len(area_ini_list)}, {len(wrong_area_ini_list)} of them are wrong')

if __name__ == '__main__':
    main()

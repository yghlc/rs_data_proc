#!/usr/bin/env python
# Filename: check_valid_tif.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 09 March, 2021
"""

import os,sys
from optparse import OptionParser
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io

machine_name = os.uname()[1]

# some folder paths
if machine_name == 'uist':
    arcticDEM_reg_tif_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
    grid_dem_diff_dir     = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/grid_dem_diffs'
elif machine_name == 'ubuntu':  # tesia
    arcticDEM_reg_tif_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
    grid_dem_diff_dir     = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/grid_dem_diffs'
else:
    raise ValueError('unknown machine:%s' % machine_name)

# check valid tif
def check_tif(tif_list):
    # command_str = 'gdalinfo ' + tif_path
    # res = os.system(command_str)
    # print(tif_path)
    invalid_list = []
    for idx, tif_path in enumerate(tif_list):
        if idx % 50 == 0:
            print('checking %d/%d'%(idx,len(tif_list)))

        try:
            src = raster_io.open_raster_read(tif_path)
        except:
            basic.outputlogMessage(' invalid tif: %s'%tif_path)
            invalid_list.append(tif_path)
    return invalid_list


def main():

    tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir, bsub_folder=False)
    invalid_tif = check_tif(tifs) 
    for tif in invalid_tif:
        print('removing %s'%tif)
        io_function.delete_file_or_dir(tif)
        

    tifs = io_function.get_file_list_by_ext('.tif', grid_dem_diff_dir, bsub_folder=False)
    invalid_tif = check_tif(tifs)
    for tif in invalid_tif:
        print('removing %s'%tif)
        io_function.delete_file_or_dir(tif)


if __name__ == '__main__':
    main()
    pass
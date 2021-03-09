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
def check_tif(tif_path):
    command_str = 'gdalinfo ' + tif_path
    res = os.system(command_str)
    if res!=0:
        print(tif_path, 'is invalid, and will be removed')
        os.system('rm '+tif_path)

def main():

    tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir, bsub_folder=False)
    for tif in tifs:
        check_tif(tif)

    tifs = io_function.get_file_list_by_ext('.tif', grid_dem_diff_dir, bsub_folder=False)
    for tif in tifs:
        check_tif(tif)


if __name__ == '__main__':
    pass
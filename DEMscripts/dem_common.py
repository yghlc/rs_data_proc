#!/usr/bin/env python
# Filename: dem_common.py 
"""
introduction: put some variables

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 March, 2021
"""

import os,sys

machine_name = os.uname()[1]

# some folder paths
if machine_name == 'uist':
    ArcticDEM_tmp_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'

elif machine_name == 'ubuntu':  # tesia
    ArcticDEM_tmp_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'

elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:   # curc
    ArcticDEM_tmp_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir'
else:
    ArcticDEM_tmp_dir = './'

tarball_dir = os.path.join(ArcticDEM_tmp_dir,'tarballs')    # strip version of ArcticDEM

arcticDEM_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir,'registration_tifs')
relative_dem_dir = os.path.join(ArcticDEM_tmp_dir,'dem_relative_8bit')

grid_dem_diffs_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs')
grid_dem_diffs_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_8bit')
grid_dem_diffs_segment_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_8bit')

dem_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'dem_slope_8bit')
dem_hillshade_dir = os.path.join(ArcticDEM_tmp_dir,'dem_hillshade')


if __name__ == '__main__':
    pass
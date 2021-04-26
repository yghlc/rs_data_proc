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
grid_dem_diffs_segment_dir = os.path.join(ArcticDEM_tmp_dir,'grid_dem_diffs_segment_results')

grid_matchtag_sum_dir = os.path.join(ArcticDEM_tmp_dir,'grid_matchtag_sum_tifs')

dem_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'dem_slope_8bit')
dem_slope_dir = os.path.join(ArcticDEM_tmp_dir,'dem_slope')
dem_hillshade_dir = os.path.join(ArcticDEM_tmp_dir,'dem_hillshade')
dem_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_tpi_8bit')

# the mosaic version of AricticDEM
arcticDEM_tile_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_reg_tifs')
arcticDEM_tile_hillshade_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_hillshade')
arcticDEM_tile_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope_8bit')
arcticDEM_tile_slope_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope')
arcticDEM_tile_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_tpi_8bit')
# dem_pattern = '*reg_dem.tif'



grid_20_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp')
dem_strip_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7/ArcticDEM_Strip_Index_Rel7.shp')

ala_north_slo_ext_shp=os.path.expanduser('~/Data/Arctic/alaska/northern_alaska_extent/ala_north_slo_extent_latlon.shp')

if __name__ == '__main__':
    pass
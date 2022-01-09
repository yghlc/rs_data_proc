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
if machine_name == 'uist-int-colorado-edu':
    ArcticDEM_tmp_dir = '/home/lhuang/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
    data_dir = '/home/lhuang/Data'

elif machine_name == 'ubuntu':  # tesia
    ArcticDEM_tmp_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
    data_dir = '/home/lihu9680/Data'

elif 'login' in machine_name or 'shas' in machine_name or 'sgpu' in machine_name:   # curc
    ArcticDEM_tmp_dir = '/scratch/summit/lihu9680/ArcticDEM_tmp_dir'
    data_dir = '/projects/lihu9680/Data'
else:
    ArcticDEM_tmp_dir = './'
    data_dir = './'

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

grd_hillshade_newest_on_top_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_newest_top_grid')

dem_headwall_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_headwall_shp')
grid_dem_headwall_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_headwall_shp_grid')
grid_hillshade_newest_HDLine_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_newest_HWLine_grid')


dem_hillshade_subImages_headwall = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade_subImages_headwall')


grid_dem_subsidence_select = os.path.join(ArcticDEM_tmp_dir,'grid_dem_subsidence_select')

# the mosaic version of AricticDEM
arcticDEM_tile_tarball_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_tarballs')
arcticDEM_tile_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_reg_tifs')
arcticDEM_tile_hillshade_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_hillshade')
arcticDEM_tile_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope_8bit')
arcticDEM_tile_slope_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope')
arcticDEM_tile_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_tpi_8bit')
# dem_pattern = '*reg_dem.tif'


# surface water mask
mask_water_dir = os.path.join(data_dir, 'global_surface_water' , 'extent_epsg3413')

grid_20_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp')
grid_20_id_raster = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km_id.tif')
dem_strip_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7/ArcticDEM_Strip_Index_Rel7.shp')
dem_tile_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Tile_Index_Rel7/ArcticDEM_Tile_Index_Rel7.shp')


# some log and information files
process_log_dir = os.path.join(ArcticDEM_tmp_dir, 'log_dir')
grid_complete_list_txt = os.path.join(process_log_dir,'grid_complete_ids.txt')  # store ids of grids that have completed
# manually exclude some grids that dont have enough data
grid_excluded_list_txt = os.path.join(process_log_dir,'grid_exclude_ids.txt')   # store ids of grids that manually exclude

grid_dem_diff_less2dem_txt = os.path.join(process_log_dir,'grid_dem_diff_less2dem_ids.txt')  # store ids of grids that has less than 2 DEM (cannot calculate DEM differnce)
grid_no_dem_txt = os.path.join(process_log_dir,'grid_no_dem_ids.txt')  # store ids of grids that don't have DEM (strip version) overlap

# store ids of grids that have overlap of DEM  (strip version), but the coverage is too smaller or all overlap DEM are invalid.
grid_no_valid_dem_txt = os.path.join(process_log_dir,'grid_no_valid_dem_ids.txt')

# some place that is really flat, cannot detect headwall based on slope from it
grid_no_headwall_txt = os.path.join(process_log_dir,'grid_no_headwall_ids.txt')

strip_dem_cover_grids_txt = os.path.join(process_log_dir,'strip_dem_cover_grids.txt') # each strip cover how many grids (ids), dict
tile_dem_cover_grids_txt = os.path.join(process_log_dir,'tile_dem_cover_grids.txt') # each tile cover how many grids (ids), dict

# rts results
grid_rts_shp_dir = os.path.join(ArcticDEM_tmp_dir, 'grid_rts_shp')


if __name__ == '__main__':
    pass
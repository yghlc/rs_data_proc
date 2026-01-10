#!/usr/bin/env python
# Filename: rm_no_need_dem_files.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 17 June, 2022
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.timeTools as timeTools


from process_largeRegion_butLimited_storage import remove_no_need_dem_files
from process_largeRegion_butLimited_storage import update_complete_grid_list

from dem_common import get_grid_id_from_path
from dem_common import arcticDEM_reg_tif_dir, grid_dem_diffs_dir

def remove_on_need_stac_dem_after_DEM_diff(bak_dir=None):
    # remove no need dem files after the DEM diff has been produced
    dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir,'*DEM_diff*.tif')
    grid_ids = [ get_grid_id_from_path(item) for item in dem_diff_list]

    # find dem files
    for g_id in grid_ids:
        tif_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,f'*_grid{g_id}_*.tif')
        if len(tif_list) > 0:
            if bak_dir is None:
                # remove them directly
                basic.outputlogMessage(f'To remove {len(tif_list)} tifs for grid: {g_id}')
                for tif in tif_list:
                    io_function.delete_file_or_dir(tif)
            else:
                # move to the backup directory
                basic.outputlogMessage(f'To move {len(tif_list)} tifs for grid: {g_id} to folder: {bak_dir}')
                for tif in tif_list:
                    io_function.movefiletodir(tif,bak_dir,overwrite=False,b_verbose=False)

    basic.outputlogMessage('completed: remove_on_need_stac_dem_after_DEM_diff')


    pass

def main():
    basic.setlogfile('rm_no_need_dem_files.log')
    # check completed list
    # update_complete_grid_list(grid_ids, task_list)  # need to update complete list on the main pre-processing workstation first.
    # this remove tarball and tif download from PGC.
    # remove_no_need_dem_files()

    backup_folder = os.path.join(arcticDEM_reg_tif_dir, f'rm_{timeTools.get_now_date_str()}')
    io_function.mkdir(backup_folder)

    # remove the file download by stac
    remove_on_need_stac_dem_after_DEM_diff(bak_dir=backup_folder)

if __name__ == '__main__':
    main()
    pass
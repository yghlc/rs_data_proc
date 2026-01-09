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


from process_largeRegion_butLimited_storage import remove_no_need_dem_files
from process_largeRegion_butLimited_storage import update_complete_grid_list

from dem_common import get_grid_id_from_path

def remove_on_need_stac_dem_after_DEM_diff():
    # remove no need dem files after the DEM diff has been produced
    from dem_common import arcticDEM_reg_tif_dir, grid_dem_diffs_dir
    dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir,'*DEM_diff*.tif')
    grid_ids = [ get_grid_id_from_path(item) for item in dem_diff_list]

    # find dem files
    for g_id in grid_ids:
        tif_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,f'*_grid{g_id}_*.tif')
        if len(tif_list) > 0:
            basic.outputlogMessage(f'To remove {len(tif_list)} for grid: {g_id}')
            for tif in tif_list:
                io_function.delete_file_or_dir(tif)




    pass

def main():
    basic.setlogfile('rm_no_need_dem_files.log')
    # check completed list
    # update_complete_grid_list(grid_ids, task_list)  # need to update complete list on the main pre-processing workstation first.
    # this remove tarball and tif download from PGC.
    # remove_no_need_dem_files()

    # remove the file download by stac
    remove_on_need_stac_dem_after_DEM_diff()

if __name__ == '__main__':
    main()
    pass
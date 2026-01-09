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


from process_largeRegion_butLimited_storage import remove_no_need_dem_files
from process_largeRegion_butLimited_storage import update_complete_grid_list

from dem_common import get_grid_id_from_path

def remove_on_need_stac_dem_after_DEM_diff():
    # remove no need dem files after the DEM diff has been produced
    from dem_common import arcticDEM_reg_tif_dir, grid_dem_diffs_dir
    dem_diff_list = io_function.get_file_list_by_pattern(grid_dem_diffs_dir,'')



    pass

def main():
    # check completed list
    # update_complete_grid_list(grid_ids, task_list)  # need to update complete list on the main pre-processing workstation first.
    # this remove tarball and tif download from PGC.
    remove_no_need_dem_files()

if __name__ == '__main__':
    main()
    pass
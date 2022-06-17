#!/usr/bin/env python
# Filename: rm_no_need_dem_files.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 17 June, 2022
"""

from process_largeRegion_butLimited_storage import remove_no_need_dem_files
from process_largeRegion_butLimited_storage import update_complete_grid_list

def main():
    # check completed list
    # update_complete_grid_list(grid_ids, task_list)  # need to update complete list on the main pre-processing workstation first.
    remove_no_need_dem_files()

if __name__ == '__main__':
    remove_no_need_dem_files()
    pass
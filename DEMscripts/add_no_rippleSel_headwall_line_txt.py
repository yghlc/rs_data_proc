#!/usr/bin/env python
# Filename: add_no_rippleSel_headwall_line_txt.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 20 July, 2022
"""

# run this script in /Miavaig/Work/lingcaoHuang/ArcticDEM_results on tesia

import os,sys
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function

from dem_common import grid_no_rippleSel_headwall_line_txt
from dem_common import get_grid_id_from_path
from identify_headwall_lines import no_result

from datetime import datetime
import time

def main():
    # read all no no_rippleSel_headwall_line grids
    no_res_grid = io_function.read_list_from_txt(grid_no_rippleSel_headwall_line_txt)

    ext_list = io_function.get_file_list_by_pattern('./','ext??_*')
    print('region folders: ')
    [print(item) for item in ext_list]

    curr_dir = os.getcwd()

    for region in ext_list:
        print('checking region: ', region)
        time.sleep(3)
        rippleSel_list = io_function.get_file_list_by_pattern(os.path.join(region, 'dem_headwall_shp_grid'),
                                                           'headwall_shps_grid*/*_rippleSel.shp')
        if len(rippleSel_list) < 1:
            print('this region may have not been completed, skip')
            continue

        folder_list = io_function.get_file_list_by_pattern(os.path.join(region,'dem_headwall_shp_grid'),'headwall_shps_grid*')
        for item in folder_list:
            print('checking folder:', item)
            shp_list = io_function.get_file_list_by_pattern(item,'*.shp')
            if len(shp_list) >= 2:
                continue
            elif len(shp_list) == 1:
                grid_id = get_grid_id_from_path(item)
                no_result_txt = os.path.join(item,no_result)
                if str(grid_id) in no_res_grid:
                    io_function.save_list_to_txt(no_result_txt, [str(grid_id)])
                    print(datetime.now(),'saving %s to %s'%(no_result,item))
            else:
                pass


if __name__ == '__main__':
    main()
    pass
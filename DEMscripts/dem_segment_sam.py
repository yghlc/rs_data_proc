#!/usr/bin/env python
# Filename: dem_segment_sam.py 
"""
introduction: using Segment anything model to segment elevation reduction using

~/codes/PycharmProjects/BigImageMapper/sam_dir/sam_predict.py
Please following the files (*.sh, *.ini) in ~/codes/PycharmProjects/BigImageMapper/sam_dir to set input files

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 February, 2025
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function

# ArcticDEM results on ygAlpha
ArcticDEM_results_dir = os.path.expanduser('~/data1/ArcticDEM_results')

def sam_sgment_a_big_region(dem_diff_dir, save_dir, tmp_dir):
    pass

def test_sam_sgment_a_big_region():
    dem_diff_dir = os.path.join(ArcticDEM_results_dir,'./ext09_for_ArcticDEM_proc','grid_dem_diffs')


def main():

    dem_diff_dir_list = io_function.get_file_list_by_pattern(ArcticDEM_results_dir,'ext??_*/*diffs*')
    print(dem_diff_dir_list)


    pass

if __name__ == '__main__':
    main()
#!/usr/bin/env python
# Filename: image_segment_test.py 
"""
introduction: pytest of image_segment.py

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 February, 2021
"""

import os,sys
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import image_segment
import raster_io
import basic_src.io_function as io_function

def test_watershed_segmentation():
    # img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/north_ext_ArcticDEM_diff_sub_1.tif'
    #                               'dem_north/north_ext_ArcticDEM_diff_sub_1_clip2.tif')
    # img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/dem_north/north_ext_ArcticDEM_diff_sub_1.tif')
    img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/'
                                  'dem_north/north_ext_ArcticDEM_diff_sub_1_clip2_8bit.tif')

    io_function.is_file_exist(img_path)
    save_dir = './'
    image_segment.segment_changes_on_dem_diff(img_path, save_dir)

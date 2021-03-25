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

import numpy as np

# def test_watershed_segmentation():
#     # img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/north_ext_ArcticDEM_diff_sub_1.tif'
#     #                               'dem_north/north_ext_ArcticDEM_diff_sub_1_clip2.tif')
#     # img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/dem_north/north_ext_ArcticDEM_diff_sub_1.tif')
#     img_path = os.path.expanduser('~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/'
#                                   'dem_north/north_ext_ArcticDEM_diff_sub_1_clip2_8bit.tif')
#
#     io_function.is_file_exist(img_path)
#     save_dir = './'
#     image_segment.segment_changes_on_dem_diff(img_path, save_dir)


def test_to_unique_label_for_superpixels():
    label_img = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/segment_parallel_sub/WR_dem_diff_DEM_diff_prj_8bit_sub_label.tif')
    out_labels, nodata = raster_io.read_raster_one_band_np(label_img)
    print('nodata',nodata)
    print('min and max labels of out_labels', np.min(out_labels), np.max(out_labels))
    new_labels = image_segment.to_unique_label_for_superpixels(out_labels)

    save_new_label = io_function.get_name_by_adding_tail(label_img,'new')
    raster_io.save_numpy_array_to_rasterfile(new_labels,save_new_label,label_img)

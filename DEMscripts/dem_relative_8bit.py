#!/usr/bin/env python
# Filename: dem_relative_8bit.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 20 March, 2021
"""

import os,sys
from optparse import OptionParser
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
import split_image
import raster_io


import multiprocessing
from multiprocessing import Pool

machine_name = os.uname()[1]

import numpy as np

# some folder paths
if machine_name == 'uist':
    ArcticDEM_tmp_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
elif machine_name == 'ubuntu':  # tesia
    ArcticDEM_tmp_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
else:
    ArcticDEM_tmp_dir = './'

arcticDEM_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir,'registration_tifs')
relative_dem_dir = os.path.join(ArcticDEM_tmp_dir,'dem_relative_8bit')


def dem_to_relative_8bit_a_patch(idx, patch, patch_count,dem_path, dst_nodata):
    # print('tile: %d / %d' % (idx + 1, patch_count))
    # patch_w = patch[2]
    # patch_h = patch[3]

    dem, nodata = raster_io.read_raster_one_band_np(dem_path, boundary=patch)
    # print(dem.shape)
    # print(dem.ndim)
    # dem_re = np.expand_dims(dem,axis=0)
    # print(dem_re.shape)

    patch_relative_dem_8bit = raster_io.image_numpy_allBands_to_8bit_hist(dem,per_min=0.02, per_max=0.98,
                                                                          src_nodata=nodata,dst_nodata=dst_nodata)

    return patch,patch_relative_dem_8bit


def dem_to_relative_dem(input,output, patch_width, patch_height, process_num):
    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    height, width, _,_ = raster_io.get_height_width_bandnum_dtype(input)
    dst_nodata = 255

    # divide the image the many small patches, then calcuate one by one, solving memory issues.
    image_patches = split_image.sliding_window(width,height, patch_width, patch_height,adj_overlay_x=0,adj_overlay_y=0)
    patch_count = len(image_patches)

    # get the difference
    dem_relative_8bit_np = np.zeros((height, width),dtype=np.uint8)

    if process_num == 1:
        for idx, patch in enumerate(image_patches):
            _, patch_rel_dem_8bit = dem_to_relative_8bit_a_patch(idx, patch, patch_count,input,dst_nodata)
            # copy to the entire image
            row_s = patch[1]
            row_e = patch[1] + patch[3]
            col_s = patch[0]
            col_e = patch[0] + patch[2]
            dem_relative_8bit_np[row_s:row_e, col_s:col_e] = patch_rel_dem_8bit
    else:
        theadPool = Pool(process_num)
        parameters_list = [(idx, patch, patch_count, input,dst_nodata) for idx, patch in enumerate(image_patches)]
        results = theadPool.starmap(dem_to_relative_8bit_a_patch, parameters_list)
        for res in results:
            patch, patch_rel_dem_8bit = res
            # copy to the entire image
            row_s = patch[1]
            row_e = patch[1] + patch[3]
            col_s = patch[0]
            col_e = patch[0] + patch[2]
            dem_relative_8bit_np[row_s:row_e, col_s:col_e] = patch_rel_dem_8bit
        theadPool.close()
    # save date diff to tif (16 bit)
    raster_io.save_numpy_array_to_rasterfile(dem_relative_8bit_np,output,input, nodata=dst_nodata,compress='lzw',tiled='yes',bigtiff='if_safer')

    return True

def main():
    basic.setlogfile('log_to_relative_dem_8bit.txt')

    if os.path.isdir(relative_dem_dir) is False:
        io_function.mkdir(relative_dem_dir)

    # 500 pixel by 500 pixel, that is 1 km by 1 km
    patch_width = 500
    patch_height = 500
    process_num = 1

    failed_tifs = []

    dem_reg_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir, '*dem_reg.tif')
    count = len(dem_reg_list)
    for idx, tif in enumerate(dem_reg_list):
        print('%d/%d convert %s to relative DEM (8bit)' % (idx + 1, count, tif))
        rel_dem_8bit = io_function.get_name_by_adding_tail(tif,'relDEM8bit')
        rel_dem_8bit = os.path.join(relative_dem_dir, os.path.basename(rel_dem_8bit))
        try:
            dem_to_relative_dem(tif,rel_dem_8bit,patch_width,patch_height,process_num)
        except:
            failed_tifs.append(tif)

    with open('to_relative_dem_failed_cases.txt','w') as f_obj:
        for item in failed_tifs:
            f_obj.writelines(item + '\n')
    pass

if __name__ == '__main__':
    main()
    pass
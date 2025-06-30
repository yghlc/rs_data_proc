#!/usr/bin/env python
# Filename: convert_Landslide4Sense.py 
"""
introduction: convert the h5 files in Landslide4Sense to three bands (RGB) or 4 band images (NIR + RGB)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 30 June, 2025
"""


import os,sys
import h5py

import numpy as np
from overrides.signature import ensure_all_kwargs_defined_in_sub

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
# import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function

# from: https://zenodo.org/records/10463239
# The Landslide4Sense dataset has three splits, training/validation/test, consisting of 3799, 245, and 800 image patches,
# respectively. Each image patch is a composite of 14 bands that include:
# Multispectral data from Sentinel-2: B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12.
# Slope data from ALOS PALSAR: B13.
# Digital elevation model (DEM) from ALOS PALSAR: B14

def extract_bands_to_8bit(h5_path, band_list=None, hdf5_key='img', save_path=None,b_8bit=False):
    # B2, B3, B4, B8 are Blue, Green, Red, NIR
    # 1, 2, 3, 7
    if band_list is None:
        band_list = [1,2,3]  # save B2, B3, B4, that is Blue, Green, Red

    if save_path is None:
        save_path = io_function.get_name_no_ext(h5_path) + '.tif'

    with h5py.File(h5_path, 'r') as hf:
        print(hf.keys())
        image = hf[hdf5_key][:]
        # print(image.shape)
        if image.ndim == 2:
            image_sel = image # no need to select bands
        else:
            image_sel = image[:,:,band_list]
            # print(image_sel.shape)
            image_sel = np.transpose(image_sel, (2, 0, 1))  # rasterio accept: band, width, height

        # convert to 8 bit
        if b_8bit:
            image_sel = raster_io.image_numpy_allBands_to_8bit_hist(image_sel,per_min=0.02, per_max=0.98)

        # save to disk
        # image_sel = np.transpose(image_sel, (2, 0, 1))  # rasterio accept: band, width, height
        raster_io.save_numpy_array_to_rasterfile_noCRS(image_sel,save_path)


def test_extract_bands_to_8bit():
    data_dir = os.path.expanduser('~/Data/public_data_AI/Landslide4Sense/TrainData')
    img_h5 = os.path.join(data_dir,'img','image_1529.h5')
    mask_h5 = os.path.join(data_dir,'mask','mask_1529.h5')

    # extract_bands_to_8bit(img_h5,b_8bit=True)
    # convert the mask
    extract_bands_to_8bit(mask_h5,hdf5_key='mask')


def convert_to_RGB_8bit():
    data_dir = os.path.expanduser('~/Data/public_data_AI/Landslide4Sense')
    out_dir = os.path.expanduser('~/Data/public_data_AI/Landslide4Sense/rgb_8bit')
    datasets = ['TrainData','ValidData','TestData']
    for split in datasets:
        img_dir = os.path.join(data_dir,split,'img')
        mask_dir = os.path.join(data_dir,split,'mask')

        save_img_dir = os.path.join(out_dir,split,'img')
        save_mask_dir = os.path.join(out_dir,split,'mask')
        if os.path.isdir(save_img_dir) is False:
            io_function.mkdir(save_img_dir)
        if os.path.isdir(save_mask_dir) is False:
            io_function.mkdir(save_mask_dir)

        img_h5_list = io_function.get_file_list_by_pattern(img_dir,'*.h5')
        print(f'Found {len(img_h5_list)} h5 files in {img_dir}')
        for h5 in img_h5_list:
            save_file_name = io_function.get_name_no_ext(h5) + '.tif'
            save_8bit_tif = os.path.join(save_img_dir, save_file_name)
            extract_bands_to_8bit(h5, b_8bit=True, save_path=save_8bit_tif)

        mask_h5_list = io_function.get_file_list_by_pattern(mask_dir, '*.h5')
        for h5 in mask_h5_list:
            print(f'Found {len(mask_h5_list)} h5 files in {mask_dir}')
            save_file_name = io_function.get_name_no_ext(h5) + '.tif'
            save_tif = os.path.join(save_mask_dir, save_file_name)
            extract_bands_to_8bit(h5, hdf5_key='mask', save_path=save_tif)



def main():
    # test_extract_bands_to_8bit()

    convert_to_RGB_8bit()

    pass


if __name__ == '__main__':
    main()

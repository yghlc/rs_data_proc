#!/usr/bin/env python
# Filename: dilate_img_valid_pixels.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 05 December, 2025
"""

import os, sys
from optparse import OptionParser

deeplabRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabRS)
import raster_io

import numpy as np
from scipy.ndimage import generic_filter

def dilation_majority_filter(image, size=3, nodata=255):
    """
    Expands valid pixels into nodata using a majority filter.

    Args:
        image: 2D numpy array.
        size: Size of the square filter window (must be odd).
        nodata: Value representing nodata.
    Returns:
        New 2D numpy array after dilation.
    """
    def majority_func(values):
        center = values[len(values) // 2]
        if center != nodata:
            return center  # keep valid pixels unchanged
        # Count only valid values
        valid = values[values != nodata]
        if valid.size == 0:
            return nodata  # all nodata in window
        # Find mode (majority)
        vals, counts = np.unique(valid, return_counts=True)
        majority_val = vals[np.argmax(counts)]
        return majority_val

    # Apply the generic filter
    footprint = np.ones((size, size), dtype=bool)
    return generic_filter(image, majority_func, footprint=footprint, mode='constant', cval=nodata)

def test_dilation_majority_filter():
    # === Example usage ===
    # Let's assume your image is loaded as a numpy array:
    # image = np.load('your_image.npy')
    # For demonstration, create a small sample array:
    image = np.array([
        [255, 1, 1, 255],
        [255, 2, 1, 255],
        [255, 255, 255, 255],
        [255, 3, 3, 255],
    ], dtype=np.uint8)

    result = dilation_majority_filter(image, size=3, nodata=255)
    print(result)

def main(options, args):

    input_path = args[0]
    save_path = args[1]
    nodata = options.nodata
    filter_size = options.filter_size
    # test_dilation_majority_filter()

    img_data, img_nodata = raster_io.read_raster_one_band_np(input_path,band=1)
    print('img_data shape:', img_data.shape)
    output_data = dilation_majority_filter(img_data,size=filter_size,nodata=nodata)
    raster_io.save_numpy_array_to_rasterfile(output_data,save_path,input_path,format='GTiff',
                                             nodata=nodata, compress='lzw', tiled='yes', bigtiff='if_safer')


if __name__ == '__main__':
    usage = "usage: %prog [options] image_path save_path"
    parser = OptionParser(usage=usage, version="1.0 2021-02-07")
    parser.description = 'convert images to unsigned 8bit '

    parser.add_option("-f", "--filter_size",
                      action="store", dest="filter_size", type=int, default=3,
                      help="filter size, such as 3 (3 by 3)")

    parser.add_option("-n", "--nodata",
                      action="store", dest="nodata", type=int, default=255,
                      help="the value for nodata")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 2:
        parser.print_help()
        sys.exit(2)

    main(options, args)

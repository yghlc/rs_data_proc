#!/usr/bin/env python
# Filename: grey_image_segment.py 
"""
introduction: segment a grey image (8 bit)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.io_function as io_function
import split_image

from image_segment import quickshift_segmentaion

import multiprocessing
from multiprocessing import Pool
import numpy as np

def segment_a_patch(idx, patch, patch_count,img_path):

    print('tile: %d / %d' % (idx + 1, patch_count))
    # read imag
    one_band_img, nodata = raster_io.read_raster_one_band_np(img_path, boundary=patch)

    # segmentation algorithm (the output of these algorithms is not alway good, need to chose the parameters carafully)
    # out_labels = watershed_segmentation(one_band_img)
    # out_labels = k_mean_cluster_segmentation(one_band_img)
    out_labels = quickshift_segmentaion(one_band_img,ratio=1.0, kernel_size=5, max_dist=20,
                           sigma=0, convert2lab=False)

    return patch, out_labels, nodata

def segment_a_grey_image(img_path, save_dir,process_num, dem_diff=None):

    out_pre = os.path.splitext(os.path.basename(img_path))[0]
    height, width, band_num, date_type = raster_io.get_height_width_bandnum_dtype(img_path)
    print('input image: height, width, band_num, date_type',height, width, band_num, date_type)

    save_labes = np.zeros((height,width),dtype=np.int32)
    # divide the image the many small patches, then calcuate one by one, solving memory issues.
    image_patches = split_image.sliding_window(width,height, 1024, 1024,adj_overlay_x=0,adj_overlay_y=0)
    patch_count = len(image_patches)

    # for idx, patch in enumerate(image_patches):
    #     out_patch,out_labels = segment_a_patch(idx, patch, patch_count,img_path)
    #     # copy to the entire image
    #     row_s = patch[1]
    #     row_e = patch[1] + patch[3]
    #     col_s = patch[0]
    #     col_e = patch[0] + patch[2]
    #     save_labes[row_s:row_e, col_s:col_e] = out_labels

    theadPool = Pool(process_num)
    parameters_list = [ (idx, patch, patch_count,img_path) for idx, patch in enumerate(image_patches)]
    results = theadPool.starmap(segment_a_patch, parameters_list)
    current_min = 0
    for res in results:
        patch, out_labels, nodata = res
        # copy to the entire image
        row_s = patch[1]
        row_e = patch[1] + patch[3]
        col_s = patch[0]
        col_e = patch[0] + patch[2]
        current_min = np.max(save_labes)
        print('current_max',current_min)
        save_labes[row_s:row_e, col_s:col_e] = out_labels + current_min

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    # save the label
    label_path = os.path.join(save_dir, out_pre + '_label.tif')
    raster_io.save_numpy_array_to_rasterfile(save_labes, label_path, img_path) # do not set nodata

    # convert the label to shapefile
    out_shp = os.path.join(save_dir, out_pre + '.shp')
    command_string = 'gdal_polygonize.py -8 %s -b 1 -f "ESRI Shapefile" %s' % (label_path, out_shp)
    res = os.system(command_string)
    if res != 0:
        sys.exit(1)

    if dem_diff is not None:
        # do something of polygon merging
        pass

def main(options, args):

    img_path = args[0]
    io_function.is_file_exist(img_path)
    save_dir = options.save_dir
    process_num = options.process_num
    org_elevation_diff = options.elevation_diff

    segment_a_grey_image(img_path,save_dir,process_num,dem_diff=org_elevation_diff)


if __name__ == "__main__":
    usage = "usage: %prog [options] image_path "
    parser = OptionParser(usage=usage, version="1.0 2021-3-1")
    parser.description = 'Introduction: segment a grey image using watershed algorithm   '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to create the mosaic")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='./',
                      help="the folder to save results")

    parser.add_option("-c", "--elevation_diff",
                      action="store", dest="elevation_diff",
                      help="the original elevation difference")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

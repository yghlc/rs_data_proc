#!/usr/bin/env python
# Filename: watershed_segment 
"""
introduction: segment a grey image

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 February, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.io_function as io_function

import numpy as np
from skimage import segmentation

def watershed_segmentation(img_2d_one_band):
    # the result is so bad when we applied this one to dem difference (float32), maybe we need choose some markers.
    out_labels = segmentation.watershed(img_2d_one_band)
    return out_labels

def k_mean_cluster_segmentation(img_2d_one_band):
    # n_segments=100, hard to decide this parameter.
    # result is bad when applied to the dem difference (float32), need to adjust the parameters.
    img_2d_one_band = img_2d_one_band.astype(np.float64)
    out_labels = segmentation.slic(img_2d_one_band)
    out_labels = out_labels.astype(np.int32)
    return out_labels

def quickshift_segmentaion(img_2d_one_band, ratio=1.0, kernel_size=5, max_dist=10, return_tree=False,
                           sigma=0, convert2lab=True, random_seed=42):
    # Segments image using quickshift clustering in Color-(x,y) space.
    # Produces an oversegmentation of the image using the quickshift mode-seeking algorithm.

    # https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_segmentations.html?highlight=segmentation
    # Quickshift has two main parameters: sigma controls the scale of the local density approximation,
    # max_dist selects a level in the hierarchical segmentation that is produced.
    # There is also a trade-off between distance in color-space and distance in image-space, given by ratio.

    # ValueError: the input array must be have a shape == (.., ..,[ ..,] 3)), got (1882, 1895, 1), set convert2lab=False

    img_2d_one_band = img_2d_one_band.astype(np.float64)
    out_labels = segmentation.quickshift(img_2d_one_band, ratio=ratio, kernel_size=kernel_size, max_dist=max_dist, return_tree=return_tree,
                           sigma=sigma, convert2lab=convert2lab, random_seed=random_seed)
    out_labels = out_labels.astype(np.int32)
    return out_labels

def segment_changes_on_dem_diff(dem_diff_tif, save_dir):

    out_pre = os.path.splitext(os.path.basename(dem_diff_tif))[0]

    # read images
    one_band_img, nodata = raster_io.read_raster_one_band_np(dem_diff_tif)

    # segmentation algorithm (the output of these algorithms is not alway good, need to chose the parameters carafully)
    # out_labels = watershed_segmentation(one_band_img)
    # out_labels = k_mean_cluster_segmentation(one_band_img)
    # out_labels = quickshift_segmentaion(one_band_img)

    # segmentation by threshold (may have too many noise)
    mean = np.nanmean(one_band_img)
    print("mean value is: %.4f"%mean)
    one_band_img = one_band_img - mean
    out_labels = np.zeros_like(one_band_img,dtype=np.uint8)
    out_labels[ np.abs(one_band_img) > 2 ] = 1

    # save the label
    label_path = os.path.join(save_dir, out_pre + '_label.tif')
    raster_io.save_numpy_array_to_rasterfile(out_labels, label_path, dem_diff_tif, nodata=0)

    # convert the label to shapefile
    out_shp = os.path.join(save_dir, out_pre + '.shp')
    command_string = 'gdal_polygonize.py -8 %s -b 1 -f "ESRI Shapefile" %s' % (label_path, out_shp)
    res = os.system(command_string)
    if res != 0:
        sys.exit(1)


def main(options, args):

    img_path = args[0]
    io_function.is_file_exist(img_path)
    save_dir = options.save_dir
    segment_changes_on_dem_diff(img_path,save_dir)

    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] image_path "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment a grey image using watershed algorithm   '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=1,
                      help="number of processes to create the mosaic")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save results")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

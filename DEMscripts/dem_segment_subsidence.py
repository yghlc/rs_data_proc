#!/usr/bin/env python
# Filename: watershed_segment 
"""
introduction: segment a grey image

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 1 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection

import cv2
import numpy as np
import pandas as pd

def post_processing_subsidence(in_shp):
    polygons = vector_gpd.read_polygons_gpd(in_shp)

    # get shapeinfo
    # poly_shapeinfo_list = []
    save_polyons = []
    for poly in polygons:
        # get INarea, INperimete, WIDTH, HEIGHT, ratio_w_h, hole_count
        # shapeinfo = vector_gpd.calculate_polygon_shape_info(poly)     # error:  'MultiPolygon' object has no attribute 'interiors'
        # poly_shapeinfo_list.append(shapeinfo)
        # if shapeinfo['INarea'] < 40:    # remove the one with area smaller than 40 m^2
        if poly.area < 90:    # remove the one with area smaller than 40 m^2
            continue
        save_polyons.append(poly)

    save_pd = pd.DataFrame({'Polygon': save_polyons})
    wkt = map_projection.get_raster_or_vector_srs_info_wkt(in_shp)
    save_shp = io_function.get_name_by_adding_tail(in_shp,'post')
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,save_shp)




def segment_subsidence_on_dem_diff(dem_diff_tif, save_dir):

    out_pre = os.path.splitext(os.path.basename(dem_diff_tif))[0]

    # read images
    one_band_img, nodata = raster_io.read_raster_one_band_np(dem_diff_tif)

    # segmentation by threshold (may have too many noise)
    # mean = np.nanmean(one_band_img)
    # print("mean value is: %.4f"%mean)
    # one_band_img = one_band_img - mean    # cannot use mean which may affect by some Outliers
    out_labels = np.zeros_like(one_band_img,dtype=np.uint8)
    out_labels[ one_band_img < -1 ] = 1     # end in a lot of noise

    # apply median filter
    out_labels = cv2.medianBlur(out_labels, 3)  # with kernal=3

    # save the label
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    label_path = os.path.join(save_dir, out_pre + '_label.tif')
    raster_io.save_numpy_array_to_rasterfile(out_labels, label_path, dem_diff_tif, nodata=0)

    # convert the label to shapefile
    out_shp = os.path.join(save_dir, out_pre + '.shp')
    command_string = 'gdal_polygonize.py -8 %s -b 1 -f "ESRI Shapefile" %s' % (label_path, out_shp)
    res = os.system(command_string)
    if res != 0:
        sys.exit(1)

    # post-processing
    post_processing_subsidence(out_shp)

def main(options, args):

    img_path = args[0]
    io_function.is_file_exist(img_path)
    save_dir = options.save_dir
    segment_subsidence_on_dem_diff(img_path,save_dir)

    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] image_path "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment subsidence based on DEM difference  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save results")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

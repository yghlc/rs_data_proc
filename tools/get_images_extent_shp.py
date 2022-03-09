#!/usr/bin/env python
# Filename: get_images_extent_shp.py 
"""
introduction: get images exetents

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 24 April, 2021
"""
import sys,os
from optparse import OptionParser

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import vector_gpd
import raster_io

import pandas as pd

def main(options, args):
    img_dir = args[0]
    output = options.output
    if os.path.isdir(img_dir):
        img_pattern = options.image_pattern
        if output is None:
            output = os.path.basename(img_dir) + '_boxes.shp'
        img_list = io_function.get_file_list_by_pattern(img_dir, img_pattern)
        if len(img_list) < 1:
            raise ValueError('No images in %s with pattern: %s'%(img_dir, img_pattern))
    else:
        # if it's a file
        img_list = [img_dir]
        if output is None:
            output = os.path.basename(img_dir) + '_bound.shp'

    print('Find %d rasters in %s'%(len(img_list),img_dir))

    # check projections?
    prj4_1st = raster_io.get_projection(img_list[0],'proj4')
    for idx in range(1, len(img_list)):
        prj4 = raster_io.get_projection(img_list[idx],'proj4')
        if prj4_1st != prj4:
            raise ValueError('projection inconsistent between %s and %s '%(img_list[0], img_list[idx]))


    img_boxes = [ raster_io.get_image_bound_box(img_path)  for img_path in img_list]
    img_box_polygons = [ vector_gpd.convert_image_bound_to_shapely_polygon(box) for box in img_boxes]

    # save to file
    wkt = map_projection.get_raster_or_vector_srs_info_proj4(img_list[0])

    save_pd = pd.DataFrame({'raster':img_list, 'Polygon':img_box_polygons})
    vector_gpd.save_polygons_to_files(save_pd,'Polygon',wkt,output)
    print('save raster extents to %s'%output)

    return

if __name__ == "__main__":

    usage = "usage: %prog [options] image_dir or image_path"
    parser = OptionParser(usage=usage, version="1.0 2021-04-24")
    parser.description = 'Introduction: get image list and their extent'

    parser.add_option("-e", "--image_pattern",
                      action="store", dest="image_pattern", default='*.tif',
                      help="the image pattern of the image file")

    parser.add_option("-o", "--output",
                      action="store", dest="output",
                      help="the output shapefile")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
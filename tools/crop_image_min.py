#!/usr/bin/env python
# Filename: crop_image_min 
"""
introduction: crop a image based on a given polygon to the mim extent
# similar to gdalwarp, but gdalwarp only crop the the extent of the polygon.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 21 February, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.RSImageProcess as RSImageProcess
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import vector_gpd

def subset_image_by_polygon_min(in_img, out_img, polygon,resample_m='bilinear',o_format='GTiff', out_res=None):
    # if same_extent:
    #     return RSImageProcess.subset_image_by_polygon_box(out_img,in_img,polygon,resample_m=resample_m, o_format=o_format, xres=out_res,yres=out_res)
    # else:
    # crop to the min extent (polygon or the image)


    return RSImageProcess.subset_image_by_polygon_box_image_min(out_img,in_img,polygon,resample_m=resample_m,o_format=o_format,
                                                                xres=out_res,yres=out_res,compress='lzw', tiled='yes', bigtiff='if_safer')


def main(options, args):

    extent_shp = args[0]
    img_path = args[1]
    save_dir = options.save_dir

    #check projection
    extent_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp)
    img_prj = map_projection.get_raster_or_vector_srs_info_proj4(img_path)
    if img_prj != extent_prj:
        raise ValueError('Project of %s and %s is different'%(extent_shp, img_path))

    out_img = io_function.get_name_by_adding_tail(img_path,'sub')
    out_img = os.path.join(save_dir, os.path.basename(out_img))

    extent_polys = vector_gpd.read_polygons_gpd(extent_shp)
    if len(extent_polys) != 1:
        raise ValueError('current only support one polygon')

    for ext_poly in extent_polys:
        subset_image_by_polygon_min(img_path, out_img, ext_poly, resample_m='bilinear', o_format='GTiff', out_res=None)




    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp image_path"
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: crop image to a min extent   '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")
    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if  len(args) < 2:
        parser.print_help()
        sys.exit(2)

    main(options, args)




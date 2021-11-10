#!/usr/bin/env python
# Filename: get_image_valid_extent.py 
"""
introduction:  get image valid region (excluding nodata pixels),
            "gdaltindex" or "bounding boxes of images" contains all pixels

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 10 November, 2021
"""
import sys,os
from optparse import OptionParser

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function

import raster_io

def get_valid_pixel_mask(in_raster, nodata=None,out_mask=None):

    if nodata is None:
        nodata = raster_io.get_nodata(in_raster)
    if nodata is None:
        raise ValueError('No data is not set in the raster or assigned')

    if out_mask is None:
        out_mask = 'tmp.tif'

    # get nodata mask from ref_raster
    command_str = 'gdal_calc.py --calc="A!=%d" --outfile=%s -A %s --NoDataValue 0 --type=Byte '%(nodata,out_mask, in_raster)

    res = os.system(command_str)
    if res != 0:
        print(res)
        sys.exit(1)

    return out_mask


def main(options, args):
    img_path = args[0]
    # if nodata is not set, it will try to read from images.
    img_nodata = options.nodata
    output_shp = options.output_shp
    if output_shp is None:
        output_shp = os.path.splitext(os.path.basename(img_path))[0] + '_valid.shp'

    #
    out_mask = options.output_mask

    valid_mask = get_valid_pixel_mask(img_path,img_nodata,out_mask=out_mask)
    # to shapefile

    raster_io.raster2shapefile(valid_mask,out_shp=output_shp,nodata=0)  # the nodata for valid_mask is 0

    if valid_mask=='tmp.tif':
        io_function.delete_file_or_dir(valid_mask)


if __name__ == "__main__":

    usage = "usage: %prog [options] raster_path "
    parser = OptionParser(usage=usage, version="1.0 2021-11-04")
    parser.description = 'Introduction: get valid region of a raster'

    parser.add_option("-n", "--nodata",
                      action="store", dest="nodata",
                      help="the nodata value of raster")

    parser.add_option("-o", "--output_shp",
                      action="store", dest="output_shp",
                      help="the output shapefile")

    parser.add_option("-m", "--output_mask",
                      action="store", dest="output_mask",
                      help="the output raster mask of valid regions, sometime, shapefile may have multiple polygons, "
                           "it will be convenient to use the raster mask ")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
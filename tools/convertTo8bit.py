#!/usr/bin/env python
# Filename: convertTo8bit 
"""
introduction: convert images to 8bit.
We also can use gdal_translate to do this, but gdal_translate has a error:
gdal_transform: scale does not conform to prescribed range https://github.com/OSGeo/gdal/issues/1813
so we set dst_min as 1, but in the output, still some 0 pixel (not nodata in the input).
In the future, GDAL may introduce -scaleclip to solve this problem.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 February, 2021
"""

import os, sys
from optparse import OptionParser


deeplabRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabRS)
import raster_io
import basic_src.io_function as io_function




def main(options, args):
    img_path = args[0]
    io_function.is_file_exist(img_path)

    img_np_allbands, src_nodata = raster_io.read_raster_all_bands_np(img_path)

    scales = options.scale
    print('input scale (src_min src_max dst_min dst_max): '+ str(scales))
    img_array_8bit = raster_io.image_numpy_allBands_to_8bit(img_np_allbands,scales,src_nodata=src_nodata,dst_nodata=options.dst_nodata)

    # save to file
    save_path = args[1]
    if options.dst_nodata is None:
        nodata = src_nodata
    else:
        nodata = options.dst_nodata
    return raster_io.save_numpy_array_to_rasterfile(img_array_8bit,save_path,img_path,
                                             nodata =nodata,compress='lzw',tiled='yes',bigtiff='if_safer')


if __name__ == '__main__':
    usage = "usage: %prog [options] image_path save_path"
    parser = OptionParser(usage=usage, version="1.0 2021-02-07")
    parser.description = 'convert images to unsigned 8bit '

    parser.add_option("-s", "--scale", nargs=4, action='append',
                      help="similar to src_min src_max dst_min dst_max, "
                           "repeat -s for multiple time for different bands if need")
    parser.add_option("-n", "--dst_nodata", action='store', type=int,
                      help="nodata for output, if None, then is the same as source file")

    parser.add_option("-u", "--hist_max_percent",
                      action="store", dest="hist_max_percent",type=float,
                      help="the upper percent for choosing the max pixel value")
    #
    parser.add_option("-l", "--hist_min_percent",
                      action="store", dest="hist_min_percent",type=float,
                      help="the lower percent for choosing the max pixel value")

    (options, args) = parser.parse_args()
    # print(options.scale)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

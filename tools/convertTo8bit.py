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
    save_path = args[1]
    io_function.is_file_exist(img_path)
    if os.path.isfile(save_path):
        print('%s exists, remove it if want to re-generate it'%save_path)
        return

    img_np_allbands, src_nodata = raster_io.read_raster_all_bands_np(img_path)
    if options.src_nodata is not None:
        src_nodata = options.src_nodata

    scales = options.scale
    if scales is not None:
        print('input scale (src_min src_max dst_min dst_max): '+ str(scales))
        img_array_8bit = raster_io.image_numpy_allBands_to_8bit(img_np_allbands,scales,src_nodata=src_nodata,dst_nodata=options.dst_nodata)
    else:
        min_percent = options.hist_min_percent
        max_percent = options.hist_max_percent
        min_max_value = options.min_max_value
        hist_bin_count = options.hist_bin_count
        img_array_8bit = raster_io.image_numpy_allBands_to_8bit_hist(img_np_allbands,min_max_value,per_min=min_percent,
                                                                     per_max=max_percent, bin_count=hist_bin_count, src_nodata=src_nodata,
                                                                     dst_nodata=options.dst_nodata)

    # save to file
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

    parser.add_option("-N", "--src_nodata", action='store', type=float,
                      help="nodata for the source data (different from the one in metadata), if None")

    parser.add_option("-n", "--dst_nodata", action='store', type=int,
                      help="nodata for output, if None, then is the same as source file")

    parser.add_option("-u", "--hist_max_percent",
                      action="store", dest="hist_max_percent",type=float,
                      help="the upper percent for choosing the max pixel value")
    #
    parser.add_option("-l", "--hist_min_percent",
                      action="store", dest="hist_min_percent",type=float,
                      help="the lower percent for choosing the max pixel value")

    parser.add_option("-c", "--hist_bin_count",
                      action="store", dest="hist_bin_count",type=int,default=10000,
                      help="the bin count for histogram when set hist_min_percent and hist_max_percent; for SAR image (float), should set > 10000")

    parser.add_option("-m", "--min_max_value",nargs=2, action="append",
                      dest="min_max_value",type=float,
                      help="if the value from hist_(min)_max_percent (less) greater than this one, it will be set as this"
                           "--min_max_value min max. (repeat to give multiple ones)")
    #


    (options, args) = parser.parse_args()
    # print(options.scale)
    # print(options.min_max_value)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

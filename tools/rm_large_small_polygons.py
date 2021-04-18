#!/usr/bin/env python
# Filename: rm_large_small_polygons 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 April, 2021
"""

import os,sys
from optparse import OptionParser
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd

import basic_src.io_function as io_function

def main(options, args):

    input_shp = args[0]
    io_function.is_file_exist(input_shp)

    min_area = options.min_area
    max_area = options.max_area
    output = os.path.basename(io_function.get_name_by_adding_tail(input_shp, 'MinMaxArea'))
    vector_gpd.remove_polygons_not_in_range(input_shp,'poly_area', min_area, max_area,output)

if __name__ == "__main__":
    usage = "usage: %prog [options] input_shp "
    parser = OptionParser(usage=usage, version="1.0 2021-2-21")
    parser.description = 'Introduction: segment subsidence based on DEM difference  '

    parser.add_option("-l", "--min_area",
                      action="store", dest="min_area", type=float, default=40,
                      help="the minimum area for each polygon")

    parser.add_option("-u", "--max_area",
                      action="store", dest="max_area", type=float, default=100000000,  # 10 km by 10 km
                      help="the maximum area for each polygon")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

#!/usr/bin/env python
# Filename: select_grids_data_overlap_vectors.py 
"""
introduction: select grid data (DEM difference raster, polygons, and others) if they overlap input polygons or points

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 December, 2023
"""

import os, sys
from optparse import OptionParser
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import pandas as pd
import vector_gpd


def find_grids_overlap_vector_shp():
    pass

def main(options, args):
    grid_indexes_shp= args[0]
    vector_shp_list = args[1:]

   # read the grid index file


   # go through each vector shp file
   for vector_shp in vector_shp_list:
       pass
   


if __name__ == '__main__':
    usage = "usage: %prog [options] grid_indexes_shp vector_shp1 vector_shp2 ... "
    parser = OptionParser(usage=usage, version="1.0 2023-12-18")
    parser.description = 'Introduction: select grid data (DEM diff, polygons etc) if they overlap input polygons or points '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='select_grid_data',
                      help="the folder to save DEMs")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

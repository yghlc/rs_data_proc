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


def find_grids_overlap_vector_shp(grid_polys, grid_attributes,vector_shp):

    overlap_grids = []
    cell_ids = []
    fileurl_list = []
    shapes = vector_gpd.read_polygons_gpd(vector_shp,b_fix_invalid_polygon=False)


    pass

def main(options, args):
    grid_indexes_shp= args[0]
    vector_shp_list = args[1:]
    if len(vector_shp_list) < 1:
        raise ValueError("No input vector files")

    #   cell_id (Integer64) = 58601
    #   tarball (String) = dem_diffs_2m_grid58601.tar.gz
    #   fileurl (String) = elevation-differences/ext11/dem_diffs_2m_grid58601.tar.gz
    #   POLYGON ((3320000 -120000,3340000 -120000,3340000 -140000,3320000 -140000,3320000 -120000))
    # read the grid index file
    grid_polys, grid_attributes = vector_gpd.read_polygons_attributes_list(grid_indexes_shp,['cell_id','fileurl'],
                                                                           b_fix_invalid_polygon=False)

    # go through each vector shp file
    for vector_shp in vector_shp_list:
        print(vector_shp)
   


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

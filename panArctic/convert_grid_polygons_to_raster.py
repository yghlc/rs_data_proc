#!/usr/bin/env python
# Filename: create_grid_polygons_perma 
"""
introduction: based on coverage of ArcticDEM, create 20 by 20 km grid

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 March, 2021
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import raster_io

import pandas as pd

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp')

def main():
    # save
    ArcticDEM_grid_20km = os.path.join(shp_dir,'ArcticDEM_grid_20km.shp')

    # ref_raster='grid_20km_bin.tif'
    polygons,ids = vector_gpd.read_polygons_attributes_list(ArcticDEM_grid_20km,'id',b_fix_invalid_polygon=False)
    save_raster = os.path.join(shp_dir, 'grid_20km_id.tif')
    # raster_io.burn_polygons_to_a_raster(ref_raster,polygons,ids,save_raster,date_type='uint16')


    # if no reference raster
    extent = vector_gpd.get_vector_file_bounding_box(ArcticDEM_grid_20km)
    # print(extent)
    res = 20000  # 20 km
    wkt_string = map_projection.get_raster_or_vector_srs_info_proj4(ArcticDEM_grid_20km)

    raster_io.burn_polygons_to_a_raster(None,polygons,ids,save_raster,date_type='uint16',
                                        xres=res,yres=res,extent=extent,ref_prj=wkt_string)



    pass

if __name__ == "__main__":
    main()
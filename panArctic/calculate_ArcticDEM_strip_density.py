#!/usr/bin/env python
# Filename: calculate_ArcticDEM_strip_density.py 
"""
introduction: calculate density of ArcticDEM strip, at resolution of 100 m by 100 m

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 17 August, 2022
"""

import os,sys
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import raster_io

import time
from datetime import datetime

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7')

def cal_ArcticDEM_strip_density_v1():
    t0 = time.time()
    ArcticDEM_strip_index = os.path.join(shp_dir,'ArcticDEM_Strip_Index_Rel7.shp')

    # ref_raster='grid_20km_bin.tif'
    # polygons = vector_gpd.read_polygons_attributes_list(ArcticDEM_grid_20km,'id',b_fix_invalid_polygon=False)
    polygons = vector_gpd.read_polygons_gpd(ArcticDEM_strip_index,b_fix_invalid_polygon=False)
    print('read ArcticDEM strip polygons, cost: ', time.time() -t0, ' seconds')

    # res = 100  # 100 m, time consuming,may take a few months
    res = 1000  # may take a few hours to one day


    # save_raster = os.path.join(shp_dir, 'ArcticDEM_strip_density_100m.tif')
    save_raster = 'ArcticDEM_strip_density_%dm.tif'%res

    t0 = time.time()
    # if no reference raster
    extent = vector_gpd.get_vector_file_bounding_box(ArcticDEM_strip_index)
    # print(extent)
    wkt_string = map_projection.get_raster_or_vector_srs_info_proj4(ArcticDEM_strip_index)
    nodata = 0

    ext_tif = 'ArcticDEM_strip_extent_%dm.tif'%res
    if os.path.isfile(ext_tif) is False:
        raster_io.burn_polygons_to_a_raster(None, polygons, 1, ext_tif, date_type='uint16',
                                                  xres=res, yres=res, extent=extent, ref_prj=wkt_string, nodata=nodata)
    else:
        print('%s already exists, skip'% ext_tif)
    print('produce extent tif, cost: ', time.time() - t0, ' seconds')
    t0 = time.time()

    density_sum = None
    for idx, poly in enumerate(polygons):
        print(datetime.now(), 'burning %d th polygon, total %d:'%(idx,len(polygons)))
        # this is very time consuming, each polygon takes around 50 seconds, on my Mac, total: 206532
        if idx == 0:
            density = raster_io.burn_polygons_to_a_raster(None,[poly],1,None,date_type='uint16',
                                            xres=res,yres=res,extent=extent,ref_prj=wkt_string, nodata=nodata)
            density_sum = density
        else:
            density = raster_io.burn_polygons_to_a_raster(None, [poly], 1, None, date_type='uint16',
                                                          xres=res, yres=res, extent=extent, ref_prj=wkt_string,
                                                          nodata=nodata)
            density_sum += density

        # # test
        # if idx > 2:
        #     break
    print('calculate density at res: %d m, cost: '%res, time.time() - t0, ' seconds')
    # t0 = time.time()

    # save to file
    raster_io.save_numpy_array_to_rasterfile(density_sum,save_raster,ext_tif)


def cal_ArcticDEM_strip_density_v2():
    t0 = time.time()
    ArcticDEM_strip_index = os.path.join(shp_dir, 'ArcticDEM_Strip_Index_Rel7.shp')

    polygons = vector_gpd.read_polygons_gpd(ArcticDEM_strip_index, b_fix_invalid_polygon=False)
    print('read ArcticDEM strip polygons, cost: ', time.time() - t0, ' seconds')
    # save_raster = os.path.join(shp_dir, 'ArcticDEM_strip_density_100m.tif')
    save_raster = 'ArcticDEM_strip_density_100m.tif'

    t0 = time.time()
    # if no reference raster
    extent = vector_gpd.get_vector_file_bounding_box(ArcticDEM_strip_index)
    # print(extent)
    res = 100  # 100 m
    wkt_string = map_projection.get_raster_or_vector_srs_info_proj4(ArcticDEM_strip_index)
    nodata = 2 ** 16 - 1

    density_sum = None
    # for idx, poly in enumerate(polygons):
    #     print(datetime.now(), 'burning %d th polygon, total %d:' % (idx, len(polygons)))
    #     # this is very time consuming, each polygon takes around 50 seconds, on my Mac, total: 206532
    #     if idx == 0:
    #         density = raster_io.burn_polygons_to_a_raster(None, [poly], 1, None, date_type='uint16',
    #                                                       xres=res, yres=res, extent=extent, ref_prj=wkt_string,
    #                                                       nodata=nodata)
    #         density_sum = density
    #     else:
    #         density = raster_io.burn_polygons_to_a_raster(None, [poly], 1, None, date_type='uint16',
    #                                                       xres=res, yres=res, extent=extent, ref_prj=wkt_string,
    #                                                       nodata=nodata)
    #         density_sum += density
    #
    #     # test
    #     if idx > 10:
    #         break
    #
    # # save to file
    # raster_io.save_numpy_array_to_rasterfile(density_sum, save_raster, ext_tif)





if __name__ == "__main__":
    cal_ArcticDEM_strip_density_v1()
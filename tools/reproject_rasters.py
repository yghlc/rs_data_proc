#!/usr/bin/env python
# Filename: reproject_rasters.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 07 March, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('/')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.RSImageProcess as RSImageProcess
import basic_src.map_projection as map_projection

def crop_mosaic_reproject_dem_diff(grid_dem_tifs, pre_name, extent_poly, o_res, new_prj, b_mosaic=False):

    # crop
    crop_tif_list = []
    # crop to the same extent
    crop_tif_dir = os.path.join('dem_diff_crop')
    if os.path.isdir(crop_tif_dir) is False:
        io_function.mkdir(crop_tif_dir)
    for tif in grid_dem_tifs:
        save_crop_path = os.path.join(crop_tif_dir, os.path.basename(io_function.get_name_by_adding_tail(tif, 'sub')) )
        if os.path.isfile(save_crop_path):
            basic.outputlogMessage('%s exists, skip cropping' % save_crop_path)
            crop_tif_list.append(save_crop_path)
        else:
            crop_tif = subset_image_by_polygon_box(tif, save_crop_path, extent_poly, resample_m='near',
                            o_format='VRT', out_res=o_res,same_extent=True,thread_num=2)
            if crop_tif is False:
                raise ValueError('warning, crop %s failed' % tif)
            crop_tif_list.append(crop_tif)
    grid_dem_tifs = crop_tif_list

    # mosaic

    if b_mosaic:
        save_mosaic = pre_name + '_DEM_diff.tif'
        resample_method = 'average'
        # create mosaic, can handle only input one file, but is slow
        result = RSImageProcess.mosaic_crop_images_gdalwarp(grid_dem_tifs, save_mosaic, resampling_method=resample_method,
                                                            o_format='Gtiff',
                                                            compress='lzw', tiled='yes', bigtiff='if_safer',
                                                            thread_num=2)
        if result is False:
            sys.exit(1)

        grid_dem_tifs = [save_mosaic]

    # reproject
    for tif in grid_dem_tifs:
        t_file = io_function.get_name_by_adding_tail(tif,'prj')
        map_projection.transforms_raster_srs(tif,new_prj,t_file,o_res,o_res,resample_m='bilinear',
                              o_format='GTiff',compress='lzw', tiled='yes', bigtiff='if_safer')


    pass


# crop, mosaic, and reproject if necessary
if extent_shp_or_ids_txt.endswith('.shp'):
    pre_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]

    # extent polygons and projection (proj4)
    extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp_or_ids_txt)
    grid_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)

    if extent_shp_prj != grid_shp_prj:
        basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
                               % (extent_shp_or_ids_txt, grid_20_shp, os.path.basename(extent_shp_or_ids_txt)))
        epsg = map_projection.get_raster_or_vector_srs_info_epsg(grid_20_shp)
        extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp_or_ids_txt, epsg)
    else:
        extent_polys = vector_gpd.read_polygons_gpd(extent_shp_or_ids_txt)

    crop_mosaic_reproject_dem_diff(grid_dem_tifs, pre_name, extent_poly, o_res, new_prj, b_mosaic=False)


def main(options, args):
    pass


if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp or grid_id_list.txt "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: re-project, crop, and get mosaic  '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    parser.add_option("-o", "--out_res",
                      action="store", dest="out_res",type=float,default=2.0,
                      help="the resolution for final output")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
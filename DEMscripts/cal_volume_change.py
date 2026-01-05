#!/usr/bin/env python
# Filename: cal_volume_change.py 
"""
introduction: calculate the volumetric change from elevation difference, or estimate the volume from single DEM

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 04 January, 2026
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function

from datetime import datetime
import raster_statistic

import numpy as np
from multiprocessing import Pool

def cal_volume(in_array, nodata,res_x, res_y,range=None, scale=1.0):
    # print(in_array)
    # print(in_array.shape)
    data_1d = in_array.flatten()
    # print(data_1d)
    total_size =  data_1d.size
    data_1d = data_1d[ data_1d != nodata]
    size_2 = data_1d.size
    data_1d = data_1d[~np.isnan(data_1d)]  # remove nan value
    size_3 = data_1d.size

    nodata_count = total_size - size_2
    nan_count = size_2 - size_3

    # apply the scale before checking range
    data_1d = data_1d * scale

    if range is not None:
        lower = range[0]
        upper = range[1]
        if lower is None:
            data_1d = data_1d[data_1d <= upper]
        elif upper is None:
            data_1d = data_1d[data_1d >= lower]
        else:
            data_1d = data_1d[np.logical_and( data_1d >= lower, data_1d <= upper ) ]
    size_4 = data_1d.size
    valid_count = total_size - nodata_count - nan_count
    not_in_range_c = size_3 - size_4

    if data_1d.size < 1:
        return 0, 0, 0, 0, 0

    # print(data_1d)
    volume = np.sum(np.abs(data_1d)*res_x*res_y)

    # print(nodata)
    # print(total_size,size_2,size_3,size_4)
    # print(valid_count, not_in_range_c)

    return volume, valid_count, not_in_range_c, np.min(data_1d), np.max(data_1d)

def cal_volume_one_polygon(idx, polygon, res,image_tiles, img_tile_polygons,nodata=None,range=range,
                           all_touched=True,tile_min_overlap=None,scale=1.0):
    ''' calculate volume for a polygon, allow multiple polygons
    '''
    # read the raster
    out_image = raster_statistic.read_raster_data_img_list(idx, polygon, image_tiles, img_tile_polygons, nodata=nodata,
                                          band=1, all_touched=all_touched, tile_min_overlap=tile_min_overlap)

    # calculate the
    return cal_volume(out_image, nodata,  res, res, range=range,scale=scale)

def cal_volumetric_change_from_dem_diff(boundary_shp,raster_file_or_files,save_path=None,range=None,
                                        tile_min_overlap=None, nodata=None,buffer=None, all_touched=True, process_num=1):
    '''
    calculate the volumetric changes of all polygons in a vector file, allow multiple rasters
    :param boundary_shp:
    :param raster_file_or_files:
    :param save_path:
    :param valid_thr:
    :param tile_min_overlap:
    :param nodata: incorrect setting of nodata will end in low values of "INrangePer".
    :param buffer:
    :param all_touched:
    :param process_num:
    :return:
    '''

    io_function.is_file_exist(boundary_shp)
    if isinstance(raster_file_or_files,str):
        io_function.is_file_exist(raster_file_or_files)
        image_tiles = [raster_file_or_files]
    elif isinstance(raster_file_or_files,list):
        image_tiles = raster_file_or_files
    else:
        raise ValueError('unsupport type for %s'%str(raster_file_or_files))
    # check projection (assume we have the same projection), check them outside this function

    # get image box
    img_tile_boxes = [raster_io.get_image_bound_box(tile) for tile in image_tiles]
    img_tile_polygons = [vector_gpd.convert_image_bound_to_shapely_polygon(box) for box in img_tile_boxes]
    polygons = vector_gpd.read_polygons_gpd(boundary_shp)
    if len(polygons) < 1:
        basic.outputlogMessage('No polygons in %s' % boundary_shp)
        return False
    # polygons_json = [mapping(item) for item in polygons]  # no need when use new verion of rasterio
    if buffer is not None:
        polygons = [poly.buffer(buffer) for poly in polygons]

    res_x, res_y = raster_io.get_xres_yres_file(image_tiles[0])
    heigth, width, bandcount, dtype = raster_io.get_height_width_bandnum_dtype(image_tiles[0])
    scale = 1.0
    if dtype == 'int16':
        scale = 0.01    # DEM diff saved in 16 bit, that is cm, multiply by 0.01 to m.

    volume_res_list = []
    in_range_per = []
    dem_diff_min = []
    dem_diff_max = []

    # process polygons one by one and the corresponding image tiles (parallel and save memory)
    # also to avoid error: daemonic processes are not allowed to have children
    if process_num == 1:

        for idx, polygon in enumerate(polygons):
            out_volume, pixel_count, not_in_range_c, min_v, max_v = cal_volume_one_polygon(idx, polygon, res_x, image_tiles, img_tile_polygons,
                                                nodata=nodata, range=range,all_touched=all_touched,
                                                tile_min_overlap=tile_min_overlap,scale=scale)
            volume_res_list.append(out_volume)
            if pixel_count > 0:
                in_range_per.append(100*(pixel_count-not_in_range_c)/pixel_count)
            else:
                in_range_per.append(0)
            dem_diff_min.append(min_v)
            dem_diff_max.append(max_v)

    elif process_num > 1:
        threadpool = Pool(process_num)
        para_list = [(idx, polygon, res_x, image_tiles, img_tile_polygons, nodata, range, all_touched,
                      tile_min_overlap,scale)
                     for idx, polygon in enumerate(polygons)]
        results = threadpool.starmap(cal_volume_one_polygon, para_list)
        # save results
        for res in results:
            out_volume, pixel_count, not_in_range_c, min_v, max_v = res
            volume_res_list.append(out_volume)
            if pixel_count > 0:
                in_range_per.append(100 * (pixel_count - not_in_range_c) / pixel_count)
            else:
                in_range_per.append(0)
            dem_diff_min.append(min_v)
            dem_diff_max.append(max_v)
        threadpool.close()
    else:
        raise ValueError('Wrong process number: %s ' % str(process_num))


    # if os.path.isfile(save_path) is False:
    #     io_function.copy_file_to_dst(boundary_shp, save_path)

    vol_att = {'vol_chg': volume_res_list, 'INrangePer': in_range_per, 'D_diff_min':dem_diff_min, 'D_diff_max':dem_diff_max}
    save_format = vector_gpd.guess_file_format_extension(save_path)
    if save_path != boundary_shp:
        vector_gpd.add_attributes_to_shp(boundary_shp, vol_att, save_as=save_path, format=save_format)
    else:
        vector_gpd.add_attributes_to_shp(boundary_shp,vol_att,format=save_format)


def test_cal_volumetric_change_from_dem_diff():
    print(datetime.now(),'Testing:', 'test_cal_volumetric_change_from_dem_diff')
    data_dir = os.path.expanduser('~/Data/rts_ArcticDEM_mapping/cal_volume_from_demDiff/')
    # poly_shp = os.path.join(data_dir,'dem_diff_sam_results_select','dem_diff_sam_selected_by_quickshift_ARTS_v6.gpkg')
    poly_shp = os.path.join(data_dir,'input_test.shp')
    dem_diff = os.path.join(data_dir,'valid_json_files_20251222_TP_dem_diff_2023_mosaic.vrt')

    # nodata = 0
    nodata = 32767

    save_path = os.path.join(data_dir,'testing_output.shp')
    cal_volumetric_change_from_dem_diff(poly_shp,dem_diff,save_path=save_path, range=[None,-2],nodata=nodata, process_num=8)


def main(options, args):

    # test_cal_volumetric_change_from_dem_diff()
    # sys.exit(0)

    boundary_shp = args[0]
    dem_diff_path = args[1]
    save_path = boundary_shp if options.save_path is None else options.save_path
    process_num = options.process_num
    # nodata on DEM difference
    bit16_nodata = 32767


    cal_volumetric_change_from_dem_diff(boundary_shp, dem_diff_path, save_path=save_path, range=[None,-2.0],
                                        tile_min_overlap=None, nodata=bit16_nodata, all_touched=True,
                                        process_num=process_num)




if __name__ == "__main__":
    usage = "usage: %prog [options] boundary_shp DEM_diff_path_or_dir "
    parser = OptionParser(usage=usage, version="1.0 2026-01-04")
    parser.description = 'Introduction: calculate volume change   '

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="The path to save the output as a new file")

    parser.add_option("-p", "--process_num",
                      action="store", dest="process_num", type=int, default=8,
                      help="number of processes to create the mosaic")




    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

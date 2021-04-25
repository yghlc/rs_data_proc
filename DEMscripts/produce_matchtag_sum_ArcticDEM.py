#!/usr/bin/env python
# Filename: produce_matchtag_sum_ArcticDEM.py
"""
introduction: produce sum of matchtag of ArcticDEM, corresponding to each DEM difference

# I divide the coverage of ArcticDEM to 20 km by 20 km grids:
~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp, each grid have a ID, in projection fo EPSG:3413, Polar Stereographic North.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 24 April, 2021
"""

import os,sys
from optparse import OptionParser
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import basic_src.io_function as io_function

import numpy as np

from dem_mosaic_crop import mosaic_crop_dem
from dem_mosaic_crop import get_dem_tif_ext_polygons
from dem_mosaic_crop import subset_image_by_polygon_box
from dem_difference import dem_diff_newest_oldest

from produce_DEM_diff_ArcticDEM import get_grid_20

#
from dem_common import grid_20_shp, arcticDEM_reg_tif_dir, grid_matchtag_sum_dir

def get_existing_matchtag_sum(matchtag_sum_dir, grid_base_name, grid_ids):

    existing_tif = []
    grid_id_no_matchtag_sum_tiff = []
    for id in grid_ids:
        # dem_diff = os.path.join(dem_diff_dir, grid_base_name + '_DEM_diff_grid%d.tif'%id)
        # if os.path.isfile(dem_diff):
        #     existing_tif.append(dem_diff)
        #     continue

        matchtag_sum_files = io_function.get_file_list_by_pattern(matchtag_sum_dir, '*_matchtag_sum_grid%d.tif'%id)
        if len(matchtag_sum_files) == 1:
            existing_tif.append(matchtag_sum_files[0])
            continue
        elif len(matchtag_sum_files) > 1:
            existing_tif.append(matchtag_sum_files[0])
            basic.outputlogMessage('warning, There are multiple matchtag_sum tif for grid: %d'%id)
            for item in matchtag_sum_files: basic.outputlogMessage(item)
            continue
        else:
            pass

        grid_id_no_matchtag_sum_tiff.append(id)
    if len(existing_tif) > 0:
        basic.outputlogMessage('%d existing grid matchtag_sum files for the input grid_ids or extent'%len(existing_tif))
    else:
        basic.outputlogMessage('no existing grid matchtag_sum files')
    return existing_tif, grid_id_no_matchtag_sum_tiff

def sum_matchtag(input_tifs, save_path):

    if len(input_tifs) < 1:
        return False
    # check band, with, height
    height, width, count, dtype = raster_io.get_height_width_bandnum_dtype(input_tifs[0])
    for idx in range(1, len(input_tifs)):
        h, w, c, type = raster_io.get_height_width_bandnum_dtype(input_tifs[idx])
        if h!=height or w!=width or c!=count or type!=dtype:
            raise ValueError('size or data type is different between %s and %s'%(input_tifs[0], input_tifs[idx]))

    if count != 1:
        raise ValueError('Matchtag should only have one band')

    sum_data = np.zeros((height,width),dtype=np.uint8)
    for tif in input_tifs:
        data = raster_io.read_raster_one_band_np(tif)
        sum_data += data

    # save to file
    raster_io.save_numpy_array_to_rasterfile(sum_data,save_path,input_tifs[0],compress='lzw',tiled='yes',bigtiff='if_safer')
    return True


def produce_matchtag_sum_grids(grid_polys, grid_ids, pre_name, matchtag_tifs,o_res,process_num=4):

    dem_ext_polys = get_dem_tif_ext_polygons(matchtag_tifs)
    matchtag_sum_tifs = []
    # mosaic and crop
    for grid_id, grid_poly in zip(grid_ids, grid_polys):

        save_dir = 'grid_%d_tmp_files'%grid_id

        # check free disk space
        work_dir = './'
        free_GB = io_function.get_free_disk_space_GB(work_dir)
        total_wait_time = 0
        while free_GB < 50 and total_wait_time < 60*60*12:
            print(' The free disk space (%.4f) is less than 50 GB, wait 60 seconds'%free_GB)
            time.sleep(60)
            total_wait_time += 60
            free_GB = io_function.get_free_disk_space_GB(work_dir)


        # get subset of tifs
        dem_poly_index = vector_gpd.get_poly_index_within_extent(dem_ext_polys, grid_poly)
        if len(dem_poly_index) < 1:
            basic.outputlogMessage('warning, no dem tifs within %d grid, skip' % grid_id)
            continue
        dem_list_sub = [matchtag_tifs[index] for index in dem_poly_index]

        mosaic_tif_list = mosaic_crop_dem(dem_list_sub, save_dir, grid_id, grid_poly, False, False,
                        process_num, 0, o_res, pre_name, resample_method='average')

        # sum matchtag
        save_matchtag_sum = os.path.join(grid_matchtag_sum_dir, pre_name + '_count%d_'%len(mosaic_tif_list)
                                         + '_matchtag_sum_grid%d.tif'%grid_id)

        if sum_matchtag(mosaic_tif_list, save_matchtag_sum):
            matchtag_sum_tifs.append(save_matchtag_sum)

    return matchtag_sum_tifs



def main(options, args):
    extent_shp_or_ids_txt = args[0]
    process_num = options.process_num
    o_res = options.out_res

    if os.path.isdir(grid_matchtag_sum_dir) is False:
        io_function.mkdir(grid_matchtag_sum_dir)

    basic.setlogfile('produce_matchtag_sum_ArcticDEM_log_%s.txt'%timeTools.get_now_time_str())

    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    grid_base_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]
    grid_polys, grid_ids = get_grid_20(extent_shp_or_ids_txt,all_grid_polys, all_ids)

    # check dem difference existence
    grid_dem_tifs, grid_ids_no_sum = get_existing_matchtag_sum(grid_matchtag_sum_dir,grid_base_name,grid_ids)
    if len(grid_ids_no_sum) > 0:
        # refine grid_polys
        if len(grid_ids) > len(grid_ids_no_sum):
            id_index = [grid_ids.index(id) for id in grid_ids_no_sum]
            grid_polys = [grid_polys[idx] for idx in id_index]

        # # download ArcticDEM and applying registration
        # tarballs, reg_tifs = download_dem_tarball(dem_strip_shp, grid_polys, arcticDEM_tarball_dir, grid_base_name,
        #                                         reg_tif_dir=arcticDEM_reg_tif_dir, poly_ids=grid_ids_no_demDiff)
        #
        # # unpack and applying registration
        # if len(tarballs) > 0:
        #     basic.outputlogMessage('Processs %d dem tarballs'%len(tarballs))
        #     out_reg_tifs = process_dem_tarball(tarballs,'./',arcticDEM_reg_tif_dir,remove_inter_data=True, apply_registration=True)
        #     basic.outputlogMessage('Get %d new registration dem tifs' % len(out_reg_tifs))
        #     reg_tifs.extend(out_reg_tifs)

        reg_tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir,bsub_folder=False)
        matchtag_tifs = [tif for tif in reg_tifs if 'matchtag' in tif]  # only keep matchtag
        # crop, sum
        out_dem_diffs = produce_matchtag_sum_grids(grid_polys, grid_ids_no_sum, grid_base_name, matchtag_tifs,o_res,
                                                   process_num=process_num)



if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp or grid_id_list.txt "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: produce matchtag sum from multiple temporal ArcticDEM  '

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

#!/usr/bin/env python
# Filename: dem_headwall_extraction_grid.py 
"""
introduction: extract headwall from slope, but process the multi-temporal DEM grid by grid (similar to the way producing DEM difference)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 14 May, 2021
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
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

from multiprocessing import Pool

# some parameters
b_mosaic_id = True          # mosaic dem with the same id
b_mosaic_date = True        # mosaic dem within one day
b_apply_matchtag = False     # don't apply matchtag, it seems that matchtag make slope worse?

# parameters for extracting headwall
min_slope = 20
min_size = 200
max_size = 50000
max_axis_width = 80
max_box_WH = 600

from produce_DEM_diff_ArcticDEM import get_grid_20
from dem_mosaic_crop import get_dem_tif_ext_polygons
from dem_mosaic_crop import mosaic_crop_dem
from dem_to_hillshade_slope_8bit import dem_to_slope
from dem_headwall_extraction import extract_headwall_from_slope

from dem_common import grid_20_shp, arcticDEM_reg_tif_dir, grid_dem_headwall_shp_dir

def get_existing_grid_headwall_shp(headwall_shp_dir, grid_base_name, grid_ids):

    existing_grid_headwall_shp = []
    grid_id_no_headwall_shp = []
    for id in grid_ids:

        headwall_shps = io_function.get_file_list_by_pattern(headwall_shp_dir, '*_grid%d.shp'%id)
        if len(headwall_shps) == 1:
            existing_grid_headwall_shp.append(headwall_shps[0])
            continue
        elif len(headwall_shps) > 1:
            existing_grid_headwall_shp.append(headwall_shps[0])
            basic.outputlogMessage('warning, There are multiple headwall shps for grid: %d'%id)
            for item in headwall_shps: basic.outputlogMessage(item)
            continue
        else:
            pass

        grid_id_no_headwall_shp.append(id)
    if len(existing_grid_headwall_shp) > 0:
        basic.outputlogMessage('%d existing grid headwall shps for the input grid_ids or extent'%len(existing_grid_headwall_shp))
    else:
        basic.outputlogMessage('no existing grid headwall shps')
    return existing_grid_headwall_shp, grid_id_no_headwall_shp

def one_dem_to_slope(tif, slope_tif_dir):
    output = os.path.join(slope_tif_dir, os.path.basename(io_function.get_name_by_adding_tail(tif, 'slope')))
    slope_tif = dem_to_slope(tif, output, '')
    return slope_tif

def dem_list_to_slope_list(dem_list, save_dir, extent_id, process_num=1):

    slope_list = []
    slope_tif_dir = os.path.join(save_dir, 'slope_sub_%d' % extent_id)
    if os.path.isdir(slope_tif_dir) is False:
        io_function.mkdir(slope_tif_dir)

    if process_num == 1:
        for idx, tif in enumerate(dem_list):
            slope_tif = one_dem_to_slope(tif,slope_tif_dir)
            if slope_tif is not False:
                slope_list.append(slope_tif)
    elif process_num > 1:

        # change backk to multi process of gdalwarp, when get mosaic, gdalwarp multi-thread cannot fully utlized CPUs
        theadPool = Pool(process_num)  # multi processes

        parameters_list = [(tif, slope_tif_dir) for idx, tif in enumerate(dem_list)]
        results = theadPool.starmap(one_dem_to_slope, parameters_list)  # need python3
        slope_list = [ out for out in results if out is not False]
    else:
        raise ValueError('Wrong process number: %s'%str(process_num))

    return slope_list

def merge_multi_headwall_shp():
    pass

def extract_headwall_grids(grid_polys, grid_ids, pre_name,reg_tifs,b_mosaic_id,
                           b_mosaic_date,keep_dem_percent, o_res,process_num=1):

    dem_ext_polys = get_dem_tif_ext_polygons(reg_tifs)
    headwall_shp_folders = []
    # mosaic and crop
    for grid_id, grid_poly in zip(grid_ids, grid_polys):

        save_dir = 'grid_%d_tmp_files' % grid_id

        # check free disk space
        work_dir = './'
        free_GB = io_function.get_free_disk_space_GB(work_dir)
        total_wait_time = 0
        while free_GB < 50 and total_wait_time < 60 * 60 * 12:
            print(' The free disk space (%.4f) is less than 50 GB, wait 60 seconds' % free_GB)
            time.sleep(60)
            total_wait_time += 60
            free_GB = io_function.get_free_disk_space_GB(work_dir)

        # get subset of tifs
        dem_poly_index = vector_gpd.get_poly_index_within_extent(dem_ext_polys, grid_poly)
        if len(dem_poly_index) < 1:
            basic.outputlogMessage('warning, no dem tifs within %d grid, skip' % grid_id)
            continue
        dem_list_sub = [reg_tifs[index] for index in dem_poly_index]


        mosaic_tif_list = mosaic_crop_dem(dem_list_sub, save_dir, grid_id, grid_poly, b_mosaic_id, b_mosaic_date,
                                          process_num, keep_dem_percent, o_res, pre_name, resample_method='average',
                                          b_mask_matchtag=b_apply_matchtag)

        # dem co-registration (cancel, the result in not good with the default setting)

        # to slope
        slope_tifs = dem_list_to_slope_list(mosaic_tif_list,save_dir,grid_id,process_num=process_num)

        # extract headwall
        multi_headwall_shp_dir = os.path.join(save_dir, 'headwall_shp_sub_%d' % grid_id)
        if os.path.isdir(multi_headwall_shp_dir) is False:
            io_function.mkdir(multi_headwall_shp_dir)
        for idx, slope in enumerate(slope_tifs):
            working_dir = os.path.join(save_dir,os.path.splitext(os.path.basename(slope))[0])
            if os.path.isdir(working_dir) is False:
                io_function.mkdir(working_dir)
            if extract_headwall_from_slope(idx, len(slope_tifs), slope, working_dir, multi_headwall_shp_dir, min_slope, min_size,
                                           max_size, max_axis_width, max_box_WH, process_num) is False:
                basic.outputlogMessage('extract headwall from %s failed'%slope)

        # merge headwall detected on different dates.


    return headwall_shp_folders



def main(options, args):
    extent_shp_or_ids_txt = args[0]
    process_num = options.process_num
    keep_dem_percent = options.keep_dem_percent
    o_res = options.out_res

    basic.setlogfile('produce_headwall_shp_ArcticDEM_log_%s.txt'%timeTools.get_now_time_str())

    if os.path.isdir(grid_dem_headwall_shp_dir) is False:
        io_function.mkdir(grid_dem_headwall_shp_dir)

        # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    grid_base_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]
    grid_polys, grid_ids = get_grid_20(extent_shp_or_ids_txt, all_grid_polys, all_ids)

    # check dem difference existence
    grid_headwall_shps, grid_id_no_headwall_shp = get_existing_grid_headwall_shp(grid_dem_headwall_shp_dir, grid_base_name, grid_ids)
    if len(grid_id_no_headwall_shp) > 0:
        # refine grid_polys
        if len(grid_ids) > len(grid_id_no_headwall_shp):
            id_index = [grid_ids.index(id) for id in grid_id_no_headwall_shp]
            grid_polys = [grid_polys[idx] for idx in id_index]


        reg_tifs = io_function.get_file_list_by_ext('.tif', arcticDEM_reg_tif_dir, bsub_folder=False)
        reg_tifs = [tif for tif in reg_tifs if 'matchtag' not in tif]  # remove matchtag
        # crop, mosacic, difference
        headwall_shp_folders = extract_headwall_grids(grid_polys, grid_id_no_headwall_shp, grid_base_name, reg_tifs,
                                               b_mosaic_id, b_mosaic_date,
                                               keep_dem_percent, o_res,process_num=process_num)




if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp or grid_id_list.txt "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: produce DEM difference from multiple temporal ArcticDEM  '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    parser.add_option("-p", "--keep_dem_percent",
                      action="store", dest="keep_dem_percent",type=float,default=10.0,
                      help="keep dem with valid percentage greater than this value")

    parser.add_option("-o", "--out_res",
                      action="store", dest="out_res",type=float,default=2.0,
                      help="the resolution for final output")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

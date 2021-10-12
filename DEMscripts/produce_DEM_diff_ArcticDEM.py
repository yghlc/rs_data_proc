#!/usr/bin/env python
# Filename: ArcticDEM_proc_grid.py
"""
introduction: produce elevation differences from  multi temporal ArcticDEM

# I divide the coverage of ArcticDEM to 20 km by 20 km grids:
~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp, each grid have a ID, in projection fo EPSG:3413, Polar Stereographic North.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 06 March, 2021
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
import basic_src.RSImageProcess as RSImageProcess



from download_arcticDEM import download_dem_tarball
from ArcticDEM_unpack_registration import process_dem_tarball
from dem_mosaic_crop import mosaic_crop_dem
from dem_mosaic_crop import get_dem_tif_ext_polygons
from dem_mosaic_crop import subset_image_by_polygon_box
from dem_difference import dem_diff_newest_oldest

# some parameters
b_mosaic_id = True
b_mosaic_date = True        # mosaic dem with one days
b_max_subsidence = False        # apply max_subsidence make results worse
b_apply_matchtag = True
b_mosaic_year = True        # mosaic dem for the same year, choose DEM close to July 1 on top.
b_mask_stripDEM_outlier = True  # mask outliers in strip DEM using the ArcticDEM tiles
b_mask_surface_water = True     # mask pixel of surface water


from dem_common import grid_20_shp,dem_strip_shp

# some folder paths
from dem_common import arcticDEM_reg_tif_dir,grid_dem_diffs_dir


def get_grid_20(extent_shp_or_id_txt, grid_polys, ids):
    '''
    get grid polygons and ids based on input extent (polygon in shpaefile) or ids (txt file)
    if "file_name_base+'_grid_ids.txt'" exists, it will read id in this file directly.
    :param extent_shp_or_id_txt:
    :param grid_polys:
    :param ids:
    :return:
    '''

    io_function.is_file_exist(extent_shp_or_id_txt)
    file_name_base = os.path.splitext(os.path.basename(extent_shp_or_id_txt))[0]

    if extent_shp_or_id_txt.endswith('.txt'):
        grid_ids = io_function.read_list_from_txt(extent_shp_or_id_txt)
        grid_ids = [int(item) for item in grid_ids ]
    else:
        shp_corresponding_grid_ids_txt = file_name_base+'_grid_ids.txt'
        if os.path.isfile(shp_corresponding_grid_ids_txt):
            print('corresponding grid ids txt file for %s exists, read grid id from txt'%extent_shp_or_id_txt)
            grid_ids = [ int(item) for item in io_function.read_list_from_txt(shp_corresponding_grid_ids_txt)]
            basic.outputlogMessage('read %d grids within the extents (%s)'
                                   % (len(grid_ids), os.path.basename(extent_shp_or_id_txt)))
        else:
            # extent polygons and projection (proj4)
            extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp_or_id_txt)
            if extent_shp_prj == '':
                raise ValueError('get proj4 of %s failed'%extent_shp_or_id_txt)
            grid_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)
            if grid_shp_prj=='':
                raise ValueError('get proj4 of %s failed' % grid_20_shp)

            if extent_shp_prj != grid_shp_prj:
                basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
                                       % (extent_shp_or_id_txt, grid_20_shp, os.path.basename(extent_shp_or_id_txt)))
                epsg = map_projection.get_raster_or_vector_srs_info_epsg(grid_20_shp)
                # print(epsg)
                # extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_shp_prj.strip())
                extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp_or_id_txt, epsg)
            else:
                extent_polys = vector_gpd.read_polygons_gpd(extent_shp_or_id_txt)

            ext_poly_count = len(extent_polys)
            if ext_poly_count < 1:
                raise ValueError('No polygons in %s'%extent_shp_or_id_txt)
            grid_index = []
            # if there are many polygons, this will take time.
            for idx,ext_poly in enumerate(extent_polys):
                print(timeTools.get_now_time_str(), idx, ext_poly_count)
                index = vector_gpd.get_poly_index_within_extent(grid_polys, ext_poly)
                grid_index.extend(index)
            grid_index = list(set(grid_index))  # remove duplicated ids
            basic.outputlogMessage('find %d grids within the extents (%s)' % (len(grid_index), os.path.basename(extent_shp_or_id_txt)) )

            grid_ids = [ ids[idx] for idx in grid_index]
            grid_ids_str = [str(item) for item in grid_ids ]
            io_function.save_list_to_txt(shp_corresponding_grid_ids_txt,grid_ids_str)

    id_index = [ids.index(id) for id in grid_ids]
    selected_grid_polys = [grid_polys[idx] for idx in id_index ]

    return selected_grid_polys, grid_ids

def get_existing_dem_diff(dem_diff_dir, grid_base_name, grid_ids):

    existing_tif = []
    grid_id_no_dem_tiff = []
    for id in grid_ids:
        # dem_diff = os.path.join(dem_diff_dir, grid_base_name + '_DEM_diff_grid%d.tif'%id)
        # if os.path.isfile(dem_diff):
        #     existing_tif.append(dem_diff)
        #     continue

        dem_diff_files = io_function.get_file_list_by_pattern(dem_diff_dir, '*_DEM_diff_grid%d.tif'%id)
        if len(dem_diff_files) == 1:
            existing_tif.append(dem_diff_files[0])
            continue
        elif len(dem_diff_files) > 1:
            existing_tif.append(dem_diff_files[0])
            basic.outputlogMessage('warning, There are multiple DEM_diff tif for grid: %d'%id)
            for item in dem_diff_files: basic.outputlogMessage(item)
            continue
        else:
            pass

        grid_id_no_dem_tiff.append(id)
    if len(existing_tif) > 0:
        basic.outputlogMessage('%d existing grid dem diff files for the input grid_ids or extent'%len(existing_tif))
    else:
        basic.outputlogMessage('no existing grid dem diff files')
    return existing_tif, grid_id_no_dem_tiff

def filter_dem_by_month(dem_list):
    allow_months=[6, 7, 8, 9]
    year_dates = [timeTools.get_yeardate_yyyymmdd(os.path.basename(item), pattern='[0-9]{8}_') for item in dem_list]
    month_list = [item.month for item in year_dates]

    old_count = len(dem_list)
    out_dem_list = []
    for dem, month in zip(dem_list,month_list):
        if month in allow_months:
            out_dem_list.append(dem)

    basic.outputlogMessage('Select %d from %d DEM within months %s'%(len(out_dem_list), old_count, str(allow_months)))

    return out_dem_list

def produce_dem_diff_grids(grid_polys, grid_ids, pre_name, reg_tifs,b_apply_matchtag, b_mosaic_id,b_mosaic_date,keep_dem_percent,o_res,process_num=4):

    dem_ext_polys = get_dem_tif_ext_polygons(reg_tifs)
    dem_diff_tifs = []
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
        dem_list_sub = [reg_tifs[index] for index in dem_poly_index]

        # filter by month  # cancel, because it removes many good data.
        # dem_list_sub = filter_dem_by_month(dem_list_sub)

        mosaic_tif_list = mosaic_crop_dem(dem_list_sub, save_dir, grid_id, grid_poly, b_mosaic_id, b_mosaic_date,
                        process_num, keep_dem_percent, o_res, pre_name, resample_method='average',b_mask_matchtag=b_apply_matchtag,
                                          b_mask_stripDEM_outlier=b_mask_stripDEM_outlier,b_mask_surface_water=b_mask_surface_water,
                                          b_mosaic_year=b_mosaic_year)


        # dem co-registration (cancel, the result in not good with the default setting)

        # dem differencing
        save_dem_diff = os.path.join(grid_dem_diffs_dir, pre_name + '_DEM_diff_grid%d.tif'%grid_id)
        save_date_diff = os.path.join(grid_dem_diffs_dir, pre_name + '_date_diff_grid%d.tif'%grid_id)

        if dem_diff_newest_oldest(mosaic_tif_list, save_dem_diff, save_date_diff, process_num,
                               b_max_subsidence=b_max_subsidence,b_save_cm=True):
            dem_diff_tifs.append(save_dem_diff)
    return dem_diff_tifs



def main(options, args):
    extent_shp_or_ids_txt = args[0]
    process_num = options.process_num
    keep_dem_percent = options.keep_dem_percent
    o_res = options.out_res
    basic.setlogfile('produce_DEM_diff_ArcticDEM_log_%s.txt'%timeTools.get_now_time_str())
    if os.path.isdir(grid_dem_diffs_dir) is  False:
        io_function.mkdir(grid_dem_diffs_dir)

    # read grids and ids
    time0 = time.time()
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')
    print('time cost of read polygons and attributes', time.time() - time0)

    # get grid ids based on input extent
    grid_base_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]
    grid_polys, grid_ids = get_grid_20(extent_shp_or_ids_txt,all_grid_polys, all_ids)

    # check dem difference existence
    grid_dem_tifs, grid_ids_no_demDiff = get_existing_dem_diff(grid_dem_diffs_dir,grid_base_name,grid_ids)
    if len(grid_ids_no_demDiff) > 0:
        # refine grid_polys
        if len(grid_ids) > len(grid_ids_no_demDiff):
            id_index = [grid_ids.index(id) for id in grid_ids_no_demDiff]
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
        reg_tifs = [tif for tif in reg_tifs if 'matchtag' not in tif]  # remove matchtag
        # crop, mosacic, difference
        out_dem_diffs = produce_dem_diff_grids(grid_polys, grid_ids_no_demDiff, grid_base_name,reg_tifs,b_apply_matchtag,b_mosaic_id,b_mosaic_date,
                                               keep_dem_percent,o_res,process_num=process_num)
        grid_dem_tifs.extend(out_dem_diffs)





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

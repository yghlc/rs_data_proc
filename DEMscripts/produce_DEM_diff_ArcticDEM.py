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

machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.RSImageProcess as RSImageProcess
import basic_src.map_projection as map_projection


from download_arcticDEM import download_dem_tarball
from ArcticDEM_unpack_registration import process_dem_tarball
from dem_mosaic_crop import mosaic_crop_dem
from dem_mosaic_crop import get_dem_tif_ext_polygons
from dem_mosaic_crop import subset_image_by_polygon_box
from dem_difference import dem_diff_newest_oldest

# some parameters
b_mosaic_id = True
b_mosaic_date = False
b_max_subsidence = False        # apply max_subsidence make results worse

grid_20_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp')
dem_strip_shp = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7/ArcticDEM_Strip_Index_Rel7.shp')

# some folder paths
if machine_name == 'uist':
    arcticDEM_tarball_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/tarballs'
    arcticDEM_reg_tif_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
    grid_dem_diff_dir     = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/grid_dem_diffs'
elif machine_name == 'ubuntu':  # tesia
    arcticDEM_tarball_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/tarballs'
    arcticDEM_reg_tif_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
    grid_dem_diff_dir     = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/grid_dem_diffs'
else:
    raise ValueError('unknown machine:%s'%machine_name)

def get_grid_20(extent_shp_or_id_txt, grid_polys, ids):

    io_function.is_file_exist(extent_shp_or_id_txt)
    file_name_base = os.path.splitext(os.path.basename(extent_shp_or_id_txt))[0]

    if extent_shp_or_id_txt.endswith('.txt'):
        grid_ids = io_function.read_list_from_txt(extent_shp_or_id_txt)
    else:
        # extent polygons and projection (proj4)
        extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp_or_id_txt)
        grid_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)

        if extent_shp_prj != grid_shp_prj:
            basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
                                   % (extent_shp_or_id_txt, grid_20_shp, os.path.basename(extent_shp_or_id_txt)))
            epsg = map_projection.get_raster_or_vector_srs_info_epsg(grid_20_shp)
            # print(epsg)
            # extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_shp_prj.strip())
            extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp_or_id_txt, epsg)
        else:
            extent_polys = vector_gpd.read_polygons_gpd(extent_shp_or_id_txt)

        if len(extent_polys) != 1:
            raise ValueError('Only support one extent polygon')
        grid_index = vector_gpd.get_poly_index_within_extent(grid_polys, extent_polys[0])
        basic.outputlogMessage('find %d grids within the extent (%s)' % (len(grid_index), os.path.basename(extent_shp_or_id_txt)) )

        grid_ids = [ ids[idx] for idx in grid_index]
        grid_ids_str = [str(item) for item in grid_ids ]
        io_function.save_list_to_txt(file_name_base+'_grid_ids.txt',grid_ids_str)

    id_index = [ids.index(id) for id in grid_ids]
    selected_grid_polys = [grid_polys[idx] for idx in id_index ]

    return selected_grid_polys, grid_ids

def get_existing_dem_diff(dem_diff_dir, grid_base_name, grid_ids):

    existing_tif = []
    grid_id_no_dem_tiff = []
    for id in grid_ids:
        dem_diff = os.path.join(dem_diff_dir, grid_base_name + '_DEM_diff_grid%d.tif'%id)
        if os.path.isfile(dem_diff):
            existing_tif.append(dem_diff)
            continue
        grid_id_no_dem_tiff.append(id)
    return existing_tif, grid_id_no_dem_tiff

def produce_dem_diff_grids(grid_polys, grid_ids, pre_name, reg_tifs, b_mosaic_id,b_mosaic_date,keep_dem_percent,o_res,process_num=4):

    dem_ext_polys = get_dem_tif_ext_polygons(reg_tifs)
    dem_diff_tifs = []
    # mosaic and crop
    for grid_id, grid_poly in zip(grid_ids, grid_polys):

        save_dir = 'grid_%d_tmp_files'%grid_id

        # get subset of tifs
        dem_poly_index = vector_gpd.get_poly_index_within_extent(dem_ext_polys, grid_poly)
        if len(dem_poly_index) < 1:
            basic.outputlogMessage('warning, no dem tifs within %d grid, skip' % grid_id)
            continue
        dem_list_sub = [reg_tifs[index] for index in dem_poly_index]

        mosaic_tif_list = mosaic_crop_dem(dem_list_sub, save_dir, grid_id, grid_poly, b_mosaic_id, b_mosaic_date,
                        process_num, keep_dem_percent, o_res, pre_name, resample_method='average')


        # dem co-registration (cancel, the result in not good with the default setting)

        # dem differencing
        save_dem_diff = os.path.join(save_dir, pre_name + '_DEM_diff_grid%d.tif'%grid_id)
        save_date_diff = os.path.join(save_dir, pre_name + '_date_diff_grid%d.tif'%grid_id)

        if dem_diff_newest_oldest(mosaic_tif_list, save_dem_diff, save_date_diff, process_num,
                               b_max_subsidence=b_max_subsidence):
            dem_diff_tifs.append(save_dem_diff)
    return dem_diff_tifs



def main(options, args):
    extent_shp_or_ids_txt = args[0]
    process_num = options.process_num
    keep_dem_percent = options.keep_dem_percent
    o_res = options.out_res

    # read grids and ids
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')

    # get grid ids based on input extent
    grid_base_name = os.path.splitext(os.path.basename(extent_shp_or_ids_txt))[0]
    grid_polys, grid_ids = get_grid_20(extent_shp_or_ids_txt,all_grid_polys, all_ids)

    # check dem difference existence
    grid_dem_tifs, grid_ids_no_demDiff = get_existing_dem_diff(grid_dem_diff_dir,grid_base_name,grid_ids)
    if len(grid_ids_no_demDiff) > 0:
        # refine grid_polys
        if len(grid_ids) > len(grid_ids_no_demDiff):
            id_index = [grid_ids.index(id) for id in grid_ids_no_demDiff]
            grid_polys = [grid_polys[idx] for idx in id_index]

        # download ArcticDEM and applying registration
        tarballs, reg_tifs = download_dem_tarball(dem_strip_shp, grid_polys, arcticDEM_tarball_dir, grid_base_name, reg_tif_dir=arcticDEM_reg_tif_dir)

        # unpack and applying registration
        if len(tarballs) > 0:
            basic.outputlogMessage('Processs %d dem tarballs'%len(tarballs))
            out_reg_tifs = process_dem_tarball(tarballs,'./',arcticDEM_reg_tif_dir,remove_inter_data=True, apply_registration=True)
            basic.outputlogMessage('Get %d new registration dem tifs' % len(out_reg_tifs))
            reg_tifs.extend(out_reg_tifs)

        # crop, mosacic, difference
        out_dem_diffs = produce_dem_diff_grids(grid_polys, grid_ids_no_demDiff, grid_base_name,reg_tifs,b_mosaic_id,b_mosaic_date,
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

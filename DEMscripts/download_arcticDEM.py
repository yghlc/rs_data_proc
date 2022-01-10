#!/usr/bin/env python
# Filename: download_arcticDEM 
"""
introduction: download arcticDEM (strip or mosaic) for a specific region

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 25 December, 2020
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.basic as basic

from urllib.parse import urlparse

import time

from dem_common import tarball_dir,arcticDEM_reg_tif_dir,arcticDEM_tile_tarball_dir,arcticDEM_tile_reg_tif_dir

from dem_common import grid_no_dem_txt, process_log_dir

def get_total_size(url_list):
    total_size = 0
    for url in url_list:
        size = io_function.get_url_file_size(url)       # bytes
        if size is not False:
            total_size += size
    return total_size/(1024.0*1024.0*1024.0)    # GB


def save_id_grid_no_dem(grid_id):
    # grid_dem_diff_less2dem_txt
    if os.path.isdir(process_log_dir) is False:
        io_function.mkdir(process_log_dir)
    # update grid_dem_diff_less2dem_txt file
    id_list = []
    if os.path.isfile(grid_no_dem_txt):
        id_list = io_function.read_list_from_txt(grid_no_dem_txt)    # no need covert to int
    id_str = str(grid_id)
    if id_str in id_list:
        return True
    else:
        # save
        id_list.append(str(grid_id))
        io_function.save_list_to_txt(grid_no_dem_txt,id_list)
        basic.outputlogMessage('Save gird id (%d) to %s' % (grid_id,grid_no_dem_txt))
        return True

def wget_file_url(url):
    cmd_str = 'wget --no-check-certificate %s' % url
    status, result = basic.exec_command_string(cmd_str)
    return status, result

def download_dem_tarball(dem_index_shp, extent_polys, save_folder, pre_name, reg_tif_dir=None, poly_ids=None,b_arcticDEM_tile=False):
    # read dem polygons and url
    dem_polygons, dem_urls = vector_gpd.read_polygons_attributes_list(dem_index_shp, 'fileurl',b_fix_invalid_polygon=False)

    basic.outputlogMessage('%d dem polygons in %s' % (len(dem_polygons), dem_index_shp))

    dem_tar_ball_list = []
    reg_tifs_list = []
    curr_dir = os.getcwd()
    b_save_grid_id_noDEM = True
    if poly_ids is None:
        poly_ids = [idx for idx in range(len(extent_polys)) ]
        b_save_grid_id_noDEM = False    # if poly_ids is not the global unique id, then don't save it.

    for count, (idx, ext_poly) in enumerate(zip(poly_ids, extent_polys)):
        basic.outputlogMessage('get data for the %d th extent (%d/%d)' % (idx, count, len(extent_polys)))

        save_txt_path = pre_name + '_dem_urls_poly_%d.txt' % idx
        if os.path.isfile(save_txt_path):
            urls = io_function.read_list_from_txt(save_txt_path)
            basic.outputlogMessage('read %d dem urls from %s' % (len(urls),save_txt_path))
        else:
            # get fileurl
            dem_poly_ids = vector_gpd.get_poly_index_within_extent(dem_polygons,ext_poly)
            basic.outputlogMessage('find %d DEM within %d th extent' % (len(dem_poly_ids), (idx )))
            urls = [dem_urls[id] for id in dem_poly_ids]

            # save to txt
            io_function.save_list_to_txt(save_txt_path, urls)
            basic.outputlogMessage('save dem urls to %s' % save_txt_path)

        if len(urls) > 0:

            # total_size_GB = get_total_size(urls)  # internet access, parallel running may cause problem. The info is not important
            # basic.outputlogMessage('the size of files will be downloaded is %.4lf GB for the %d th extent '%(total_size_GB,(idx+1)))
            # time.sleep(5)   # wait 5 seconds

            # download them using wget one by one
            for ii, url in enumerate(urls):
                tmp = urlparse(url)
                filename = os.path.basename(tmp.path)
                save_dem_path = os.path.join(save_folder,filename)
                if reg_tif_dir is not None:
                    tar_base = os.path.basename(filename)[:-7]
                    # file_pattern = ['*dem_reg.tif', '*reg_dem.tif'] # Arctic strip and tile (mosaic) version
                    if b_arcticDEM_tile:
                        reg_tifs = io_function.get_file_list_by_pattern(reg_tif_dir, tar_base + '*reg_dem.tif')
                    else:
                        reg_tifs = io_function.get_file_list_by_pattern(reg_tif_dir,tar_base+'*dem_reg.tif')
                    if len(reg_tifs) > 0:
                        basic.outputlogMessage('warning, unpack and registrated tif for %s already exists, skip downloading' % filename)
                        reg_tifs_list.append(reg_tifs[0])
                        continue

                if os.path.isfile(save_dem_path):
                    basic.outputlogMessage('warning, %s already exists, skip downloading'%filename)
                else:
                    # download the dem
                    basic.outputlogMessage('starting downloading %d th DEM (%d in total)'%((ii+1),len(urls)))
                    os.chdir(save_folder)
                    status, result = wget_file_url(url)
                    # try new url if it exists
                    if status != 0 and "302 Moved Temporarily" in result:
                        str_list = result.split()
                        loc_idx = str_list.index('Location:')
                        new_url = str_list[loc_idx + 1]  # find the new URL
                        print('try to download the file using new url: %s'%new_url)
                        status, result = wget_file_url(new_url)

                    if status != 0:
                        print(result)
                        # sys.exit(status)
                        basic.outputlogMessage('failed to download DEM: %s'%url)
                    os.chdir(curr_dir)

                dem_tar_ball_list.append(save_dem_path)

        else:
            basic.outputlogMessage('Warning, can not find DEMs within %d th extent'%(idx))
            if b_save_grid_id_noDEM:
                save_id_grid_no_dem(idx)

    return dem_tar_ball_list, reg_tifs_list

def main(options, args):

    extent_shp = args[0]
    dem_index_shp = args[1]
    save_folder = options.save_dir
    b_arcticDEM_tile = False
    if save_folder is None:
        if 'Tile' in os.path.basename(dem_index_shp):
            save_folder =arcticDEM_tile_tarball_dir
            reg_tif_dir = arcticDEM_tile_reg_tif_dir
            b_arcticDEM_tile = True
        else:
            save_folder = tarball_dir
            reg_tif_dir = arcticDEM_reg_tif_dir

    pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    pre_name += '_Tile' if 'Tile' in os.path.basename(dem_index_shp) else '_Strip'

    # extent polygons and projection (proj4)
    extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp)
    dem_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(dem_index_shp)

    if extent_shp_prj != dem_shp_prj:
        basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
                               %(extent_shp,dem_index_shp,os.path.basename(extent_shp)))
        epsg = map_projection.get_raster_or_vector_srs_info_epsg(dem_index_shp)
        # print(epsg)
        # extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_shp_prj.strip())
        extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,epsg)
    else:
        extent_polys = vector_gpd.read_polygons_gpd(extent_shp)

    # read 'grid_id' if the extent shp is from grid shp file, if not, grid_id_list will be None
    grid_id_list = vector_gpd.read_attribute_values_list(extent_shp,'grid_id')

    if len(extent_polys) < 1:
        raise ValueError('No polygons in %s'%extent_shp)
    else:
        basic.outputlogMessage('%d extent polygons in %s'%(len(extent_polys),extent_shp))

    download_dem_tarball(dem_index_shp,extent_polys,save_folder,pre_name,reg_tif_dir=reg_tif_dir,
                         poly_ids=grid_id_list,b_arcticDEM_tile=b_arcticDEM_tile)


if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp dem_indexes_shp"
    parser = OptionParser(usage=usage, version="1.0 2020-12-25")
    parser.description = 'Introduction: download ArcticDEM within an extent  '
    # parser.add_option("-x", "--save_xlsx_path",
    #                   action="store", dest="save_xlsx_path",
    #                   help="save the sence lists to xlsx file")

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",
                      help="the folder to save DEMs")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

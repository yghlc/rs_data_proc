#!/usr/bin/env python
# Filename: get_ArcticDEM_filesize.py
"""
introduction: get ArcticDEM filesize  for a specific region

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

from datetime import datetime
import time
import math

from dem_common import arcticDEM_reg_tif_dir,arcticDEM_tile_reg_tif_dir

from ArcticDEM_unpack_registration import is_ArcticDEM_tiles

machine_name = os.uname()[1]
import pandas as pd

from download_arcticDEM import get_total_size


def get_one_url_file_size(url, ii, total):
    try:
        print('%d/%d'%(ii,total))
        size = io_function.get_url_file_size(url)
        if size is not False:
            return size/(1024.0*1024.0*1024.0)    # GB
    except:
        print(datetime.now(),'get file size of %s failed'%url)
        return None


def get_file_size_dem_tarball(dem_index_shp, extent_polys, pre_name, xlsx_size_path,poly_ids=None):
    # read dem polygons and url
    dem_polygons, dem_urls = vector_gpd.read_polygons_attributes_list(dem_index_shp, 'fileurl',b_fix_invalid_polygon=False)

    if os.path.isfile(xlsx_size_path):
        size_table = pd.read_excel(xlsx_size_path)
        save_idx_list = size_table['index'].to_list()
        save_url_list = size_table['fileurl'].to_list()
        save_size_list = size_table['filesize'].to_list()
    else:
        save_idx_list = [item for item in range(len(dem_urls))]  # index list
        save_url_list = dem_urls
        save_size_list = [float('nan')]*len(save_idx_list)

    basic.outputlogMessage('%d dem polygons in %s' % (len(dem_polygons), dem_index_shp))

    if poly_ids is None:
        poly_ids = [ item for item in range(len(extent_polys)) ]

    for count, (idx, ext_poly) in enumerate(zip(poly_ids, extent_polys)):
        basic.outputlogMessage('get ArcticDEM filesize for the %d th extent (%d/%d)' % (idx, count, len(extent_polys)))

        save_filesize_txt = pre_name + '_dem_FileSize_poly_%d.txt' % idx
        if os.path.isfile(save_filesize_txt):
            basic.outputlogMessage('%s exists, skip'%save_filesize_txt)
            continue

        # get fileurl
        dem_poly_idx_list = vector_gpd.get_poly_index_within_extent(dem_polygons,ext_poly)
        basic.outputlogMessage('find %d DEM within %d th extent' % (len(dem_poly_idx_list), (idx )))
        urls = [dem_urls[id] for id in dem_poly_idx_list]
        url_size_list = [ save_size_list[id] for id in dem_poly_idx_list ]

        if len(urls) > 0:
            total_count = len(urls)
            for ii, (url, fileS,url_idx) in enumerate(zip(urls,url_size_list,dem_poly_idx_list)):
                # remove url start with /mnt and end with .tif
                if url.startswith('/mnt') and url.endswith('.tif'):
                    basic.outputlogMessage("error: not a valid url: %s"%url)
                    continue
                if math.isnan(fileS) is False:
                    continue
                url_size_GB = get_one_url_file_size(url,ii,total_count)
                url_size_list[ii] = url_size_GB
                save_size_list[url_idx] = url_size_GB

            url_size_list_noNone = [item for item in url_size_list if math.isnan(item) is False ]

            if len(url_size_list_noNone) != len(url_size_list):
                basic.outputlogMessage('There are %d None value in url_size_list'%(len(url_size_list) - len(url_size_list_noNone)))

            total_size_GB = sum(url_size_list_noNone)

            basic.outputlogMessage('the size of files will be downloaded is %.4lf GB for the %d th extent '%(total_size_GB,(idx+1)))
            with open(save_filesize_txt,'w') as f_obj:
                f_obj.writelines('%d DEM files, total size is  %.6lf GB \n'%(len(urls),total_size_GB))
        else:
            basic.outputlogMessage('Warning, can not find DEMs within %d th extent'%(idx))


    # save table
    save_dict = {'index':save_idx_list, 'filesize':save_size_list, 'fileurl':save_url_list}

    save_dict_pd = pd.DataFrame(save_dict)
    with pd.ExcelWriter(xlsx_size_path) as writer:
        save_dict_pd.to_excel(writer, sheet_name='url_file_size')

    return None

def main(options, args):

    extent_shp = args[0]
    dem_index_shp = args[1]

    pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    pre_name += '_Tile' if 'Tile' in os.path.basename(dem_index_shp) else '_Strip'

    xlsx_size_path = os.path.splitext(os.path.basename(dem_index_shp))[0] + '_fileSize.xlsx'
    print('xlsx file for saving file size',xlsx_size_path)

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

    get_file_size_dem_tarball(dem_index_shp,extent_polys,pre_name,xlsx_size_path,poly_ids=grid_id_list)


if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp dem_indexes_shp"
    parser = OptionParser(usage=usage, version="1.0 2022-2-19")
    parser.description = 'Introduction: get the file size of ArcticDEM within an extent  '
    parser.add_option("-x", "--save_xlsx_path",
                      action="store", dest="save_xlsx_path",
                      help="save the DEM lists to xlsx file")

    # parser.add_option("-m", "--max_process_num",
    #                   action="store", dest="max_process_num", type=int,default=16,
    #                   help="the maximum number of processes for downloading in parallel")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

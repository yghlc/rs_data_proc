#!/usr/bin/env python
# Filename: get_ArcticDEM_filesize.py
"""
introduction: statistics the file informatoin of ArcticDEM files from dem_index_shp or gpkg files

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 25 January, 2026
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.basic as basic


machine_name = os.uname()[1]
import pandas as pd
import geopandas as gpd

from tqdm import tqdm


def statistics_on_dems_in_gpkg(gpkg_dir, save_xlsx_path=None):
    # statistics the dem files information from gpkg files in a directory
    gpkg_list = io_function.get_file_list_by_pattern(gpkg_dir,'*.gpkg')
    basic.outputlogMessage(f'Find {len(gpkg_list)} gpkg files in {gpkg_dir}')

    records = []

    for gpkg in tqdm(gpkg_list, desc="Processing gpkg files"):
        gpkg_gpd = gpd.read_file(gpkg)


        if 'pgc:is_xtrack' not in gpkg_gpd.columns:
            basic.outputlogMessage(f'No is_xtrack field in {gpkg}, skip it')
            continue
        # get the count of is_xtrack is true and false
        count_xtrack_true = gpkg_gpd['pgc:is_xtrack'].sum()
        count_xtrack_false = len(gpkg_gpd) - count_xtrack_true
        true_percent = count_xtrack_true / len(gpkg_gpd) * 100.0
        false_percent = count_xtrack_false / len(gpkg_gpd) * 100.0
        record = {'gpkg_file': os.path.basename(gpkg),
                  'total_dems': len(gpkg_gpd),  
                    'count_xtrack_true': count_xtrack_true,
                    'count_xtrack_false': count_xtrack_false,
                    'true_percent': true_percent,
                    'false_percent': false_percent
                    }
        records.append(record)

    if len(records) < 1:
        basic.outputlogMessage('No records found from gpkg files')
        return

    result_df = pd.DataFrame.from_records(records)
    basic.outputlogMessage(f'Total {len(result_df)} dem records from {len(gpkg_list)} gpkg files')

    if save_xlsx_path is not None:
        result_df.to_excel(save_xlsx_path, index=False)
        basic.outputlogMessage(f'Saved the statistics to {save_xlsx_path}')

def main(options, args):

    if len(args) > 1:
        extent_shp = args[0]
        dem_index_shp = args[1]
        # need to use the functions in get_ArcticDEM_filesize.py
        return
    
        # pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
        # pre_name += '_Tile' if 'Tile' in os.path.basename(dem_index_shp) else '_Strip'

        # xlsx_size_path = os.path.splitext(os.path.basename(dem_index_shp))[0] + '_fileSize.xlsx'
        # print('xlsx file for saving file size',xlsx_size_path)

        # # extent polygons and projection (proj4)
        # extent_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(extent_shp)
        # dem_shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(dem_index_shp)

        # if extent_shp_prj != dem_shp_prj:
        #     basic.outputlogMessage('%s and %s do not have the same projection, will reproject %s'
        #                         %(extent_shp,dem_index_shp,os.path.basename(extent_shp)))
        #     epsg = map_projection.get_raster_or_vector_srs_info_epsg(dem_index_shp)
        #     # print(epsg)
        #     # extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_shp_prj.strip())
        #     extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,epsg)
        # else:
        #     extent_polys = vector_gpd.read_polygons_gpd(extent_shp)

        # # read 'grid_id' if the extent shp is from grid shp file, if not, grid_id_list will be None
        # grid_id_list = vector_gpd.read_attribute_values_list(extent_shp,'grid_id')

        # if len(extent_polys) < 1:
        #     raise ValueError('No polygons in %s'%extent_shp)
        # else:
        #     basic.outputlogMessage('%d extent polygons in %s'%(len(extent_polys),extent_shp))

    else:
        # just statistics from gpkg file (saved when downloading DEM files using STAC)
        gpkg_dir = options.gpkg_dir
        save_xlsx_path = options.save_xlsx_path if options.save_xlsx_path is not None else os.path.basename(gpkg_dir)+'_dem_info_statistic.xlsx'

        statistics_on_dems_in_gpkg(gpkg_dir, save_xlsx_path=save_xlsx_path)



if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp dem_index_shp "
    parser = OptionParser(usage=usage, version="1.0 2026-1-25")
    parser.description = 'Introduction: get and statistics ArcticDEM information within an extent '

    parser.add_option("-x", "--save_xlsx_path",
                      action="store", dest="save_xlsx_path",
                      help="save the DEM lists to xlsx file")
    
    parser.add_option("-d", "--gpkg_dir",
                      action="store", dest="gpkg_dir",
                      help="the directory of gpkg files, if set, will ignore dem_index_shps")

    # parser.add_option("-m", "--max_process_num",
    #                   action="store", dest="max_process_num", type=int,default=16,
    #                   help="the maximum number of processes for downloading in parallel")


    (options, args) = parser.parse_args()
    # if len(sys.argv) < 2 or len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)

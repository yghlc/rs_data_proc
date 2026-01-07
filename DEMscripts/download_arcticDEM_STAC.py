#!/usr/bin/env python
# Filename: download_arcticDEM_STAC.py 
"""
introduction: download ArcticDEM using STAC (The SpatioTemporal Asset Catalog) for a specific region

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 06 January, 2026
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
# import basic_src.map_projection as map_projection
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools
# import basic_src.basic as basic

from datetime import datetime

import time
pgc_stac = "https://stac.pgc.umn.edu/api/v1/"

def get_collection_for_search(collection_id='arcticdem-strips-s2s041-2m'):

    cat = pystac_client.Client.open(pgc_stac)
    # print("Catalog Title: {}".format(cat.title))
    # print(cat)
    #
    # for collection in cat.get_collections():
    #     print(collection)
    #
    # collection_search = cat.collection_search(q="arcticdem",)
    # for result in collection_search.collections():
    #     print(result.id, f"{result.description}", sep="\n")
    #     print("\n")

    ######################################################################
    # collection id for ArcticDEM:
    # arcticdem-mosaics-v4.1-2m (ArcticDEM DEM mosaics, version 4.1, 2m resolution)
    # arcticdem-mosaics-v3.0-2m (ArcticDEM DEM mosaics, version 3.0, 2m resolution)
    # arcticdem-mosaics-v3.0-32m (ArcticDEM DEM mosaics, version 3.0, 32m resolution)
    # arcticdem-mosaics-v4.1-10m (ArcticDEM DEM mosaics, version 4.1, 10m resolution)
    # arcticdem-strips-s2s041-2m (ArcticDEM time-stamped strip DEMs, s2s version 4.1, 2m resolution)
    # arcticdem-mosaics-v3.0-10m (ArcticDEM DEM mosaics, version 3.0, 10m resolution)
    # arcticdem-mosaics-v4.1-32m (ArcticDEM DEM mosaics, version 4.1, 32m resolution)
    ######################################################################

    # STAC item
    data_collection = cat.get_collection(collection_id)
    # print(data_collection)
    return data_collection

def download_dem_within_polygon(client,collection_id, poly_extent, date_start='2008-01-01', date_end='2026-12-31',
                                search_save='tmp.gpkg', save_crs_code=3413,save_dir='data_save'):

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    # else:
    #     if len(os.listdir(save_dir)) > 0:
    #         print(f'the saved folder: {save_dir} exists and not empty, skip downloading')
    #         return

    bbox = vector_gpd.get_polygon_bounding_box(poly_extent)
    search = client.search(
        collections=[collection_id],
        bbox=bbox,
        datetime=f"{date_start}/{date_end}"
    )

    # fetch the items that fit our search parameters
    items = list(search.items())
    item_count = len(items)
    print(f'Found {item_count} items')
    if item_count < 1:
        return

    items_gdf = gpd.GeoDataFrame.from_features(search.item_collection().to_dict(), crs="epsg:4326").to_crs(save_crs_code)

    items_gdf.to_file(search_save)
    print(f'saved search results to {search_save}')

    stack = stackstac.stack(items, epsg=save_crs_code, bounds_latlon=bbox)
    # print(stack)
    # attr_list = ["band","time", 'data_type', 'nodata','unit', 'epsg' ]
    # for att in attr_list:
    #     print('Attribute:',att,':')
    #     print(stack[att].values)
    bands_to_save = ['dem','mask']
    bands_to_save_idx = [stack['band'].values.tolist().index(item) for item in bands_to_save]
    data_types = [ stack['data_type'].values.tolist()[idx] for idx in bands_to_save_idx]
    # print(data_types)


    for img_time in stack['time'].values:
        # print(img_time,type(img_time))
        # img_time is "np.datetime64",  astype(str) to str, astype(datetime) to datetime.
        dt_obj = timeTools.isoformat_str_2_datetime(img_time.astype(str))
        # dt_obj = img_time.astype(datetime)
        dt_str = timeTools.datetime2str(dt_obj,format='%Y%m%d_%H%M%S')
        # print(dt_str)
        for band, d_type in zip(bands_to_save, data_types):
            img_save_path = os.path.join(save_dir,f'{collection_id}_{dt_str}_{band}.tif')
            if os.path.isfile(img_save_path):
                print(f'warning, {img_save_path} exist, skip downloading')
                continue

            selected = stack.sel(band=band, time=img_time)
            # Ensure the DataArray has spatial metadata
            selected.rio.write_crs(stack.attrs['crs'], inplace=True)  # or use arr.rio.crs if available

            # Save to GeoTIFF
            selected = selected.astype(d_type)
            selected.rio.to_raster(img_save_path, compress="LZW")
            print(f'saved geotiff to {img_save_path}')

        # break   # for testing



def main(options, args):
    extent_shp = args[0]
    print(extent_shp)
    # STAC need lat/lon for searching?
    ext_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,'EPSG:4326')
    # print(ext_polys)
    base_name = io_function.get_name_no_ext(extent_shp)
    data_dir = './'

    collection_id = 'arcticdem-strips-s2s041-2m'
    client = pystac_client.Client.open(pgc_stac)

    # b_arcticDEM_tile = False
    # collection_client = get_collection_for_search(collection_id='arcticdem-strips-s2s041-2m')
    # print(collection_client)

    for idx, poly in enumerate(ext_polys):
        search_save = base_name + f'_poly_{idx+1}.gpkg'
        save_dir = os.path.join(data_dir, base_name + f'_poly_{idx+1}')
        download_dem_within_polygon(client,collection_id,  poly, date_start='2008-01-01',
                                    date_end='2026-12-31', search_save=search_save, save_dir=save_dir)



if __name__ == '__main__':
    # import need library
    import pystac_client
    import geopandas as gpd
    import stackstac
    import xarray as xr
    import rioxarray


    usage = "usage: %prog [options] extent_shp "
    parser = OptionParser(usage=usage, version="1.0 2026-01-06")
    parser.description = 'Introduction: download ArcticDEM within an extent  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",
                      help="the folder to save DEMs")

    parser.add_option("-m", "--max_process_num",
                      action="store", dest="max_process_num", type=int, default=8,
                      help="the maximum number of processes for downloading in parallel")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

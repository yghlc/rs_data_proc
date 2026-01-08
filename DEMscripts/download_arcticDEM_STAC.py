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
import basic_src.basic as basic

from datetime import datetime

import time
pgc_stac = "https://stac.pgc.umn.edu/api/v1/"


from dem_common import arcticDEM_reg_tif_dir,arcticDEM_tile_reg_tif_dir

from download_arcticDEM import save_id_grid_no_dem

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

def get_collection_examples_meta(client,collection_id, poly_extent, ext_id, date_start='2008-01-01', date_end='2026-12-31',
                                search_save='tmp.gpkg', save_crs_code=3413):

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
        basic.outputlogMessage(
            f'Warning, can not find DEMs within {ext_id} th extent with: {collection_id} from {date_start} to {date_end}  ')
        return False

    items_gdf = gpd.GeoDataFrame.from_features(search.item_collection().to_dict(), crs="epsg:4326").to_crs(
        save_crs_code)

    items_gdf.to_file(search_save)
    basic.outputlogMessage(f'saved search results to {search_save}')

    stack = stackstac.stack(items, epsg=save_crs_code, bounds_latlon=bbox)
    print(stack)
    attr_list = ["band","time", 'data_type', 'nodata','unit', 'epsg', 'id']
    for att in attr_list:
        print('Attribute:',att,':', type(stack[att].values))
        print(stack[att].values)
    sys.exit(0)

def get_bands_to_save(collection_id):
    # we can check these by printing some examples and QGIS STAC plugin
    if collection_id == 'arcticdem-strips-s2s041-2m':
        return ['dem','mask','matchtag']    # ['dem' 'mask' 'matchtag' 'hillshade' 'hillshade_masked']
    elif collection_id == 'arcticdem-mosaics-v4.1-2m':
        return ['dem']  # ['dem' 'mad' 'count' 'maxdate' 'mindate' 'datamask' 'hillshade']
    else:
        raise ValueError(f'Unknown correction id: {collection_id}')

def save_one_image_to_local(stack,selected,d_type,img_save_path):
    # Ensure the DataArray has spatial metadata
    selected.rio.write_crs(stack.attrs['crs'], inplace=True)  # or use arr.rio.crs if available

    # Save to GeoTIFF
    selected = selected.astype(d_type)
    selected.rio.to_raster(img_save_path, compress="LZW")
    basic.outputlogMessage(f'saved geotiff to {img_save_path}')

def download_dem_within_polygon(client,collection_id, poly_extent, ext_id, date_start='2008-01-01', date_end='2026-12-31',
                                search_save='tmp.gpkg', save_crs_code=3413,save_dir='data_save',b_unique_grid=False):

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

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
        basic.outputlogMessage(f'Warning, can not find DEMs within {ext_id} th extent with: {collection_id} from {date_start} to {date_end}  ')
        return False

    items_gdf = gpd.GeoDataFrame.from_features(search.item_collection().to_dict(), crs="epsg:4326").to_crs(save_crs_code)

    items_gdf.to_file(search_save)
    basic.outputlogMessage(f'saved search results to {search_save}')

    stack = stackstac.stack(items, epsg=save_crs_code, bounds_latlon=bbox)
    # print(stack)
    # attr_list = ["band","time", 'data_type', 'nodata','unit', 'epsg', 'id']
    # for att in attr_list:
    #     print('Attribute:',att,':', type(stack[att].values))
    #     print(stack[att].values)
    # sys.exit(0)
    bands_to_save = get_bands_to_save(collection_id)
    bands_to_save_idx = [stack['band'].values.tolist().index(item) for item in bands_to_save]
    data_types = [ stack['data_type'].values.tolist()[idx] for idx in bands_to_save_idx]

    # file name use image ID
    for img_time, img_id in zip(stack['time'].values,stack['id'].values):
        # print(img_id, type(img_id)) # numpy.str_
        for band, d_type in zip(bands_to_save, data_types):
            if b_unique_grid:
                img_save_path = os.path.join(save_dir, f'{img_id}_grid{ext_id}_{band}.tif')
            else:
                img_save_path = os.path.join(save_dir,f'{img_id}_sub{ext_id}_{band}.tif')
            if os.path.isfile(img_save_path):
                basic.outputlogMessage(f'warning, {img_save_path} exist, skip downloading')
                continue

            selected = stack.sel(band=band, time=img_time)
            save_one_image_to_local(stack, selected, d_type, img_save_path)

        # break # for testing

    return True


def download_dem_stac(client,collection_id,extent_polys, output_dir, date_start,date_end,pre_name,
                      poly_ids=None,save_prj_code=3413,b_unique_grid=False):

    # download data through STAC, not need to unpack
    b_save_grid_id_noDEM = True
    if poly_ids is None:
        poly_ids = [idx for idx in range(len(extent_polys)) ]
        b_save_grid_id_noDEM = False    # if poly_ids is not the global unique id, then don't save it.

    for count, (idx, ext_poly) in enumerate(zip(poly_ids, extent_polys)):
        if b_unique_grid:
            search_save = os.path.join(output_dir, pre_name + f'_grid{idx}.gpkg')
        else:
            search_save = os.path.join(output_dir, pre_name + f'_poly_{idx}.gpkg')
        # save_dir = os.path.join(output_dir, pre_name + f'_poly_{idx+1}')

        # # for testing, output information, will quit after running this
        # get_collection_examples_meta(client,collection_id, ext_poly, idx, date_start=date_start,
        #                             date_end=date_end, search_save=search_save,save_crs_code=save_prj_code)

        res = download_dem_within_polygon(client,collection_id, ext_poly, idx, date_start=date_start,
                                    date_end=date_end, search_save=search_save, save_dir=output_dir,save_crs_code=save_prj_code,
                                          b_unique_grid=b_unique_grid)

        # if res is False, mean no data found
        if res is False and b_save_grid_id_noDEM is True:
            save_id_grid_no_dem(idx)


def main(options, args):
    extent_shp = args[0]
    print('input extent:', extent_shp)
    b_arcticDEM_mosaic = options.b_arcticDEM_mosaic # default is False

    global max_task_count
    max_task_count = options.max_process_num

    # reg_tif_dir is the folder to save tifs
    if b_arcticDEM_mosaic:
        reg_tif_dir = arcticDEM_tile_reg_tif_dir
    else:
        reg_tif_dir = arcticDEM_reg_tif_dir

    # use the user specific save_dir if it's set
    if options.save_dir is not None:
        reg_tif_dir = options.save_dir
    if os.path.isdir(reg_tif_dir) is False:
        io_function.mkdir(reg_tif_dir)
    reg_tif_dir = os.path.abspath(reg_tif_dir)  # change to absolute path

    # STAC need lat/lon for searching?
    ext_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,'EPSG:4326')


    if b_arcticDEM_mosaic:
        collection_id = 'arcticdem-mosaics-v4.1-2m'  # ArcticDEM mosaic at 2 m
    else:
        collection_id = 'arcticdem-strips-s2s041-2m'  # ArcticDEM strip at 2 m
    if options.collection_id is not None:
        collection_id = options.collection_id

    pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    pre_name += '_Mosaic' if b_arcticDEM_mosaic else '_Strip'

    # read 'grid_id' if the extent shp is from grid shp file, if not, grid_id_list will be None
    grid_id_list = vector_gpd.read_attribute_values_list(extent_shp,'grid_id')
    b_unique_grid = False if grid_id_list is None else True
    if len(ext_polys) < 1:
        raise ValueError('No polygons in %s'%extent_shp)
    else:
        basic.outputlogMessage('read %d extent polygons in %s for downloading using STAC'%(len(ext_polys),extent_shp))

    date_start = options.date_start
    date_end = options.date_end

    # projection for ArcticDEM products
    save_prj_code = 3413

    client = pystac_client.Client.open(pgc_stac)

    # b_arcticDEM_tile = False
    # collection_client = get_collection_for_search(collection_id='arcticdem-strips-s2s041-2m')
    # print(collection_client)

    download_dem_stac(client, collection_id, ext_polys, reg_tif_dir, date_start, date_end, pre_name,
                      poly_ids=grid_id_list, save_prj_code=save_prj_code,b_unique_grid=b_unique_grid)




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

    parser.add_option("-c", "--collection_id",
                      action="store", dest="collection_id",
                      help="the collection id, set this allow to download other products from PGC")

    parser.add_option("-s", "--date_start",
                      action="store", dest="date_start", type=str, default='2008-01-01',
                      help="the start date for downloading data")

    parser.add_option("-e", "--date_end",
                      action="store", dest="date_end", type=str, default='2026-12-31',
                      help="the end date for downloading data")


    parser.add_option("-p", "--max_process_num",
                      action="store", dest="max_process_num", type=int, default=8,
                      help="the maximum number of processes for downloading in parallel")

    parser.add_option("-m", "--b_arcticDEM_mosaic",
                      action="store_true", dest="b_arcticDEM_mosaic",default=False,
                      help="if set, will download the mosaic of ArcticDEM, other wise, download the strips")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

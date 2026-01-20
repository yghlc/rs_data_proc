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
import multiprocessing
from multiprocessing import Process

import time
pgc_stac = "https://stac.pgc.umn.edu/api/v1/"


from dem_common import arcticDEM_reg_tif_dir,arcticDEM_tile_reg_tif_dir

from download_arcticDEM import save_id_grid_no_dem

machine_name = os.uname()[1]
# the maximum number of processes for downloading in parallel
max_task_count = 1
download_tasks = []

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

def save_one_image_to_local(stack,selected,d_type,img_save_path,nodata_value=None,crop_poly=None,
                            out_res=None,save_txt=None):
    # for the case running this function in sub-process
    if multiprocessing.current_process().name != 'MainProcess':
        import rioxarray

        # Ensure the DataArray has spatial metadata
    selected.rio.write_crs(stack.attrs['crs'], inplace=True)  # or use arr.rio.crs if available

    if out_res is not None:
        # in the output we to resolution of 2.01, 2.023, ect, be consistent, reproject to 2.0
        selected = selected.rio.reproject(selected.rio.crs,resolution=out_res)

    # Crop to polygon if provided
    if crop_poly is not None:
        selected = selected.rio.clip([crop_poly], stack.attrs['crs'])

    # Set nodata value if provided
    if nodata_value is not None:
        selected = selected.rio.write_nodata(nodata_value)

    # Save to GeoTIFF
    selected = selected.astype(d_type)
    selected.rio.to_raster(img_save_path, compress="LZW")
    basic.outputlogMessage(f'saved geotiff to {img_save_path}')
    # for the case running this function in sub-process
    # print('debugging:',multiprocessing.current_process().name)
    # if multiprocessing.current_process().name != 'MainProcess':
    #     with open(save_file_list_txt_a_poly,'a') as f_obj:
    #         f_obj.writelines(img_save_path + '\n')
    # else:
    #     return img_save_path
    if save_txt is not None:
        with open(save_txt, 'a') as f_obj:
            f_obj.writelines(img_save_path + '\n')

    return img_save_path

def filter_search_results_by_polygon(items, search_result_dict, poly_latlon, poly_prj, prj_crs_code, min_overlap_per=0.1):

    # search_result_dict contain two elements: 'type', 'features' 

    items_gdf = gpd.GeoDataFrame.from_features(search_result_dict, crs="epsg:4326").to_crs(prj_crs_code)
    print('Total items before filtering:', len(items_gdf), len(items), len(search_result_dict['features']))
    select_items = []
    select_search_result_dict = {'type': search_result_dict['type'], 'features': []}
    for item, search_dict, item_poly in zip(items,search_result_dict['features'],items_gdf.geometry):
        # print(item)
        # print(search_dict)
        # print(item_poly)
        inter_geo = item_poly.intersection(poly_prj)
        if inter_geo.is_empty:
            # print('No intersection, skip this item')
            continue
        inter_area = inter_geo.area
        poly_area = poly_prj.area
        overlap_per = inter_area / poly_area
        # print(f'overlap area: {inter_area}, polygon area: {poly_area}, overlap percentage: {overlap_per}')
        if overlap_per >= min_overlap_per:
            select_items.append(item)
            select_search_result_dict['features'].append(search_dict)
            # print('Selected this item')
        else:
            # print('Not enough overlap, skip this item')
            pass
    print(f'Total items after filtering minimum {100*min_overlap_per}% overlap:', len(select_items), len(select_search_result_dict['features']))

    return select_items, select_search_result_dict



def download_dem_within_polygon(client,collection_id, poly_latlon, poly_prj, ext_id, date_start='2008-01-01', date_end='2026-12-31',
                                search_save='tmp.gpkg', save_crs_code=3413,save_dir='data_save',b_unique_grid=False,out_res=None):

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    bbox = vector_gpd.get_polygon_bounding_box(poly_latlon)
    search = client.search(
        collections=[collection_id],
        bbox=bbox,
        datetime=f"{date_start}/{date_end}"
    )

    # fetch the items that fit our search parameters
    items = list(search.items())
    item_count = len(items)
    print(f'Found {item_count} items')

    search_result_dict = search.item_collection().to_dict()
    # io_function.save_dict_to_txt_json(f'search_result_poly_{ext_id}.json', search_result_dict) # for debugging

    # filter the results by checking geometry if they are overlap more than 10%  of the poly_extent
    items, search_result_dict = filter_search_results_by_polygon(items, search_result_dict, poly_latlon, poly_prj, save_crs_code, min_overlap_per=0.05)

    # sys.exit(1)

    if item_count < 1:
        basic.outputlogMessage(f'Warning, can not find DEMs within {ext_id} th extent with: {collection_id} from {date_start} to {date_end}  ')
        return False

    items_gdf = gpd.GeoDataFrame.from_features(search_result_dict, crs="epsg:4326").to_crs(save_crs_code)

    items_gdf.to_file(search_save)
    basic.outputlogMessage(f'saved search results to {search_save}')
    # sys.exit(1)

    if item_count > 300:
        basic.outputlogMessage(f'error, too many ({item_count}) DEMs within {ext_id} th extent with: {collection_id} from {date_start} to {date_end}')
        raise ValueError(f'error, too many ({item_count}) DEMs within {ext_id} th extent with: {collection_id} from {date_start} to {date_end}')

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
    nodata_list = [ stack['nodata'].values.tolist()[idx] for idx in bands_to_save_idx]

    save_file_list_txt = os.path.join(os.path.dirname(search_save), io_function.get_name_no_ext(search_save) + '.txt')
    if os.path.isfile(save_file_list_txt):
        io_function.delete_file_or_dir(save_file_list_txt)

    # file name use image ID
    saved_file_list = []
    for img_time, img_id in zip(stack['time'].values,stack['id'].values):
        # print(img_id, type(img_id)) # numpy.str_
        for band, d_type,nodata in zip(bands_to_save, data_types,nodata_list):
            if b_unique_grid:
                img_save_path = os.path.join(save_dir, f'{img_id}_grid{ext_id}_{band}.tif')
            else:
                img_save_path = os.path.join(save_dir,f'{img_id}_sub{ext_id}_{band}.tif')
            if os.path.isfile(img_save_path):
                basic.outputlogMessage(f'warning, {img_save_path} exist, skip downloading')
                saved_file_list.append(img_save_path)
                continue

            # check if enough disk space available
            io_function.wait_until_enough_disk_space(save_dir, min_disk_GB=5)

            # in rare case, img_time is duplicate, end in multime bands with same time stamp
            # caused problem in the later DEM proessing. 
            selected = stack.sel(band=band, time=img_time)
            # to select by band and id instead of time, however, id is not in coordinate
            # print('debuging 1:',selected.dims, selected.shape, ', using: the band and time:', band, img_time)

            # if more than one image, then filter by the ID. 
            if selected.ndim != 2:
                # drop=True to drop the unused dimension, non-matching entries and reduce dimensionality
                selected = selected.where(selected['id'] == img_id, drop=True)
                # print('debuging 2:',selected.dims, selected.shape, ', using: the band and time:', band, img_time)
                # remove the single-dimensional entries from the shape of an array.
                selected = selected.squeeze()
                # print('debuging 3:',selected.dims, selected.shape, ', using: the band and time:', band, img_time)

            if selected.size < 1:
                raise ValueError(f'Selection failed: empty data for band {band}, time {img_time}, id {img_id}')
            if selected.ndim != 2:
                raise ValueError(f'Selection failed: not 2D data for band {band}, time {img_time}, id {img_id}, Dims: {selected.dims}, Shape: {selected.shape}') 


            if max_task_count == 1:
            ########### download one by one ######################
                res_path = save_one_image_to_local(stack, selected, d_type, img_save_path, nodata_value=nodata, crop_poly=poly_prj,
                                                   out_res=out_res,save_txt=save_file_list_txt)
                saved_file_list.append(res_path)
            else:
            ################ parallel download ###############
                basic.check_exitcode_of_process(download_tasks)  # if there is one former job failed, then quit
                while True:
                    job_count = basic.alive_process_count(download_tasks)
                    if job_count >= max_task_count:
                        print(datetime.now(),f'You are running {max_task_count} or more tasks in parallel on {machine_name}, wait ')
                        time.sleep(30)  #
                        continue
                    break
                while True:
                    avail_per = basic.get_available_cpu_mem_per()
                    if avail_per < 30:
                        print(datetime.now(),f'available CPU memory is less than 30%: {avail_per}, wail 10 seconds')
                        time.sleep(10)
                    else:
                        break

                # start the processing

                sub_process = Process(target=save_one_image_to_local, args=(stack, selected, d_type, img_save_path,nodata,
                                                                            poly_prj,out_res,save_file_list_txt))
                sub_process.start() # start a process, don't wait
                download_tasks.append(sub_process)

                basic.close_remove_completed_process(download_tasks)
                time.sleep(1)  # wait a little bit, avoid request data at the same time.

    # wait until all task complete
    while True:
        job_count = basic.alive_process_count(download_tasks)
        if job_count > 0:
            print(machine_name, datetime.now(), 'wait until all task are completed, alive task account: %d ' % job_count)
            time.sleep(30)  #
        else:
            break

    if max_task_count != 1:
        if os.path.isfile(save_file_list_txt):
            tmp_list = io_function.read_list_from_txt(save_file_list_txt)
            saved_file_list.extend(tmp_list)
        else:
            # save the list
            io_function.save_list_to_txt(save_file_list_txt,saved_file_list)

    # checking file count
    if item_count*len(bands_to_save) != len(saved_file_list):
        raise ValueError(f'STAC downloading for polygon {ext_id} is not completed')


    return True


def download_dem_stac(client,collection_id,extent_polys_latlon,extent_polys_prj , output_dir, date_start,date_end,pre_name,
                      poly_ids=None,save_prj_code=3413,b_unique_grid=False, out_res=None):

    # download data through STAC, not need to unpack
    b_save_grid_id_noDEM = True
    if poly_ids is None:
        poly_ids = [idx for idx in range(len(extent_polys_latlon)) ]
        b_save_grid_id_noDEM = False    # if poly_ids is not the global unique id, then don't save it.

    for count, (idx, ext_poly_ll, poly_prj) in enumerate(zip(poly_ids, extent_polys_latlon, extent_polys_prj)):
        if b_unique_grid:
            search_save = os.path.join(output_dir, pre_name + f'_grid{idx}.gpkg')
        else:
            search_save = os.path.join(output_dir, pre_name + f'_poly_{idx}.gpkg')
        # save_dir = os.path.join(output_dir, pre_name + f'_poly_{idx+1}')

        # # for testing, output information, will quit after running this
        # get_collection_examples_meta(client,collection_id, ext_poly, idx, date_start=date_start,
        #                             date_end=date_end, search_save=search_save,save_crs_code=save_prj_code)

        res = download_dem_within_polygon(client,collection_id, ext_poly_ll,poly_prj, idx, date_start=date_start,
                                    date_end=date_end, search_save=search_save, save_dir=output_dir,save_crs_code=save_prj_code,
                                          b_unique_grid=b_unique_grid,out_res=out_res)

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
    ext_polys_latlon = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,'EPSG:4326')
    prj_code = vector_gpd.get_projection(extent_shp)
    if prj_code != 3413:
        basic.outputlogMessage(f'Map projection of {extent_shp} is not EPSG:3413, will convert to it')
        ext_polys_3413 = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:3413')
    else:
        ext_polys_3413 = vector_gpd.read_polygons_gpd(extent_shp)


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
    if len(ext_polys_latlon) < 1:
        raise ValueError('No polygons in %s'%extent_shp)
    else:
        basic.outputlogMessage('read %d extent polygons in %s for downloading using STAC'%(len(ext_polys_latlon),extent_shp))

    date_start = options.date_start
    date_end = options.date_end

    # projection for ArcticDEM products
    save_prj_code = 3413
    out_resolution = options.out_res

    client = pystac_client.Client.open(pgc_stac)

    # b_arcticDEM_tile = False
    # collection_client = get_collection_for_search(collection_id='arcticdem-strips-s2s041-2m')
    # print(collection_client)

    download_dem_stac(client, collection_id, ext_polys_latlon, ext_polys_3413, reg_tif_dir, date_start, date_end, pre_name,
                      poly_ids=grid_id_list, save_prj_code=save_prj_code,b_unique_grid=b_unique_grid,out_res=out_resolution)




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

    parser.add_option("-r", "--out_res",
                      action="store", dest="out_res", type=float, default=2.0,
                      help="the spatial resolution for output (2.0)")

    parser.add_option("-p", "--max_process_num",
                      action="store", dest="max_process_num", type=int, default=4,
                      help="the maximum number of processes for downloading in parallel")

    parser.add_option("-m", "--b_arcticDEM_mosaic",
                      action="store_true", dest="b_arcticDEM_mosaic",default=False,
                      help="if set, will download the mosaic of ArcticDEM, other wise, download the strips")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

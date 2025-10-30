#!/usr/bin/env python
# Filename: s2_download_gee.py 
"""
introduction: download sentinel-2 imagery from Google Earth Engine for large regions

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 11 July, 2023
"""
import sys,os
from optparse import OptionParser
import time
from datetime import datetime, timezone
from datetime import timedelta

import multiprocessing
from multiprocessing import Pool
# import numpy.array_api
import numpy as np
import re

# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools
from collections import defaultdict

# input earth engine, used: ~/programs/anaconda3/envs/gee/bin/python, or change to ee by "source activate gee"
# need shapely, geopandas, gdal
import ee
import random

parameters_list = None


from gee_common import export_one_imagetoDrive,wait_all_task_finished, maximum_submit_tasks, active_task_count, \
    shapely_polygon_to_gee_polygon, reproject,directly_save_image_to_local

# image specification (more in ChangeDet_DL/dataTools/get_timelapse_img_gee.py)
# on Feb 6, update COPERNICUS/S2_SR to COPERNICUS/S2_SR_HARMONIZED (https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR)
img_speci = {'sentinel2_rgb_sr': {'product': 'COPERNICUS/S2_SR_HARMONIZED', 'bands': ['B4', 'B3', 'B2'], 'res': 10},
            # B8 is the NIR band
            'sentinel2_Nrgb_sr': {'product': 'COPERNICUS/S2_SR_HARMONIZED', 'bands': ['B4', 'B3', 'B2','B8'], 'res': 10},
            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2
            'landsat9_rgb':{'product':'LANDSAT/LC09/C02/T1_L2', 'bands':['SR_B4', 'SR_B3', 'SR_B2'], 'res':30},
            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_TOA
            'landsat9_pan':{'product':'LANDSAT/LC09/C02/T1_TOA', 'bands':['B8'], 'res':15},

            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2
            'landsat8_rgb':{'product':'LANDSAT/LC08/C02/T1_L2', 'bands':['SR_B4', 'SR_B3', 'SR_B2'], 'res':30},
             # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2
            'landsat8_Nrgb':{'product':'LANDSAT/LC08/C02/T1_L2', 'bands':['SR_B5','SR_B4', 'SR_B3', 'SR_B2'], 'res':30},
            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_TOA
             'landsat8_pan':{'product':'LANDSAT/LC08/C02/T1_TOA', 'bands':['B8'], 'res':15},

             # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_TOA
            'landsat7_pan':{'product':'LANDSAT/LE07/C02/T1_TOA', 'bands':['B8'], 'res':15},
            'landsat7_Nrgb':{'product':'LANDSAT/LE07/C02/T1_TOA', 'bands':['B4','B3','B2','B1'], 'res':30},
            }


def initialize_ee(project):
    # avoid too many rapid requests to Google's API when parallel call this function
    # MaxRetryError: HTTPSConnectionPool(host='oauth2.googleapis.com', port=443): Max retries exceeded with url:
    # /token (Caused by SSLError(SSLError(1, '[SSL] record layer failure (_ssl.c:2559)')))"
    # wait random 0.1 to 2 seconds

    # random_wait = np.random.randint(1, 5)     # alway genreated when call in the parallel function

    # Generate a random wait time between 0.1 and 2.0 seconds
    random_wait = random.uniform(0.1, 10)

    time.sleep(random_wait)
    print(f'Initializing Earth Engine with project: {project}, please wait {random_wait:.1f} seconds')
    ee.Initialize(project=project)

def get_product_info(product):
    if "S2" in product:
        product_info = product.split('/')[-1]   # will get S2_SR_HARMONIZED
    elif "LANDSAT" in product:
        product_info = product.split('/')[-3] + '_' +product.split('/')[-2] + '_' + product.split('/')[-1]  # will get LC09_C02_T1_L2
    else:
        product_info = product.replace('/','_')
    return product_info

def get_export_dir_name(region_name,product, start_date,end_date):
    # cannot have a sub folder in Google Drive.
    date_range_str = re.sub(r'\D','',start_date) + '_' + re.sub(r'\D','',end_date) # only keep digits
    product_info = get_product_info(product)
    export_dir = region_name + '_' + product_info + '_' +  date_range_str + '_images'
    return export_dir

def get_save_file_name(region_name, product, ext_id, start_date,end_date, b_not_mosaic):
    product_info = get_product_info(product)
    date_range_str = re.sub(r'\D', '', start_date) + '_' + re.sub(r'\D', '', end_date)
    if b_not_mosaic:
        # save_file_name = region_name + '_' + product_info + f'_img{ext_id}'
        save_file_name = f'img{ext_id}'        # avoid file name is too long
    else:
        save_file_name = region_name + '_' + product_info + '_' + date_range_str + f'_grid{ext_id}'

    return save_file_name


def to_toUint16(image, product):
    # if omit "toUint16()" , it output float64
    if 'S2' in product:
        image = image.toUint16()
    elif 'LANDSAT' in product:
        # for Landsat, the top of atmosphere reflectance is scaled by a factor of 10000
        if "T1_TOA" in product:
            image = image.multiply(10000).toUint16()
        else:
            image = image.toUint16()
    else:
        raise ValueError('%s not supported yet in to_toUint16()'%(product))
    return image

def remove_downloaded_tasks(region_name, product, start_date, end_date,extent_polygons, extent_ids, b_not_mosaic):
    # remove already donwloaded IDs before run parallel, to improve efficiency
    export_dir = get_export_dir_name(region_name, product, start_date, end_date)

    remain_extent_polygons = []
    remain_extent_ids = []
    for idx, (extent, ext_id) in enumerate(zip(extent_polygons, extent_ids)):
        save_file_name = get_save_file_name(region_name, product, ext_id, start_date, end_date, b_not_mosaic)
        local_record = os.path.join(export_dir, save_file_name + '.submit')
        if os.path.isfile(local_record):
            print(f'{idx+1}/{len(extent_polygons)} task: {local_record} already be submitted to GEE or downloaded, skip')
        else:
            remain_extent_polygons.append(extent)
            remain_extent_ids.append(ext_id)

    return remain_extent_polygons, remain_extent_ids

def gee_download_images(region_name,start_date, end_date, ext_id, extent, product, resolution, projection,
                        band_names, cloud_cover_thr=0.3, crop=False, b_vis=False, wait_all_finished=True,b_save2local=False,
                        b_not_mosaic=False,max_download_count=3):

    start = ee.Date(start_date)  # '%Y-%m-%d'
    finish = ee.Date(end_date)

    polygon_bound = extent
    # cloud_cover = 'CLOUDY_PIXEL_PERCENTAGE' #  for Sentinel-2
    if 'S2' in product:
        cloud_cover = 'CLOUDY_PIXEL_PERCENTAGE'
    elif 'LANDSAT' in product:
        cloud_cover = 'CLOUD_COVER'
    else:
        raise ValueError('%s not supported yet in cloud mask'%(product))

    # got error on Narval: requests.exceptions.SSLError: None: Max retries exceeded with url: /token (Caused by None)
    # the internet there is very good, send too many request during a short time, got rejected

    export_dir = get_export_dir_name(region_name,product, start_date,end_date)
    save_file_name = get_save_file_name(region_name, product, ext_id, start_date,end_date, b_not_mosaic)
    # print('export_dir:',export_dir)
    # print('save_file_name:', save_file_name)
    # sys.exit()

    if b_vis:
        save_file_name = save_file_name + '_8bit'
    # checking file existence before the query.
    save_file_path = os.path.join(export_dir, save_file_name + '.tif')
    local_record = os.path.join(export_dir, save_file_name + '.submit')
    # if b_save2local:
    #     if os.path.isfile(save_file_path):
    #         print('%s already exists, skipping downloading' % save_file_path)
    #         return False
    # else:
    if os.path.isfile(local_record):
        print('task %s already be submitted to GEE or downloaded, skip' % local_record)
        return False


    filtercollection = ee.ImageCollection(product). \
        filterBounds(polygon_bound). \
        filterDate(start, finish). \
        filter(ee.Filter.lt(cloud_cover, cloud_cover_thr * 100)). \
        sort(cloud_cover, True)

    # # map(cloud_mask). \
    # filter(ee.Filter.calendarRange(month_range[0], month_range[1], 'month')). \

    # check count  # getInfo can get python number (not ee.Number)
    count = filtercollection.size().getInfo()
    if count < 1:
        basic.outputlogMessage(f'No results for the polygon (id: {ext_id} ) with {product}')
        return False

    print(f'Find {count} image for the polygon (id: {ext_id} )' )


    # save some task record in the local folder
    if os.path.isdir(export_dir) is False:
        io_function.mkdir(export_dir)

    valid_pixel_percent_thr = 60

    # first_image = ee.Image(filtercollection.first()).select(band_names)
    # # reproject?
    # # default is UTM, reproject to the projection of extent_shp
    # first_image = reproject(first_image,projection,resolution)
    if b_not_mosaic:
        # does not create mosiac, but download many images
        ## Make a list of images
        img_list = filtercollection.toList(filtercollection.size())

        # export all
        n = 0
        tasklist = []
        while True:
            try:
                dtype = np.uint16
                an_img = ee.Image(img_list.get(n)).select(band_names)
                # image_info = an_img.getInfo()
                # print('image_info',image_info)
                time_start = an_img.get('system:time_start').getInfo()
                acquisition_time = datetime.fromtimestamp(time_start / 1000, tz=timezone.utc).strftime('%Y%m%d-%H%M%S')

                an_img = reproject(an_img, projection, resolution)
                if b_vis:
                    # s2_mosaic_rgb_8bit = mosaic.visualize(bands=band_names, min=0, max=2000)
                    an_img = an_img.visualize(bands=band_names, min=0, max=2000)
                    dtype = np.uint8

                # image_info = an_img.getInfo()
                # print('image_info',image_info)

                # change file name3
                # time_start = an_img.get('system:time_start').getInfo()  # not available
                # acquisition_time = datetime.fromtimestamp(time_start / 1000, tz=timezone.utc).strftime('%Y%m%d-%H%M%S')
                save_file_path_ii = io_function.get_name_by_adding_tail(save_file_path,f'm{n}_{acquisition_time}')
                save_file_name_ii = os.path.basename(save_file_path_ii)

                if b_save2local:
                    an_img_array = an_img.sampleRectangle(extent, defaultValue=0)
                    an_img_features = an_img_array.getInfo()  # the actual download
                    b_saved = directly_save_image_to_local(save_file_path_ii, dtype, an_img, an_img_features,
                                                 valid_pixel_percent_thr=valid_pixel_percent_thr)
                    if b_saved is False:
                        n += 1
                        continue

                else:
                    # print('testing:',save_file_name_ii)
                    task = export_one_imagetoDrive(an_img, export_dir, save_file_name_ii, extent, resolution,
                                                   wait2finished=wait_all_finished)
                    tasklist.append(task)

                # save a record in the local dir
                with open(local_record, 'w') as f_obj:
                    f_obj.writelines('downloaded or submitted to GEE on %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

                n += 1
                # print('%s: Start %dth task to download images covering %dth polygon' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n, polygon_idx))
                if n >= max_download_count:
                    break
            except Exception as e:
                error = str(e).split(':')
                if error[0] == 'List.get':
                    break
                else:
                    raise e

        # wait all task filished
        if len(tasklist) > 0:
            wait_all_task_finished(tasklist)


    else:

        # if omit "toUint16()" , it output float64
        dtype = np.uint16
        mosaic = filtercollection.select(band_names,band_names).median()  # create mosaic using median values
        mosaic = to_toUint16(mosaic, product) # convet to uint16
        mosaic = reproject(mosaic, projection, resolution)
        # Error: Image.visualize: Expected a string or list of strings for field 'bands'. (Error code: 3)
        # rgbVis = {'min': 0, 'max': 2000, 'bands': band_names }
        # s2_mosaic_rgb_8bit = mosaic.visualize(rgbVis)
        if b_vis:
            # s2_mosaic_rgb_8bit = mosaic.visualize(bands=band_names, min=0, max=2000)
            mosaic = mosaic.visualize(bands=band_names, min=0, max=2000)
            dtype = np.uint8


        if b_save2local:
            mosaic_array = mosaic.sampleRectangle(extent, defaultValue=0)
            mosaic_features = mosaic_array.getInfo()  # the actual download
            directly_save_image_to_local(save_file_path,dtype,mosaic,mosaic_features)
            # save a record in the local dir
            with open(local_record, 'w') as f_obj:
                f_obj.writelines('Downloaded on %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        else:
            # only export first image
            # export_one_imagetoDrive(first_image,export_dir,polygon_idx,crop_region, img_speci['res'])
            task = export_one_imagetoDrive(mosaic, export_dir, save_file_name, extent, resolution, wait2finished=wait_all_finished)
            # task = export_one_imagetoDrive(s2_mosaic_rgb_8bit, export_dir, save_file_name+'_8bit', extent, resolution, wait2finished=wait_all_finished)

            # save a record in the local dir
            with open(local_record, 'w') as f_obj:
                f_obj.writelines('submitted to GEE on %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

            return task




def parallel_gee_download_images_to_local(idx, total_count, region_name,start_date, end_date, ext_id, extent, product, resolution, projection,
                        bands, cloud_cover_thr=0.3, crop=False, b_vis=False, wait_all_finished=True,b_save2local=True,
                                          b_not_mosaic=False,max_download_count=3, gee_project=None):
    # ee.Initialize(project=gee_project)
    print(f'{idx + 1}/{total_count} Downloading images for a polygon (id: {ext_id}) ')

    extent_gee = shapely_polygon_to_gee_polygon(extent)

    gee_download_images(region_name, start_date, end_date, ext_id, extent_gee, product, resolution,
                               projection, bands, cloud_cover_thr=cloud_cover_thr,
                               crop=crop, b_vis=b_vis, wait_all_finished=wait_all_finished,b_save2local=b_save2local,
                                b_not_mosaic=b_not_mosaic,max_download_count=max_download_count)

def gee_download_sentinel2_image(extent_shp, region_name,id_column_name, start_date, end_date,cloud_cover_thr,
                                 b_save2local=False, process_num=8,b_not_mosaic=False,max_download_count=3, gee_project=None,
                                 image_type='sentinel2_Nrgb_sr'):

    # checking input shapefile
    projection = map_projection.get_raster_or_vector_srs_info_epsg(extent_shp)

    # product = img_speci['sentinel2_rgb_sr']['product']
    # resolution = img_speci['sentinel2_rgb_sr']['res']
    # bands = img_speci['sentinel2_rgb_sr']['bands']

    product = img_speci[image_type]['product']
    resolution = img_speci[image_type]['res']
    bands = img_speci[image_type]['bands']


    b_crop = True  # crop images to input extent
    b_visualize = False   # to 8bit, for visualization

    extent_polygons = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:4326')
    extent_ids = vector_gpd.read_attribute_values_list(extent_shp, id_column_name)
    if extent_ids is None:
        extent_ids = [ str(item) for item in range(len(extent_polygons))]
    else:
        if isinstance(extent_ids[0], (float, np.floating)): # only check the first one
            extent_ids = [str(int(item)) for item in extent_ids]

    all_task_list = []
    if b_save2local:
        if process_num == 1:
            for idx, (extent, ext_id) in enumerate(zip(extent_polygons, extent_ids)):
                print(f'{idx + 1}/{len(extent_polygons)} Downloading images for a polygon (id: {ext_id})')

                extent_gee = shapely_polygon_to_gee_polygon(extent)

                gee_download_images(region_name, start_date, end_date, ext_id, extent_gee, product, resolution,
                                    projection, bands, cloud_cover_thr=cloud_cover_thr,
                                    crop=b_crop, b_vis=b_visualize, wait_all_finished=True,
                                    b_save2local=b_save2local,b_not_mosaic=b_not_mosaic,max_download_count=max_download_count)

        else:
            # save to local disks
            theadPool = Pool(process_num,initializer=initialize_ee, initargs=(gee_project,))  # multi processes ,initializer=initialize_ee()

            # remove already downloaded IDs before run parallel, to improve efficiency
            extent_polygons_2, extent_ids_2 = remove_downloaded_tasks(region_name, product, start_date, end_date, extent_polygons, extent_ids, b_not_mosaic)

            parameters_list = [(idx,len(extent_polygons_2), region_name,start_date, end_date,ext_id, extent,product, resolution, projection,
                                bands,cloud_cover_thr,b_crop,b_visualize, True, True,b_not_mosaic,max_download_count, gee_project)
                               for idx, (extent, ext_id) in enumerate(zip(extent_polygons_2, extent_ids_2))]
            results = theadPool.starmap(parallel_gee_download_images_to_local, parameters_list)  # need python3
            theadPool.close()

    else:
        # exporting to Google Drive
        for idx, (extent, ext_id) in enumerate(zip(extent_polygons, extent_ids)):

            # # for test (two grids cover Willow River)
            # if ext_id not in [1166,1241]: # 1241
            #     continue

            active_tasks = active_task_count(all_task_list)
            while active_tasks > maximum_submit_tasks:
                print('%s : %d (>%d) tasks already in the queue, wait 60 seconds' % (
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'), active_tasks, maximum_submit_tasks))
                time.sleep(60)
                active_tasks = active_task_count(all_task_list)

            print(f'{idx+1}/{len(extent_polygons)} Downloading images for a polygon (id: {ext_id})')

            extent_gee = shapely_polygon_to_gee_polygon(extent)

            task = gee_download_images(region_name, start_date, end_date, ext_id, extent_gee, product, resolution,
                                       projection, bands, cloud_cover_thr=cloud_cover_thr,
                                       crop=b_crop, b_vis=b_visualize, wait_all_finished=False,
                                       b_not_mosaic=b_not_mosaic,max_download_count=max_download_count)

            if task is False:
                continue

            all_task_list.append(task)

        # wait until all finished
        wait_all_task_finished(all_task_list)

    if b_save2local:
        # save downloaded images list

        # save labels for each download images if "class_int" exists
        class_int_list = vector_gpd.read_attribute_values_list(extent_shp, 'class_int')

        print('To save image and labels (if class_int exists) for each downloaded image ')

        # get image folder (the same as in "gee_download_images")
        export_dir = get_export_dir_name(region_name, product, start_date, end_date)

        save_tif_list = []
        save_tif_list_txt = export_dir+'_list.txt'
        total_count = len(extent_ids)

        if os.path.isdir(export_dir):
            if class_int_list is None:
                class_int_list_str = [' ']*len(extent_ids)  # class label is empty
            else:
                class_int_list_str = [str(item) for item in class_int_list]

            ### this is too slow, each time, need to get_file_list_by_pattern ##
            ##
            # # for poly_id, c_label in zip(extent_ids,class_int_list_str):
            # for idx, (poly_id, c_label) in enumerate(zip(extent_ids, class_int_list_str), start=1):
            #     if idx%100 ==0:
            #         print(datetime.now(),f"Progress: {idx/total_count:.2f}% ({idx}/{total_count})")
            #     # class_int_list
            #     save_file_name = region_name + '_' + product_info[-1] + '_img%d' % poly_id
            #     save_file_list = io_function.get_file_list_by_pattern(export_dir,save_file_name + '*.tif')
            #     for s_file in save_file_list:
            #         save_tif_list.append(s_file+' '+c_label)
            ##################################################################


            all_tif_files = io_function.get_file_list_by_pattern(export_dir, '*.tif')
            file_groups = defaultdict(list)
            for file_path in all_tif_files:
                file_basename = os.path.basename(file_path)
                poly_id = re.findall(r'img\d+', file_basename)[0][3:]
                file_groups[poly_id].append(file_basename)
            for idx, (poly_id, c_label) in enumerate(zip(extent_ids, class_int_list_str), start=1):
                if idx % 1000 == 0:
                    print(datetime.now(), f"Progress: {idx / total_count:.2%} ({idx}/{total_count})")
                for s_file in file_groups.get(poly_id, []):  # Use .get() to avoid KeyError
                    save_tif_list.append(f"{s_file} {c_label}")


            io_function.save_list_to_txt(save_tif_list_txt ,save_tif_list)

        else:
            print(f'Warning, folder: {export_dir} does not exist')


def main(options, args):


    # all images will save to Google Drive first
    # currently, manually download them from Google Drive

    extent_shp = args[0]
    io_function.is_file_exist(extent_shp)

    region_name = options.region_name
    if region_name is None:
        region_name = os.path.splitext(os.path.basename(extent_shp))[0]
    cloud_cover_thr = options.cloud_cover

    id_column_name = options.id_column
    b_save2local = options.b_save2local
    process_num = options.process_num

    # global max_download_count, b_not_mosaic
    # max_download_count = 3
    # b_not_mosaic = False
    max_download_count = options.max_count
    b_not_mosaic = options.b_not_mosaic
    image_type = options.image_type


    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # for each computer, need to run "earthengine authenticate" first.
    # ee.Initialize()
    # after Oct 2024, GEE need to link to a Google project
    google_cloud_project = "gee-project-99319"
    if process_num == 1:
        ee.Initialize(project=google_cloud_project)
    else:
        print('Will Initialize gee in the parallel Pool')


    start_date, end_date = options.start_date, options.end_date
    gee_download_sentinel2_image(extent_shp,region_name, id_column_name,start_date, end_date,cloud_cover_thr,
                                 b_save2local=b_save2local, process_num=process_num,b_not_mosaic=b_not_mosaic,
                                 max_download_count=max_download_count,gee_project=google_cloud_project,
                                 image_type=image_type)




if __name__ == '__main__':

    usage = "usage: %prog [options] extent_shp "
    parser = OptionParser(usage=usage, version="1.0 2023-07-11")
    parser.description = 'Introduction: get Sentinel-2 imagery from Google Earth Engine for a study area '
    parser.add_option("-s", "--start_date", default='2020-06-01',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2020-06-01")
    parser.add_option("-e", "--end_date", default='2020-09-30',
                      action="store", dest="end_date",
                      help="the end date for inquiry, with format year-month-day, e.g., 2020-09-30")

    # parser.add_option("-d", "--dates_list_txt",
    #                   action="store", dest="dates_list_txt",
    #                   help="a list a date, if this is set, the start_date and end_date will be ignored")

    parser.add_option("-c", "--cloud_cover",
                      action="store", dest="cloud_cover", type=float, default = 0.1,
                      help="the could cover threshold, only accept images with cloud cover less than the threshold")

    parser.add_option("-n", "--region_name",
                      action="store", dest="region_name",
                      help="the name of the area, which we download SAR data for")

    parser.add_option("-i", "--id_column",
                      action="store", dest="id_column",default='id',
                      help="the name of unique ID column")

    parser.add_option('-r',"--resolution", action='store', dest='resolution',default=10, type=float,
                      help='the resolution of output raster')

    parser.add_option('-p',"--projection_epsg", action='store', dest='projection_epsg',
                      help='the projection of output raster')

    parser.add_option("", "--b_save2local",
                      action="store_true", dest="b_save2local", default=False,
                      help="if set, will save images to local disk, not exporting to Google Drive (for small images less than 30 MB)")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=8,
                      help="number of processes to download images to local disks")

    parser.add_option("", "--max_count",
                      action="store", dest="max_count", type=int,default=3,
                      help="the maximum count of images for a polygon to download, default is 3")
    parser.add_option("", "--b_not_mosaic",
                      action="store_true", dest="b_not_mosaic", default=False,
                      help="if set, it will not create mosaic for each polygon, but download all available images)")

    parser.add_option("", "--image_type",
                      action="store", dest="image_type", default='sentinel2_Nrgb_sr',
                      help="the image types want to download (contain product, bands, resolution), more find 'img_speci' in the head of this script")

    # multiprocessing.set_start_method('spawn')
    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

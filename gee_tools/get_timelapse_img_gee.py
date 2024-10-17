#!/usr/bin/env python
# Filename: get_timelapse_img_gee 
"""
introduction: get time-lapse images using Google Earth Engine (GEE)



authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 08 December, 2019
"""

import sys,os
from optparse import OptionParser
from datetime import datetime


# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import basic_src.basic as basic
import basic_src.io_function as io_function

import vector_gpd
from shapely.geometry import mapping        # transform to GeJSON format

import math
import time

# input earth engine, used: ~/programs/anaconda3/envs/gee/bin/python, or change to ee by "source activate gee"
# need shapely, geopandas, gdal
import ee

# re-project
from functools import partial
import pyproj
from shapely.ops import transform

shp_polygon_projection = None
month_range = [7,8]
max_download_count  = None

# image specification
img_speci = {# https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_L2
            'landsat9_rgb':{'product':'LANDSAT/LC09/C02/T1_L2', 'bands':['B4', 'B3', 'B2'], 'res':30},
            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC09_C02_T1_TOA
            'landsat9_pan':{'product':'LANDSAT/LC09/C02/T1_TOA', 'bands':['B8'], 'res':15},
            # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_SR
            'landsat8_rgb':{'product':'LANDSAT/LC08/C01/T1_SR', 'bands':['B4', 'B3', 'B2'], 'res':30},
              # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1_TOA
             'landsat8_pan':{'product':'LANDSAT/LC08/C01/T1_TOA', 'bands':['B8'], 'res':15},
              # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_SR
              'landsat7_rgb':{'product':'LANDSAT/LE07/C02/T1_TOA', 'bands':['B4', 'B3', 'B2'], 'res':30},
              # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1_TOA
              'landsat7_pan':{'product':'LANDSAT/LE07/C02/T1_TOA', 'bands':['B8'], 'res':15},
              # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1_SR
               'landsat5_rgb':{'product':'LANDSAT/LT05/C01/T1_SR', 'bands':['B3', 'B2', 'B1'], 'res':30},
              # https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2
                'sentinel2_rgb':{'product':'COPERNICUS/S2', 'bands':['B4', 'B3', 'B2'], 'res':10},
              # https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR#description
                 'sentinel2_rgb_sr': {'product': 'COPERNICUS/S2_SR', 'bands': ['B4', 'B3', 'B2'], 'res': 10}
            }

def get_projection_proj4(geo_file):
    import basic_src.map_projection as map_projection
    return map_projection.get_raster_or_vector_srs_info_proj4(geo_file)

def meters_to_degress_onEarth(distance):
    return (distance/6371000.0)*180.0/math.pi

def reproject_shapely_polygon(in_polygon_shapely,src_prj, new_projection):
    project = partial(
        pyproj.transform,
        pyproj.Proj(src_prj),  # source coordinate system
        pyproj.Proj(new_projection))  # destination coordinate system

    out_polygon = transform(project, in_polygon_shapely)  # apply projection
    # polygon_shapely = transform(project, polygon_shapely)
    return out_polygon

def get_image_name(image_info,product):

    if image_info['type'] != 'Image':
        raise ValueError('input is a information of an image')
    # elif 'LC08/C01/T1_SR' in product:
    #     return maskL578clouds_SR
    # elif 'LE07/C01/T1_SR' in product:
    #     return maskL578clouds_SR
    # elif 'LT05/C01/T1_SR' in product:
    #     return maskL578clouds_SR
    # elif 'LC08/C01/T1_TOA' in product:
    #     return maskL578clouds
    # elif 'LE07/C01/T1_TOA' in product:
    #     return maskL578clouds
    elif 'S2' in product:
        name = 'SPACECRAFT_NAME'
        satellite = image_info['properties'][name]
        # return name
    elif 'T1_SR' in product:
        name = 'SATELLITE'
        satellite = image_info['properties'][name]
        # return name
    elif 'LC08/C01/T1_TOA' in product:
        satellite = 'Landsat_8_TOA'
    elif 'LE07/C01/T1_TOA' in product:
        satellite = 'Landsat_7_TOA'
    elif 'LE07/C02/T1_TOA' in product:
        satellite = 'Landsat_7_TOA'
    else:
        raise ValueError('%s not supported yet'%(product))

    # satellite = image_info['properties'][name]

    #  timestamp in miliseconds (like in JavaScript), but fromtimestamp() expects Unix timestamp, in seconds
    fromtimestamp = image_info['properties']['system:time_start'] / 1000  #
    acquire_time = datetime.fromtimestamp(fromtimestamp)
    acquire_time_str = acquire_time.strftime('%Y-%m-%d')
    bands_str = '_'.join([item['id'] for item in image_info['bands']])

    out_name = satellite +'_' + bands_str + '_' +acquire_time_str
    return out_name

def get_crop_region(polygon_shapely,img_crs, buffer_size):


    # re-projection if necessary (no need to reproject, GEE only accept lat lon to crop)
    # but need to do conduct buffer operation in XY coordinates

    if img_crs != 'epsg:4326':
        polygon_prj = reproject_shapely_polygon(polygon_shapely,shp_polygon_projection, img_crs)
        expansion_polygon = polygon_prj.buffer(buffer_size)
        # projection back
        expansion_polygon = reproject_shapely_polygon(expansion_polygon,img_crs, shp_polygon_projection)
    else:
        # change buffer area from meters to degree
        buffer_size = meters_to_degress_onEarth(buffer_size)
        expansion_polygon = polygon_shapely.buffer(buffer_size)

    # defined a crop region (# this cause an  Internal error )
    # expansion_polygon_json = mapping(expansion_polygon)
    # crop_region = ee.Geometry(expansion_polygon_json)

    polygon_env = expansion_polygon.envelope
    x, y = polygon_env.exterior.coords.xy
    crop_region_list = [min(x),min(y), max(x), max(y)]
    crop_region = ee.Geometry.Rectangle(crop_region_list)

    return crop_region

# /**
#  * Function to mask clouds using the Sentinel-2 QA band
#  * @param {ee.Image} image Sentinel-2 image
#  * @return {ee.Image} cloud masked Sentinel-2 image
#  */
def maskS2clouds(image):
    qa = image.select('QA60')
    #// Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11
    #// Both flags should be set to zero, indicating clear conditions.
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    # return image.updateMask(mask).divide(10000)  # divide (10000) would lose some properties
    return image.updateMask(mask)

def maskL578clouds_SR(image):
    qa = image.select('pixel_qa')
    # The cloud layer is represented as the 5 place, the cloud layer confidence is 6-7, and the cloud shadow is the 3 place
    # Select the pixels that have clouds and the cloud confidence is medium, and are covered by cloud shadows.
    cloudBitMask = 1 << 5
    cloudConfi = 1 << 7
    cloudShadow = 1 <<3
    mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cloudShadow).eq(0))
    #.Or(qa.bitwiseAnd(cloudConfi).eq(0))
    return image.updateMask(mask)
    # cloud = qa.bitwiseAnd(1 << 3).And (qa.bitwiseAnd(1 << 5)).Or (qa.bitwiseAnd(1 << 7))
    # # remove boundary pixels
    # mask = image.mask().reduce(ee.Reducer.min())
    #
    # return image.updateMask(cloud.Not()).updateMask(mask)

def maskL578clouds(image):
    qa = image.select('BQA')
    # The cloud layer is represented as the 4 place, the cloud layer confidence is 5-6, and the cloud shadow is the 7-8 place
    # Select the pixels that have clouds and the cloud confidence is medium, and are covered by cloud shadows.
    cloud = qa.bitwiseAnd(1 << 4).And(qa.bitwiseAnd(1 << 6)).Or(qa.bitwiseAnd(1 << 8))
    # remove boundary pixels
    mask = image.mask().reduce(ee.Reducer.min())

    return image.updateMask(cloud.Not()).updateMask(mask)
    # cloudShadowBitMask = 1 << 3
    # cloudsBitMask = 1 << 5
    # #// 3 cloudshadow, 5 cloud
    # qa = image.select('pixel_qa')
    # #// Both flags should be set to zero
    # mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And (qa.bitwiseAnd(cloudsBitMask).eq(0));
    # #//
    # return image.updateMask(mask)
    # #.divide(10000).select("B[0-9]*").copyProperties(image, ["system:time_start"]);

def get_cloud_mask_function(product):

    if 'S2' in product:
        return maskS2clouds
    elif 'LC08/C01/T1_SR' in product:
        return maskL578clouds_SR
    elif 'LE07/C01/T1_SR' in product:
        return maskL578clouds_SR
    elif 'LT05/C01/T1_SR' in product:
        return maskL578clouds_SR
    elif 'LC08/C01/T1_TOA' in product:
        return maskL578clouds
    elif 'LE07/C01/T1_TOA' in product:
        return maskL578clouds
    else:
        raise ValueError('%s not supported yet',product)

def export_one_imagetoDrive(select_image, save_folder,polygon_idx, crop_region, res, product, wait2finished=True):

    # map(cloud_mask).
    image_info = select_image.getInfo()
    save_file_name = get_image_name(image_info,product) + '_poly_%d'%polygon_idx

    local_record_folder = 'local_records'
    if os.path.isdir(local_record_folder) is False:
        io_function.mkdir(local_record_folder)

    local_record = os.path.join(os.path.join(local_record_folder, save_file_name+'.submit'))
    if os.path.isfile(local_record):
        print('task %s already be submitted to GEE, skip' % local_record)
        return False

    task = ee.batch.Export.image.toDrive(image=select_image,
                                         region=crop_region,
                                         description=save_file_name,
                                         folder=save_folder,
                                         maxPixels = 1e10,      # default is 1e8
                                         scale=res)
    # region=crop_region,

    task.start()

    # save a record in the local dir
    with open(local_record, 'w') as f_obj:
        f_obj.writelines('submitted to GEE on %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    # print(task.status())
    # print(ee.batch.Task.list())
    if wait2finished:
        import time
        print('%s: Start transferring %s to Drive..................' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),save_file_name))
        while task.active():
            print('%s: Transferring %s to Drive_..................'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),save_file_name))
            time.sleep(20)
        print('%s: Done with the Export to the Drive_'%datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        return True
    else:
        return task

def wait_all_task_finished(all_tasks, polygon_idx):
    all_tasks = [ item for item in all_tasks if isinstance(item, bool) is False]
    all_count = len(all_tasks)

    finished_count = 0
    # print('%s: Start transferring %d images covering %d th polygon to Drive..........' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),all_count,polygon_idx))
    while finished_count < all_count:
        finished_count = 0
        for task in all_tasks:
            if task.active():
                continue
            else:
                finished_count += 1
        time.sleep(20)
        print('%s: Transferring %d images to Drive, finished %d ones..........' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),all_count, finished_count) )

    print('%s: Done with the Export to the Drive'%datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    return True

# quick test
def environment_test():
    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # ee.Initialize()       # initialize in previous place
    image = ee.Image('srtm90_v4')
    print(image.getInfo())

# def test_download():
#     ee.Initialize()
#     from pyproj import Proj, transform
#     image = ee.Image(
#         'LANDSAT/LE07/C01/T1/LE07_149035_20010930')  # change the name of this image if testing on another image is required
#     projI = image.select(
#         'B1').projection().getInfo()  # Exporting the whole image takes time, therefore, roughly select the center region of the image
#     crs = projI['crs']
#     outProj = Proj(init='epsg:4326')
#     inProj = Proj(init=crs)
#     kk = projI['transform']
#     print(projI)
#     lon1, lat1 = transform(inProj, outProj, kk[2] + (2000 * 30), kk[5] - (2000 * 30))
#     lon2, lat2 = transform(inProj, outProj, kk[2] + (5000 * 30),
#                            kk[5] - (5000 * 30))  # change these coordinates to increase or decrease image size
#
#     bounds = [lon1, lat2, lon2, lat1]
#     print(bounds, lon1, lat1)
#
#     # imageL7=imageL7.select(['B1','B2','B3','B4','B5','B6_VCID_2','B7','B8'])
#     imageL7 = image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_2', 'B7', 'B8'])  # bug fix -Lingcao
#     geometry = ([lon1, lat1], [lon1, lat2], [lon2, lat2], [lon2, lat1])
#     config = {
#         'description': 'Landsat07image',
#         'region': geometry,
#         'scale': 15,  # the image is exported with 15m resolution
#         'fileFormat': 'GeoTIFF'
#     }
#     #        'maxPixels': 1e12
#
#     exp = ee.batch.Export.image.toDrive(imageL7, **config);
#
#     exp.start()  # It takes around 5-10 minutes for 6000 * 6000 * 8 image to be exported
#     print(exp.status())
#     print(ee.batch.Task.list())
#     import time
#     while exp.active():
#         print('Transferring Data to Drive..................')
#         time.sleep(30)
#     print('Done with the Export to the Drive')

def gee_download_time_lapse_images(start_date, end_date, cloud_cover_thr, img_speci, polygon_shapely, polygon_idx, save_dir, buffer_size):
    '''
    python time lapse image using google earth engine
    :param start_date: start date, e.g., 2019-12-31
    :param end_date: e.g., 2019-12-31
    :param cloud_cover_thr: e.g., 0.3
    :param img_speci: defined product, bands, resolution e.g., 'LANDSAT/LC08/C01/T1' or 'COPERNICUS/S2 or COPERNICUS/S2_SR'
    :param polygon_bound: the extent
    :param polygon_idx: the index of the polygon in the original folder
    :param buffer_size: buffer_size
    :return:
    '''

    start = ee.Date(start_date)  # '%Y-%m-%d'
    finish = ee.Date(end_date)

    # get polygon bounding box (use the envelope)
    polygon_env = polygon_shapely.envelope
    x, y = polygon_env.exterior.coords.xy
    polygon_bound = ee.Geometry.Polygon([[x[0],y[0]],
                                         [x[1], y[1]],
                                         [x[2], y[2]],
                                         [x[3], y[3]]])

    # cloud_mask = get_cloud_mask_function(img_speci['product'])

    product = img_speci['product']
    print(product)

    if 'S2' in product:
        cloud_cover = 'CLOUDY_PIXEL_PERCENTAGE'
        # return cloud_cover
    elif 'LANDSAT' in product:
        cloud_cover = 'CLOUD_COVER'
        # return cloud_cover
    else:
        raise ValueError('%s not supported yet in cloud mask'%(product))
    filtercollection = ee.ImageCollection(img_speci['product']). \
        filterBounds(polygon_bound). \
        filterDate(start, finish). \
        filter(ee.Filter.calendarRange(month_range[0], month_range[1], 'month')). \
        filter(ee.Filter.lt(cloud_cover, cloud_cover_thr*100)). \
        sort(cloud_cover, True)   # True: ascending

    # # map(cloud_mask). \

    # check count  # getInfo can get python number (not ee.Number)
    count = filtercollection.size().getInfo()
    if count < 1:
        basic.outputlogMessage('No results for %dth polygon with %s'%(polygon_idx,str(img_speci)))
        return False

    if max_download_count is None:
        print('Image count %d, try to download all'%count)
    else:
        print('Image count %d, try to %s images (maximum), less could cover first ' %(count, str(max_download_count)))
    # print(filtercollection)                 # print serialized request instructions
    # print(filtercollection.getInfo())       # print object information
    # polygon_idx
    # export_dir = 'gee_saved' # os.path.join(save_dir,'sub_images_of_%d_polygon'%polygon_idx)
    # export_dir = ee.String(os.path.join(save_dir,'images_of_%d_polygon'%polygon_idx))
    # cannot have a sub folder in Google Drive.
    export_dir = 'images_of_%d_polygon'%polygon_idx

    # Make a list of images
    img_list = filtercollection.toList(filtercollection.size())

    first_image = ee.Image(filtercollection.first()).select(img_speci['bands'])

    projI = first_image.select(img_speci['bands'][0]).projection().getInfo()  # Exporting the whole image takes time, therefore, roughly select the center region of the image
    img_crs = projI['crs']

    crop_region = get_crop_region(polygon_shapely,img_crs, buffer_size)


    # only export first image
    # export_one_imagetoDrive(first_image,export_dir,polygon_idx,crop_region, img_speci['res'])

    # export all
    n = 0
    tasklist = []
    while True:
        try:
            img = ee.Image(img_list.get(n)).select(img_speci['bands'])
            task = export_one_imagetoDrive(img, export_dir, polygon_idx, crop_region, img_speci['res'], img_speci['product'], wait2finished=False)
            tasklist.append(task)
            n += 1
            print('%s: Start %dth task to download images covering %dth polygon'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n, polygon_idx))
            if max_download_count is not None and n >= max_download_count:
                break
        except Exception as e:
            error = str(e).split(':')
            if error[0] == 'List.get':
                break
            else:
                raise e

    # wait all task filished
    wait_all_task_finished(tasklist,polygon_idx)


    return True


def download_time_series_for_a_polygon(start_date, end_date, cloud_cover_thr, image_type, polygon_shapely, polygon_idx, save_dir, buffer_size):

    if image_type is None:
        # download all image products

        for image_type in img_speci.keys():
            gee_download_time_lapse_images(start_date, end_date, cloud_cover_thr,
                                       img_speci[image_type], polygon_shapely, polygon_idx, save_dir, buffer_size)
    else:
        if image_type not in img_speci.keys():
            raise ValueError('%s not in the key of img_speci: %s' % (image_type, str(img_speci.keys())))

        gee_download_time_lapse_images(start_date, end_date, cloud_cover_thr,
                                       img_speci[image_type], polygon_shapely, polygon_idx, save_dir, buffer_size)


def main(options, args):


    # initialize earth engine environment
    #ee.Initialize()
    # after Oct 2024, GEE need to link to a Google project
    # ee.Authenticate()
    ee.Initialize(project='gee-project-99319')

    # environment_test()  # for each computer, need to run "earthengine authenticate" first.
    # environment_test()
    # test_download()
    # return False

    polygons_shp = args[0]

    # all images will save to Google Drive first, then move them to the below folder.
    time_lapse_save_folder = args[1]  # folder for saving downloaded images

    months_range_str = options.month_range
    global month_range
    month_range = [ int(item) for item in months_range_str.split(',') ]

    global max_download_count
    max_download_count = options.max_count

    # check training polygons
    assert io_function.is_file_exist(polygons_shp)
    os.system('mkdir -p ' + time_lapse_save_folder)


    # # check these are EPSG:4326 projection
    global shp_polygon_projection
    shp_polygon_projection = get_projection_proj4(polygons_shp).strip()

    if shp_polygon_projection != '+proj=longlat +datum=WGS84 +no_defs':
        raise ValueError('Only accept %s, please check projection of %s'%(shp_polygon_projection,polygons_shp))

    # no need to convert, keep meters.
    # if shp_polygon_projection == '+proj=longlat +datum=WGS84 +no_defs':
    #     crop_buffer = meters_to_degress_onEarth(options.buffer_size)
    # else:
    #     crop_buffer = options.buffer_size
    crop_buffer = options.buffer_size

    # read polygons, not json format, but shapely format
    polygons = vector_gpd.read_polygons_json(polygons_shp, no_json=True)

    for idx, geom in enumerate(polygons):

        basic.outputlogMessage('downloading and cropping images for %dth polygon, total: %d polygon(s)' %
                               (idx, len(polygons)))

        download_time_series_for_a_polygon(options.start_date, options.end_date, options.cloud_cover,
                                   options.image_type, geom, idx, time_lapse_save_folder, crop_buffer)

        #break


    pass



if __name__ == "__main__":

    usage = "usage: %prog [options] polygon_shp save_dir"
    parser = OptionParser(usage=usage, version="1.0 2019-12-08")
    parser.description = 'Introduction: get time-lapse images using Google Earth Engine (GEE) '
    parser.add_option("-s", "--start_date",default='2016-01-01',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2016-01-01")
    parser.add_option("-e", "--end_date",default='2019-12-31',
                      action="store", dest="end_date",
                      help="the end date for inquiry, with format year-month-day, e.g., 2019-12-31")
    parser.add_option("-m", "--month_range",default='7,8',
                      action="store", dest="month_range",
                      help="selected months for inquiring images, the start and end month")
    parser.add_option("-c", "--cloud_cover",
                      action="store", dest="cloud_cover", type=float, default = 0.1,
                      help="the could cover threshold, only accept images with cloud cover less than the threshold")
    parser.add_option("-b", "--buffer_size",
                      action="store", dest="buffer_size", type=int, default = 3000,
                      help="the buffer size to crop image in meters")
    parser.add_option("-n", "--max_count",
                      action="store", dest="max_count", type=int,
                      help="the maximum count of images for a polyon to download, default is None (all available images)")
    parser.add_option("-i", "--image_type",
                      action="store", dest="image_type",
                      help="the image types want to download (contain product, bands, resolution), more find 'img_speci' in the head of this script")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    basic.setlogfile('get_timelapse_img_gee_%s.log' % str(datetime.date(datetime.now())))

    main(options, args)
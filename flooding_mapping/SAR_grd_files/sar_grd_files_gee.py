#!/usr/bin/env python
# Filename: sar_grd_files_gee.py
"""
introduction: download SAR GRD files from Google Earth Engine

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 29 May, 2021
"""


import sys,os
from optparse import OptionParser
import time
from datetime import datetime
from datetime import timedelta

# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools


gee_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'gee_tools')
sys.path.insert(0, gee_dir)
from gee_common import export_one_imagetoDrive,wait_all_task_finished, maximum_submit_tasks, active_task_count
from gee_common import environment_test, reproject
from snow_cover import get_projection_proj4,get_extent_xy_from_shp

import ee

out_projection = None
out_res = 10

def get_sentinel_sar_product_id(image_info):
    if image_info['type'] != 'Image':
        raise ValueError('The information is not Image type, but is %s'%image_info['type'])
    id_str = image_info['id'].split('/')
    image_id = id_str[2]
    return image_id

def reproject_image(image):
    # reproject to EPSG:3413 (polar scene), with resolution of 500
    reprj_img = reproject(image, out_projection, out_res)
    return reprj_img

def mask_sar_grd_edge(image):
    # mask edge
    edge = image.lt(-30.0)
    maskedImage = image.mask().And(edge.Not())
    return image.updateMask(maskedImage)

def image_float(image):
    # the sar images is float 64, convert to float 32
    return image.toFloat()


def gee_download_sentinel1_sar_grd(region_name, product, polarization, mode,orbit_pass, start_date, end_date,
                                   extent_4points, save_dir,resolution, projection,crop=False,wait_all_finished=True):
    '''
    downloand Sentinel-1 SAR grd file from Google Earth Engine
    https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD
    :param region_name:
    :param product:
    :param polarization:
    :param mode:
    :param orbit_pass:
    :param start_date:
    :param end_date:
    :param extent_4points:
    :param save_dir:
    :param resolution:
    :param projection:
    :param crop: if download individual images, set to False, not crop to extent_4points
    :param wait_all_finished:
    :return:
    '''

    res = resolution   #  resolution

    # date range
    start = ee.Date(start_date)  # '%Y-%m-%d'
    finish = ee.Date(end_date)

    # extent
    x = extent_4points[0]
    y = extent_4points[1]
    polygon_bound = ee.Geometry.Polygon([[x[0],y[0]],
                                         [x[1], y[1]],
                                         [x[2], y[2]],
                                         [x[3], y[3]]])

    sar_grd_collections = ee.ImageCollection(product). \
        filterBounds(polygon_bound). \
        filterDate(start, finish). \
        filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization)). \
        filter(ee.Filter.eq('instrumentMode', mode)). \
        filter(ee.Filter.eq('orbitProperties_pass',orbit_pass)). \
        select(polarization). \
        map(mask_sar_grd_edge). \
        map(image_float)

    # the default projection is UTM, try to reproject to EPSG:2163, but the result is wrong, so cancel reproject
    # if projection is not None:
    #     sar_grd_collections = sar_grd_collections.map(reproject_image)

    # check count  # getInfo can get python number (not ee.Number)
    count = sar_grd_collections.size().getInfo()
    if count < 1:
        basic.outputlogMessage('No %s results for extent: %s polygon, polarization: %s, mode: %s, orbit_pass: %s '
                               %(product, str(extent_4points),polarization, mode,orbit_pass))
        return False
    basic.outputlogMessage('%s SAR GRD count %d, date %s to %s, polarization: %s, mode: %s, orbit_pass: %s'
                           %(product,count,start_date,end_date,polarization, mode,orbit_pass))

    export_dir = 'S1_%s'%region_name
    # download them if not exist
    # Make a list of images
    img_list = sar_grd_collections.toList(sar_grd_collections.size()) # ee.List
    # print(img_list)
    # for idx in range(count):
    #     img = ee.Image(img_list.get(idx)) #.select(img_speci['bands'])
    #     img_info = img.getInfo()
    #     # io_function.save_dict_to_txt_json('img_%d.json'%idx,img_info)
    #     select_image = img.select(band_name)
    #     image_info = select_image.getInfo()
    #     # print(image_info)
    #     print(get_snow_product_id(image_info))


    # save some task record in the local folder
    if os.path.isdir(export_dir) is False:
        io_function.mkdir(export_dir)


    # export all
    n = 0
    tasklist = []
    while True:
        try:
            img = ee.Image(img_list.get(n)).select(polarization)
            image_info = img.getInfo()
            # print(image_info)
            save_file_name = get_sentinel_sar_product_id(image_info)
            local_record = os.path.join(os.path.join(export_dir,save_file_name))
            if os.path.isfile(local_record):
                print('task %s already submit to GEE, skip'%local_record)
                n += 1
                continue

            if crop:
                crop_region = polygon_bound
            else:
                crop_region = None
            task = export_one_imagetoDrive(img, export_dir, save_file_name, crop_region, res, wait2finished=False)
            tasklist.append(task)
            n += 1

            # save a record in the local dir
            with open(local_record,'w')  as f_obj:
                f_obj.writelines('submitted to GEE on %s'%datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            print('%s: Start %dth task to download'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n))
        except Exception as e:
            error = str(e).split(':')
            if error[0] == 'List.get':
                break
            else:
                raise e

    if wait_all_finished:
        # wait all task filished
        wait_all_task_finished(tasklist)

    return tasklist

def test_gee_download_sentinel1_sar_grd():
    # note, it's weird that when run in "pytest -s ", we cannot import ee,
    # but run directly use python is ok.

    extent_shp = os.path.expanduser('~/Data/flooding_area/Houston/extent/houston_ext1.shp')

    region_name = 'Houston'
    polarization = 'VV' # ['VV','HH','HV','VH']
    mode = 'IW'
    orbit_pass = 'DESCENDING' # ['','ASCENDING']
    product = 'COPERNICUS/S1_GRD'

    start_date= '2017-08-20'
    end_date= '2017-09-01'
    save_dir = './'
    extent_4points = get_extent_xy_from_shp(extent_shp)

    resolution = 10
    projection = 'EPSG:2163' # US National Atlas Equal Area, projection for Houston region.
    global out_projection
    global out_res
    out_projection = projection
    out_res = resolution

    gee_download_sentinel1_sar_grd(region_name, product, polarization, mode, orbit_pass, start_date, end_date,
                                   extent_4points, save_dir,resolution, projection,crop=False, wait_all_finished=True)

    pass



def main(options, args):

    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # for each computer, need to run "earthengine authenticate" first.
    ee.Initialize()

    # all images will save to Google Drive first
    # currently, manually download them from Google Drive
    save_folder = args[0]

    extent_shp = options.extent_shp
    if extent_shp is None:
        raise ValueError('must provide a exent shp')

    region_name = options.region_name
    if region_name is None:
        region_name = os.path.splitext(os.path.basename(extent_shp))[0]

    # checking input shapefile
    # if extent_shp is not None:
    assert io_function.is_file_exist(extent_shp)
    shp_polygon_projection = get_projection_proj4(extent_shp).strip()
    if shp_polygon_projection != '+proj=longlat +datum=WGS84 +no_defs':
        raise ValueError('Only accept %s, please check projection of %s' % (shp_polygon_projection, extent_shp))

    resolution = options.resolution
    global out_res
    out_res = resolution

    projection = options.projection_epsg
    global out_projection
    out_projection = projection

    extent_4points = get_extent_xy_from_shp(extent_shp)

    if os.path.isdir(save_folder) is False:
        io_function.mkdir(save_folder)

    # Sentinel-1 SAR GRD
    # https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD


    # product_list = ['COPERNICUS/S1_GRD']

    # https://developers.google.com/earth-engine/tutorials/community/sar-basics

    # The Level-1 GRD processing assumes an ellipsoid Earth surface model, thus Level-1 GRD needs to be processed to create the Level-2 geocoded,
    # calibrated backscattering coefficients which end up in the GEE COPERNICUS/S1_GRD_FLOAT collection.

    # The S1_GRD_FLOAT collection, and its log-scaled COPERNICUS/S1_GRD computed equivalent,
    # contains "analysis-ready" images, as they have square pixel spacing, with pixel values as FLOAT32,
    # and a UTM projection in the auto-determined zone. Thus, they can be combined with other images, used in feature extractions and reductions, etc.

    product_list = ['COPERNICUS/S1_GRD_FLOAT']
    Polarisation_list = ['VV','HH','HV','VH']
    instrumentMode_list = ['IW']
    orbit_pass_list = ['DESCENDING','ASCENDING']

    b_crop = False  # crop images to input extent

    dates_list = options.dates_list_txt
    if dates_list is None:
        start_date, end_date = options.start_date, options.end_date
        for product in product_list:
            for polarization in Polarisation_list:
                for mode in instrumentMode_list:
                    for orbit_pass in orbit_pass_list:
                        gee_download_sentinel1_sar_grd(region_name, product, polarization, mode, orbit_pass, start_date,
                                                       end_date,
                                                       extent_4points, save_folder, resolution, projection, crop=b_crop,
                                                       wait_all_finished=True)
    else:
        all_task_list = []
        dates_list = io_function.read_list_from_txt(dates_list)
        for idx, date_str in enumerate(dates_list):

            active_tasks = active_task_count(all_task_list)
            while active_tasks > maximum_submit_tasks:
                print('%s : %d (>%d) tasks already in the queue, wait 60 seconds' % (
                datetime.now().strftime('%Y-%m-%d %H:%M:%S'), active_tasks, maximum_submit_tasks))
                time.sleep(60)
                active_tasks = active_task_count(all_task_list)

            print('%d/%d Get snow cover for %s' % (idx, len(dates_list), date_str))
            start_date = date_str
            end_date = timeTools.str2date(start_date, format='%Y-%m-%d') + timedelta(days=1)
            end_date = timeTools.date2str(end_date, format='%Y-%m-%d')

            for product in product_list:
                for polarization in Polarisation_list:
                    for mode in instrumentMode_list:
                        for orbit_pass in orbit_pass_list:
                            tasks = gee_download_sentinel1_sar_grd(region_name, product, polarization, mode, orbit_pass,
                                                           start_date,
                                                           end_date,
                                                           extent_4points, save_folder, resolution, projection,
                                                           crop=b_crop,
                                                           wait_all_finished=False)
                            if tasks is False:
                                continue
                            all_task_list.extend(tasks)

        # wait until all finished
        wait_all_task_finished(all_task_list)



if __name__ == '__main__':

    # ee.Initialize()
    # # environment_test()
    # test_gee_download_sentinel1_sar_grd()
    # sys.exit(0)

    usage = "usage: %prog [options] save_dir "
    parser = OptionParser(usage=usage, version="1.0 2021-05-29")
    parser.description = 'Introduction: get SAR grd files from Google Earth Engine for a study area '
    parser.add_option("-s", "--start_date", default='2016-01-01',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2016-01-01")
    parser.add_option("-e", "--end_date", default='2019-12-31',
                      action="store", dest="end_date",
                      help="the end date for inquiry, with format year-month-day, e.g., 2019-12-31")

    parser.add_option("-d", "--dates_list_txt",
                      action="store", dest="dates_list_txt",
                      help="a list a date, if this is set, the start_date and end_date will be ignored")

    parser.add_option("", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the area extent")

    parser.add_option("-n", "--region_name",
                      action="store", dest="region_name",
                      help="the name of the area, which we download SAR data for")

    parser.add_option('-r',"--resolution", action='store', dest='resolution',default=10, type=float,
                      help='the resolution of output raster')

    parser.add_option('-p',"--projection_epsg", action='store', dest='projection_epsg',
                      help='the projection of output raster')


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


#!/usr/bin/env python
# Filename: snow_cover 
"""
introduction: download snow cover from Google Earth Engine

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 May, 2021
"""

import sys,os
from optparse import OptionParser
from datetime import datetime
from datetime import timedelta

# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools

import vector_gpd

# input earth engine, used: ~/programs/anaconda3/envs/gee/bin/python, or change to ee by "source activate gee"
# need shapely, geopandas, gdal
import ee

from gee_common import export_one_imagetoDrive,wait_all_task_finished

# re-project
from functools import partial
import pyproj
from shapely.ops import transform

def reproject_to_3413_500m(image):
    # reproject to EPSG:3413 (polar scene), with resolution of 500
    reprj_img = image.reproject('EPSG:3413', None, 500)
    return reprj_img

def get_projection_proj4(geo_file):
    import basic_src.map_projection as map_projection
    return map_projection.get_raster_or_vector_srs_info_proj4(geo_file)


def get_snow_product_id(image_info):
    if image_info['type'] != 'Image':
        raise ValueError('The information is not Image type, but is %s'%image_info['type'])
    id_str = image_info['id'].split('/')
    snow_id = '_'.join(id_str[2:])
    return snow_id

# quick test
def environment_test():
    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # ee.Initialize()       # initialize in the previous step
    # for each computer, need to run "earthengine authenticate" first.

    image = ee.Image('srtm90_v4')
    print(image.getInfo())

def gee_download_modis_snow(region_name, snow_product,band_name,start_date, end_date, extent_4points, save_dir,wait_all_finished=True):
    '''
    download snow cover for a date range for a region
    :param snow_products:
    :param band_names:
    :param start_date:
    :param end_date:
    :param extent_4points: [x[1-4], y[1-4]]
    :param save_dir:
    :return:
    '''

    res = 500   #  resolution

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

    snow_cover_collections = ee.ImageCollection(snow_product). \
        filterBounds(polygon_bound). \
        filterDate(start, finish). \
        map(reproject_to_3413_500m)

    # check count  # getInfo can get python number (not ee.Number)
    count = snow_cover_collections.size().getInfo()
    if count < 1:
        basic.outputlogMessage('No %s results for extent: %s polygon '%(snow_product, str(extent_4points)))
        return False
    basic.outputlogMessage('%s snow covers count %d, date %s to %s'%(snow_product,count,start_date,end_date))

    export_dir = 'snow_cover_%s'%region_name
    # download them if not exist
    # Make a list of images
    img_list = snow_cover_collections.toList(snow_cover_collections.size()) # ee.List
    # print(img_list)
    # for idx in range(count):
    #     img = ee.Image(img_list.get(idx)) #.select(img_speci['bands'])
    #     img_info = img.getInfo()
    #     # io_function.save_dict_to_txt_json('img_%d.json'%idx,img_info)
    #     select_image = img.select(band_name)
    #     image_info = select_image.getInfo()
    #     # print(image_info)
    #     print(get_snow_product_id(image_info))

    # export all
    n = 0
    tasklist = []
    while True:
        try:
            img = ee.Image(img_list.get(n)).select(band_name)
            image_info = img.getInfo()
            save_file_name = region_name + '_' + get_snow_product_id(image_info) + '_' + band_name
            crop_region = polygon_bound

            task = export_one_imagetoDrive(img, export_dir, save_file_name, crop_region, res, wait2finished=False)
            tasklist.append(task)
            n += 1
            print('%s: Start %dth task to download snow cover'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n))
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


def get_extent_xy_from_shp(input_shp):
    polygons = vector_gpd.read_polygons_gpd(input_shp)
    if len(polygons) != 1:
        raise ValueError('The extent only support one polygon')
    x,y = vector_gpd.get_polygon_envelope_xy(polygons[0])
    return x,y

def test_get_extent_xy_from_shp():
    # note, it's weird that when run in "pytest -s ", we cannot import ee,
    # but run directly use python is ok.

    print('\ntest_get_extent_xy_from_shp')
    WR_exent = os.path.expanduser('~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent_latlon.shp')
    x,y = get_extent_xy_from_shp(WR_exent)
    print(x)
    print(y)


def main(options, args):
    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # for each computer, need to run "earthengine authenticate" first.
    ee.Initialize()

    # all images will save to Google Drive first, then move them to the below folder
    # not yet implemented, currently, manually download them from Google Drive
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

    extent_4points = get_extent_xy_from_shp(extent_shp)


    if os.path.isdir(save_folder) is False:
        io_function.mkdir(save_folder)

    # snow_products = ['MODIS/006/MYD10A1', 'MODIS/006/MOD10A1']
    snow_products = ['MODIS/006/MYD10A1']   # , 'MODIS/006/MOD10A1'


    band_names = ['NDSI']
    # there are nine bands.
    #       "id": "NDSI_Snow_Cover",
    #       "id": "NDSI_Snow_Cover_Basic_QA",
    #       "id": "NDSI_Snow_Cover_Algorithm_Flags_QA",
    #       "id": "NDSI",
    #       "id": "Snow_Albedo_Daily_Tile",
    #       "id": "orbit_pnt",
    #       "id": "granule_pnt",
    #       "id": "NDSI_Snow_Cover_Class",
    #       "id": "Snow_Albedo_Daily_Tile_Class",

    dates_list = options.dates_list_txt
    if dates_list is None:
        start_date, end_date = options.start_date, options.end_date
        for product in snow_products:
            for band_name in band_names:
                gee_download_modis_snow(region_name, product, band_name, start_date, end_date, extent_4points, save_folder)
    else:
        all_task_list = []
        dates_list = io_function.read_list_from_txt(dates_list)
        for date_str in dates_list:
            print('Get snow cover for %s'%date_str)
            start_date= date_str
            end_date = timeTools.str2date(start_date,format='%Y-%m-%d') + timedelta(days=1)
            end_date = timeTools.date2str(end_date,format='%Y-%m-%d')

            for product in snow_products:
                for band_name in band_names:
                    tasks = gee_download_modis_snow(region_name, product, band_name, start_date, end_date, extent_4points,
                                            save_folder,wait_all_finished=False)
                    all_task_list.extend(tasks)

        # wait until all finished
        wait_all_task_finished(all_task_list)






    pass



if __name__ == "__main__":

    usage = "usage: %prog [options] save_dir "
    parser = OptionParser(usage=usage, version="1.0 2019-12-08")
    parser.description = 'Introduction: get snow cover of a study area '
    parser.add_option("-s", "--start_date",default='2016-01-01',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2016-01-01")
    parser.add_option("-e", "--end_date",default='2019-12-31',
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
                      help="the name of the area, which we download snow cover for")

    # parser.add_option("", "--extent_xy",
    #                   action="store", dest="extent_xy",
    #                   help="a list of xy points")

    # parser.add_option("-b", "--buffer_size",
    #                   action="store", dest="buffer_size", type=int, default = 3000,
    #                   help="the buffer size to crop image in meters")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    basic.setlogfile('get_timelapse_img_gee_%s.log' % str(datetime.date(datetime.now())))

    main(options, args)


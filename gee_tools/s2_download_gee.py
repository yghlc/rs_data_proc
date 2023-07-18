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
from datetime import datetime
from datetime import timedelta

# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.timeTools as timeTools

# input earth engine, used: ~/programs/anaconda3/envs/gee/bin/python, or change to ee by "source activate gee"
# need shapely, geopandas, gdal
import ee

from gee_common import export_one_imagetoDrive,wait_all_task_finished, maximum_submit_tasks, active_task_count, \
    shapely_polygon_to_gee_polygon, reproject

# image specification (more in ChangeDet_DL/dataTools/get_timelapse_img_gee.py)
img_speci = {'sentinel2_rgb_sr': {'product': 'COPERNICUS/S2_SR', 'bands': ['B4', 'B3', 'B2'], 'res': 10}
            }

def gee_download_images(region_name,start_date, end_date, ext_id, extent, product, resolution, projection,
                        band_names, cloud_cover_thr=0.3, crop=False, b_vis=False, wait_all_finished=True):

    start = ee.Date(start_date)  # '%Y-%m-%d'
    finish = ee.Date(end_date)

    polygon_bound = extent
    cloud_cover = 'CLOUDY_PIXEL_PERCENTAGE' #  for Sentinel-2

    filtercollection = ee.ImageCollection(product). \
        filterBounds(polygon_bound). \
        filterDate(start, finish). \
        filter(ee.Filter.lt(cloud_cover, cloud_cover_thr * 100)). \
        sort('CLOUD_COVER', True)

    # # map(cloud_mask). \
    # filter(ee.Filter.calendarRange(month_range[0], month_range[1], 'month')). \

    # check count  # getInfo can get python number (not ee.Number)
    count = filtercollection.size().getInfo()
    if count < 1:
        basic.outputlogMessage('No results for the polygon (id: %d ) with %s' % (ext_id, str(img_speci)))
        return False

    print('Find %d image for the polygon (id: %d )' %(count, ext_id))
    # cannot have a sub folder in Google Drive.
    product_info = product.split('/')
    export_dir = region_name + '_' + product_info[-1] + '_images'
    save_file_name = region_name + '_' + product_info[-1] + '_grid%d'%ext_id

    # save some task record in the local folder
    if os.path.isdir(export_dir) is False:
        io_function.mkdir(export_dir)

    # Make a list of images
    # img_list = filtercollection.toList(filtercollection.size())

    # first_image = ee.Image(filtercollection.first()).select(band_names)
    # # reproject?
    # # default is UTM, reproject to the projection of extent_shp
    # first_image = reproject(first_image,projection,resolution)

    # if omit "toUint16()" , it output float64
    mosaic = filtercollection.select(band_names,band_names).median().toUint16()  # create mosaic using median values
    mosaic = reproject(mosaic, projection, resolution)
    # Error: Image.visualize: Expected a string or list of strings for field 'bands'. (Error code: 3)
    # rgbVis = {'min': 0, 'max': 2000, 'bands': band_names }
    # s2_mosaic_rgb_8bit = mosaic.visualize(rgbVis)
    if b_vis:
        # s2_mosaic_rgb_8bit = mosaic.visualize(bands=band_names, min=0, max=2000)
        mosaic = mosaic.visualize(bands=band_names, min=0, max=2000)
        save_file_name = save_file_name + '_8bit'

    local_record = os.path.join(os.path.join(export_dir, save_file_name+'.submit'))
    if os.path.isfile(local_record):
        print('task %s already be submitted to GEE, skip' % local_record)
        return False

    # only export first image
    # export_one_imagetoDrive(first_image,export_dir,polygon_idx,crop_region, img_speci['res'])
    task = export_one_imagetoDrive(mosaic, export_dir, save_file_name, extent, resolution, wait2finished=wait_all_finished)
    # task = export_one_imagetoDrive(s2_mosaic_rgb_8bit, export_dir, save_file_name+'_8bit', extent, resolution, wait2finished=wait_all_finished)

    # save a record in the local dir
    with open(local_record, 'w') as f_obj:
        f_obj.writelines('submitted to GEE on %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    return task

def gee_download_sentinel2_image(extent_shp, region_name, start_date, end_date,cloud_cover_thr):

    # checking input shapefile
    projection = map_projection.get_raster_or_vector_srs_info_epsg(extent_shp)

    product = img_speci['sentinel2_rgb_sr']['product']
    resolution = img_speci['sentinel2_rgb_sr']['res']
    bands = img_speci['sentinel2_rgb_sr']['bands']

    b_crop = True  # crop images to input extent
    b_visualize = True   # to 8bit, for visualization

    extent_polygons = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:4326')
    extent_ids = vector_gpd.read_attribute_values_list(extent_shp, 'id')
    extent_ids = [int(item) for item in extent_ids]

    all_task_list = []

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

        print('%d/%d Downloading Sentinel-2 for a polygon (id: %d ) ' % (idx+1, len(extent_polygons), ext_id))

        extent_gee = shapely_polygon_to_gee_polygon(extent)

        task = gee_download_images(region_name, start_date, end_date, ext_id, extent_gee, product, resolution,
                                   projection, bands, cloud_cover_thr=cloud_cover_thr,
                                   crop=b_crop, b_vis=b_visualize, wait_all_finished=False)

        if task is False:
            continue

        all_task_list.append(task)

    # wait until all finished
    wait_all_task_finished(all_task_list)

def main(options, args):

    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # for each computer, need to run "earthengine authenticate" first.
    ee.Initialize()

    # all images will save to Google Drive first
    # currently, manually download them from Google Drive

    extent_shp = args[0]
    io_function.is_file_exist(extent_shp)

    region_name = options.region_name
    if region_name is None:
        region_name = os.path.splitext(os.path.basename(extent_shp))[0]
    cloud_cover_thr = options.cloud_cover

    start_date, end_date = options.start_date, options.end_date
    gee_download_sentinel2_image(extent_shp,region_name,start_date, end_date,cloud_cover_thr)




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

    parser.add_option('-r',"--resolution", action='store', dest='resolution',default=10, type=float,
                      help='the resolution of output raster')

    parser.add_option('-p',"--projection_epsg", action='store', dest='projection_epsg',
                      help='the projection of output raster')


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

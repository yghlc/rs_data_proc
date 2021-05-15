#!/usr/bin/env python
# Filename: snow_cover 
"""
introduction: downloaa snow cover from Google Earth Engine

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 May, 2021
"""

import sys,os
from optparse import OptionParser
from datetime import datetime

# path for DeeplabforRS
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import basic_src.basic as basic
import basic_src.io_function as io_function

import vector_gpd

# input earth engine, used: ~/programs/anaconda3/envs/gee/bin/python, or change to ee by "source activate gee"
# need shapely, geopandas, gdal
import ee

# re-project
from functools import partial
import pyproj
from shapely.ops import transform

def reproject_to_3413_500m(image):
    # reproject to EPSG:3413 (polar scene), with resolution of 500
    reprj_img = image.reproject('EPSG:3413', None, 500)
    return reprj_img

# quick test
def environment_test():
    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # ee.Initialize()       # initialize in the previous step
    # for each computer, need to run "earthengine authenticate" first.

    image = ee.Image('srtm90_v4')
    print(image.getInfo())


def main(options, args):
    ee.Initialize()

    # all images will save to Google Drive first, then move them to the below folder.
    save_folder = args[0]

    extent_shp = options.extent_shp




    pass



if __name__ == "__main__":

    usage = "usage: %prog [options] save_dir polygon_shp"
    parser = OptionParser(usage=usage, version="1.0 2019-12-08")
    parser.description = 'Introduction: get snow cover of a study area '
    parser.add_option("-s", "--start_date",default='2016-01-01',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2016-01-01")
    parser.add_option("-e", "--end_date",default='2019-12-31',
                      action="store", dest="end_date",
                      help="the end date for inquiry, with format year-month-day, e.g., 2019-12-31")

    parser.add_option("", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the area extent")

    # parser.add_option("-b", "--buffer_size",
    #                   action="store", dest="buffer_size", type=int, default = 3000,
    #                   help="the buffer size to crop image in meters")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    basic.setlogfile('get_timelapse_img_gee_%s.log' % str(datetime.date(datetime.now())))

    main(options, args)


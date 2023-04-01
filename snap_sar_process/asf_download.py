#!/usr/bin/env python
# Filename: asf_download.py 
"""
introduction: download data from The Alaska Satellite Facility (ASF)

add time: 31 October, 2022

original version: https://github.com/yghlc/Sentinel-1-Pre-Processing/blob/main/asf_download.py
modified by Lingcao Huang, March 2023

"""

import os,sys
from optparse import OptionParser
from datetime import datetime
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

from vector_gpd import read_shape_gpd_to_NewPrj
from vector_gpd import shapefile_to_ROIs_wkt

# install by: conda install -c conda-forge asf_search (or conda activate sar)
import asf_search as asf

from netrc import netrc

from basic_src.io_function import read_list_from_txt

dataset2platform = {'SENTINEL1':asf.PLATFORM.SENTINEL1,
                    'ALOS': asf.PLATFORM.ALOS,
                    'RADARSAT': asf.PLATFORM.RADARSAT
                    } # have more, but need to asf_search first to see

def get_user_password_netrc():
    # Set up authentication using .netrc file
    urs = 'urs.earthdata.nasa.gov'  # Address to call for authentication
    netrcDir = os.path.expanduser("~/.netrc")
    user = netrc(netrcDir).authenticators(urs)[0]
    passwd = netrc(netrcDir).authenticators(urs)[2]
    return user, passwd

def download_data_from_asf_list(file_list_txt, save_dir, username, password):
    ## ROI
    print(datetime.now(),'Searching... ... ...')

    file_ids_list = read_list_from_txt(file_list_txt)
    if len(file_ids_list) < 1:
        raise ValueError('No file list in %s'%file_list_txt)
    results = asf.granule_search(file_ids_list)
    print(datetime.now(),'Found %s results' % (len(results)))
    session = asf.ASFSession()
    session.auth_with_creds(username, password)
    print(datetime.now(),'Downloading... ... ...')

    download_dir = save_dir

    if not os.path.isdir(download_dir):
        os.makedirs(download_dir)

    # it will skip files that have been downloaded
    results.download(path=download_dir, session=session)  # download results to a path
    print(datetime.now(),'Finished Download')

    ## Save results to an output log
    log_filename = os.path.join(download_dir, "download_data.json")
    print(' ')
    print(datetime.now(),'Saving log results to ', log_filename)
    stdoutOrigin = sys.stdout
    # sys.stdout = open (download_dir + region + "_download_log.txt", "w")
    sys.stdout = open(log_filename, "w")
    print(results)
    sys.stdout.close()
    sys.stdout = stdoutOrigin


def download_data_from_asf(idx,roi_count,roi_wkt, save_dir, start_date, end_date, processingLevel, username, password,
                           beamMode='IW',platform=asf.PLATFORM.SENTINEL1,flightDirection='DESCENDING'):
    ## ROI
    print(datetime.now(),'Searching... ... ...')
    if isinstance(platform, list):
        platform_list = platform
    else:
        platform_list = [platform]

    print(datetime.now(), 'Input search parameters:')
    print('roi_wkt',roi_wkt)
    print('save_dir, start_date, end_date', save_dir, start_date, end_date)
    print('processingLevel (filetype), beamMode, platform, flightDirection', processingLevel, beamMode, platform,flightDirection)

    results = asf.geo_search(platform=platform_list, intersectsWith=roi_wkt, start=start_date,
                             end=end_date,
                             beamMode=beamMode, processingLevel=processingLevel,flightDirection=flightDirection)
    print(datetime.now(),'Found %s results' % (len(results)))
    session = asf.ASFSession()
    session.auth_with_creds(username, password)
    print(datetime.now(),'Downloading... ... ...')

    if roi_count == 1:
        download_dir = save_dir
    elif roi_count > 1:
        download_dir = os.path.join(save_dir,'roi_%d'%idx)
    else:
        raise ValueError('There is zero ROI')
    if not os.path.isdir(download_dir):
        os.makedirs(download_dir)

    # it will skip files that have been downloaded
    results.download(path=download_dir, session=session)  # download results to a path
    print(datetime.now(),'Finished Download')

    ## Save results to an output log
    log_filename = os.path.join(download_dir, "download_data.json")
    print(' ')
    print(datetime.now(),'Saving log results to ', log_filename)
    stdoutOrigin = sys.stdout
    # sys.stdout = open (download_dir + region + "_download_log.txt", "w")
    sys.stdout = open(log_filename, "w")
    print(results)
    sys.stdout.close()
    sys.stdout = stdoutOrigin

def main(options, args):
    extent_shp = args[0]
    assert os.path.isfile(extent_shp)

    save_dir = options.save_dir
    start_date = options.start_date
    end_date = options.end_date
    user_name = options.username
    password = options.password
    flightDirection = options.flightdirection.upper()

    if user_name is None or password is None:
        print('Get user name and password from the .netrc file')
        user_name, password = get_user_password_netrc()

    print(datetime.now(), 'download data from ASF, start_date: %s, end_date: %s, user: %s, \nwill save to %s'%(start_date,end_date,user_name,save_dir))

    processingLevel = options.filetype_product.upper()

    platform = dataset2platform[options.dataset_platform.upper()]


    if extent_shp.endswith('.txt'):
        print(datetime.now(), "the input is a TXT file")
        file_list_txt = extent_shp
        download_data_from_asf_list(file_list_txt, save_dir, user_name, password)

    else:
        # shapefile to  ROI
        ROIs_wkt = shapefile_to_ROIs_wkt(extent_shp)
        for idx, roi_wkt in enumerate(ROIs_wkt):
            # download data
            download_data_from_asf(idx, len(ROIs_wkt), roi_wkt, save_dir, start_date, end_date, processingLevel, user_name,
                                   password,
                                   beamMode='IW', platform=platform, flightDirection=flightDirection)



if __name__ == "__main__":

    usage = "usage: %prog [options] extent_shp or file_ids.txt"
    parser = OptionParser(usage=usage, version="1.0 2022-10-31")
    parser.description = 'Introduction: download data from the Alaska Satellite Facility  '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='asf_data',
                      help="the folder to save downloaded data")

    parser.add_option("-s", "--start_date",default='2018-04-30',
                      action="store", dest="start_date",
                      help="start date for inquiry, with format year-month-day, e.g., 2018-05-23")
    parser.add_option("-e", "--end_date",default='2018-06-30',
                      action="store", dest="end_date",
                      help="the end date for inquiry, with format year-month-day, e.g., 2018-05-23")

    parser.add_option("", "--dataset",
                      action="store", dest="dataset_platform", default='SENTINEL1',
                      help="The dataset want to download (Satellite)")

    parser.add_option("", "--flightdirection",
                      action="store", dest="flightdirection", default='DESCENDING',
                      help="The flight direction of SAR imagery, Ascending or Descending")

    parser.add_option("", "--filetype",
                      action="store", dest="filetype_product", default='GRD_HD',
                      help="The data product want to download, such as GRD_HD or SLC")


    parser.add_option("-u", "--username",
                      action="store", dest="username",
                      help="Earth Data account")
    parser.add_option("-p", "--password",
                      action="store", dest="password",
                      help="password for the earth data account")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)


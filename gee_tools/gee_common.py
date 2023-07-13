#!/usr/bin/env python
# Filename: gee_common.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 May, 2021
"""

import os,sys
import time
from datetime import datetime
import ee


# if submit too many task, it will end in errors, like the
# ee.ee_exception.EEException: Too many tasks already in the queue (3000). Please wait for some of them to complete.
# so, set a maximum_submit_tasks
maximum_submit_tasks = 2000

# quick test
def environment_test():
    # https://www.earthdatascience.org/tutorials/intro-google-earth-engine-python-api/
    # ee.Initialize()       # initialize in the previous step
    # for each computer, need to run "earthengine authenticate" first.

    image = ee.Image('srtm90_v4')
    print(image.getInfo())

def reproject(image, new_prj_epsg, resolution):
    # reproject to EPSG:3413 (polar scene), with resolution of 500
    reprj_img = image.reproject(new_prj_epsg, None, resolution)
    return reprj_img

def export_one_imagetoDrive(select_image, save_folder,save_file_name, crop_region, res,maxPixels=1e9, wait2finished=True):


    task = ee.batch.Export.image.toDrive(image=select_image,
                                         region=crop_region,
                                         description=save_file_name,
                                         folder=save_folder,
                                         maxPixels=maxPixels,
                                         scale=res)
    # region=crop_region,

    task.start()
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

def active_task_count(tasks):
    count = 0
    for task in tasks:
        if task.active():
            count += 1
    return count

def wait_all_task_finished(all_tasks):

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

def shapely_polygon_to_gee_polygon(polygon):
    x, y = polygon.exterior.coords.xy
    # # extent
    # x = extent_4points[0]
    # y = extent_4points[1]
    # polygon_bound = ee.Geometry.Polygon([[x[0],y[0]],
    #                                      [x[1], y[1]],
    #                                      [x[2], y[2]],
    #                                      [x[3], y[3]]])
    return ee.Geometry.Polygon([ item for item in zip(x, y)])

# // reproject Polar stere
# function reproject(image) {
#   // print(image);
#   // print(image.get('system:index'));
#   var reprojected = image
#     .reproject('EPSG:3413', null, 500);
#   var dateString = ee.String(image.get('system:index'));//.getInfo().replace(/_/g,"-");
#   // dateString = dateString.replace('_','-');  // bug: only replace the first '_'
#   // print(dateString)
#   // rename to the band description to date of snow cover
#   return reprojected.rename(dateString);//.rename('ndsi'.concat(dateString)); //
# }
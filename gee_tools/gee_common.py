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

from collections import OrderedDict
import numpy as np
import rasterio
from rasterio.transform import Affine
from rasterio.enums import ColorInterp
import pyproj

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
        if task is None:
            continue
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

def directly_save_image_to_local(save_file_path, dtype, image, image_features):
    # for small image < 32M ? we can save it directly to local file

    # outside this function, download data by:
    # patch = image.select(*bands).sampleRectangle(region, defaultValue=0)
    # features = patch.getInfo()  # the actual download
    file_name = save_file_path
    if os.path.isfile(file_name):
        print('%s already exists, skipping downloading'%file_name)
        return

    metadata = image.getInfo()

    coords0 = np.array(image_features["geometry"]["coordinates"][0])    # these are lat/lon?
    # Define the source and target coordinate reference systems
    src_crs = pyproj.CRS.from_string('EPSG:4326')  # Source CRS: WGS84 (latitude and longitude)
    dst_crs_str = metadata['bands'][0]['crs']
    dst_crs = pyproj.CRS.from_string(dst_crs_str)  # Target CRS
    dst_res = metadata['bands'][0]['crs_transform'][0]  # Target resolution
    # Create a Proj transformer to convert coordinates
    transformer = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    # Reproject the point array to EPSG:3413
    lon = coords0[:, 0]
    lat = coords0[:, 1]
    x_new, y_new = transformer.transform(lon, lat)
    # points_news = np.column_stack((x_new, y_new))


    # raster data
    raster = OrderedDict()
    for meta in metadata['bands']:
        band = meta['id']
        img = np.atleast_3d(image_features["properties"][band])
        raster[band] = img.astype(dtype)


    img_all = np.concatenate( [raster[key] for key in raster.keys()] , axis=2)
    height, width, channels = img_all.shape

    transform = (dst_res, 0, min(x_new)-dst_res/2.0, 0, -dst_res, max(y_new) + dst_res/2.0) # # (resX, 0, X_min, 0, -resY, Y_max)
    # transform = (dst_res, 0, x_new[0], 0, -dst_res, y_new[0]) # # (resX, 0, X_min, 0, -resY, Y_max)
    profile = {
        "driver": "GTiff",
        "width": width,
        "height": height,
        "count": channels,
        "crs": dst_crs_str,
        "transform": transform,
        "dtype": img_all.dtype,
        "compress": "lzw",
        "predictor": 2,
        # "tiled": True,  # Set tile size
        # "blockxsize": 256,  # Set block size in X direction
        # "blockysize": 256,  # Set block size in Y direction
    }

    with rasterio.open(file_name, "w", **profile) as f:
        f.write(img_all.transpose(2, 0, 1))
        for i, key in enumerate(raster.keys()):
            f.set_band_description(i+1, key)





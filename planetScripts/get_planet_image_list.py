#!/usr/bin/env python
# Filename: get_planet_image_list 
"""
introduction: get images list of downloaded Planet image in time period and spatial extent
# also create a shapefile to show the extent of download files

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 16 July, 2020
"""

import sys,os
from optparse import OptionParser

# import thest two to make sure load GEOS dll before using shapely
import shapely
from shapely.geometry import mapping # transform to GeJSON format
from shapely.geometry import shape

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import vector_gpd

from datetime import datetime
import json
import time


from xml.dom import minidom


import pandas as pd

def read_cloud_cover(metadata_path):
    xmldoc = minidom.parse(metadata_path)
    nodes = xmldoc.getElementsByTagName("opt:cloudCoverPercentage")
    cloud_per = float(nodes[0].firstChild.data)
    # print(cloud_per)
    return cloud_per

def read_acquired_date(metadata_path):
    xmldoc = minidom.parse(metadata_path)
    nodes = xmldoc.getElementsByTagName("ps:acquisitionDateTime")
    #acquisitionDatetime = float(nodes[0].firstChild.data)
    # 2016-08-23T03:27:21+00:00
    # acquisitionDatetime = datetime.strptime(nodes[0].firstChild.data, '%y-%m-%dT%H:%M:%S')
    acquisitionDate = pd.to_datetime(nodes[0].firstChild.data).date()
    # acquisitionDateTime = pd.to_datetime(nodes[0].firstChild.data).to_pydatetime()
    # print(cloud_per)
    return acquisitionDate

def get_geojson_list_overlap_a_polygon(polygon, geojson_list):
    out_geojson_list = []
    out_shape_list = []
    for geojson_file in geojson_list:
        # print(geojson_file)
        with open(geojson_file) as json_obj:
            geom = json.load(json_obj)
        img_ext = shape(geom)
        inter = polygon.intersection(img_ext)

        if inter.is_empty is False:
            out_shape_list.append(img_ext)
            out_geojson_list.append(geojson_file)
    return out_geojson_list, out_shape_list

def get_Planet_SR_image_list_overlap_a_polygon(polygon,geojson_list, cloud_cover_thr, save_list_path=None):
    '''
    get planet surface reference (SR) list overlap a polygon (within or overlap part of the polygon)
    :param polygon: polygon in the shapely format
    :param geojson_list:
    :param save_list_path: save the list to a txt file.
    :return:
    '''

    image_path_list = []
    cloud_cover_list = []
    for geojson_file in geojson_list:
        # print(geojson_file)
        with open(geojson_file) as json_obj:
            geom = json.load(json_obj)
        img_ext = shape(geom)

        inter = polygon.intersection(img_ext)
        if inter.is_empty is False:
            img_dir = os.path.splitext(geojson_file)[0]
            sr_img_paths = io_function.get_file_list_by_pattern(img_dir,'*_SR.tif')
            if len(sr_img_paths) < 0:
                basic.outputlogMessage('warning, no SR image in %s, try to find Analytic images'%img_dir)
                sr_img_paths = io_function.get_file_list_by_pattern(img_dir, '*AnalyticMS.tif')
                
            meta_data_paths = io_function.get_file_list_by_pattern(img_dir,'*_metadata.xml')
            if len(sr_img_paths) != len(meta_data_paths):
                # raise ValueError('the count of metadata files and images is different')
                basic.outputlogMessage('warning: the count of metadata files and images is different for %s'%img_dir)
                continue

            if len(sr_img_paths) < 1:
                basic.outputlogMessage('warning, no Planet SR image in the %s'%img_dir)
            elif len(sr_img_paths) > 1:
                basic.outputlogMessage('warning, more than one Planet SR image in the %s'%img_dir)
            else:
                # check cloud cover
                cloud_cover = read_cloud_cover(meta_data_paths[0])
                if cloud_cover > cloud_cover_thr:
                    continue

                # add image
                image_path_list.append(sr_img_paths[0])
                cloud_cover_list.append(cloud_cover)

    if save_list_path is not None:
        with open(save_list_path,'w') as f_obj:
            for image_path in image_path_list:
                f_obj.writelines(image_path + '\n')

    return image_path_list, cloud_cover_list

def read_a_meta_of_scene(scene_folder_or_geojson,scene_id_list):

    # get scene id
    if os.path.isfile(scene_folder_or_geojson): # geojson file
       scene_id = os.path.splitext(os.path.basename(scene_folder_or_geojson))[0]
       geojson_path = scene_folder_or_geojson
       scene_folder = os.path.splitext(scene_folder_or_geojson)[0]
    else:
        # scene_folder
        scene_id = os.path.basename(scene_folder_or_geojson)
        geojson_path = scene_folder_or_geojson + '.geojson'
        scene_folder = scene_folder_or_geojson

    # if already exists
    if scene_id in scene_id_list:
        return None,None,None,None,None,None,None,None,None

    print(scene_id)

    # get metadata path
    cloud_cover = 101
    acquisitionDate = datetime(1970,1,1)
    metadata_paths = io_function.get_file_list_by_pattern(scene_folder,'*metadata.xml')
    if len(metadata_paths) < 1:
        basic.outputlogMessage('warning, there is no metadata file in %s'%scene_folder)
    elif len(metadata_paths) > 1:
        basic.outputlogMessage('warning, there are more than one metadata files in %s' % scene_folder)
    else:
        # read metadata
        metadata_path = metadata_paths[0]
        cloud_cover = read_cloud_cover(metadata_path)
        acquisitionDate =  read_acquired_date(metadata_path)

    assets = io_function.get_file_list_by_pattern(scene_folder,'*')
    asset_count = len(assets)
    asset_files = sorted([ os.path.basename(item) for item in assets])
    asset_files =','.join(asset_files)

    image_type = 'analytic'  # 'analytic_sr' (surface reflectance) or 'analytic'
    sr_tif = io_function.get_file_list_by_pattern(scene_folder,'*_SR.tif')
    if len(sr_tif) == 1:
        image_type = 'analytic_sr'

    # consider as the downloading time
    if os.path.isfile(geojson_path):
        modified_time = io_function.get_file_modified_time(geojson_path)
    else:
        geojson_path = ''
        modified_time = io_function.get_file_modified_time(scene_folder)


    return scene_id,cloud_cover,acquisitionDate,geojson_path,scene_folder,asset_count,image_type,asset_files,modified_time


def read_extent_shapefile_epgs4326(exent_shp):
    '''
    read extent polygon, need in latlon projections
    :param exent_shp:
    :return:
    '''

    if exent_shp is None:
        return None

    shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(exent_shp).strip()
    if shp_prj != '+proj=longlat +datum=WGS84 +no_defs':
        raise ValueError('only support the projections of longlat (EPSG:4326)')

    extent_polygons = vector_gpd.read_polygons_gpd(exent_shp)
    return extent_polygons

def get_meta_dict(scene_geojson_or_folder_list):

    scene_id_list = []
    cloud_cover_list = []
    acqui_date_list = []  # acquisitionDate
    geojson_file_list = []
    scene_folder_list = []
    asset_count_list = []
    asset_files_list = []
    image_type_list = []  # 'analytic_sr' (surface reflectance) or 'analytic'
    modife_time_list = []

    scene_without_asset = []  # find scene folders without asset

    for a_scene_file_dir in scene_geojson_or_folder_list:
        # print(id)
        scene_id, cloud_cover, acquisitionDate, geojson_path, scene_folder, asset_count, image_type, asset_files, modified_time = \
            read_a_meta_of_scene(a_scene_file_dir, scene_id_list)

        if scene_id is None:
            continue

        scene_id_list.append(scene_id)
        cloud_cover_list.append(cloud_cover)
        acqui_date_list.append(acquisitionDate)
        geojson_file_list.append(geojson_path)
        scene_folder_list.append(scene_folder)
        asset_count_list.append(asset_count)
        asset_files_list.append(asset_files)
        image_type_list.append(image_type)
        modife_time_list.append(modified_time)

        if asset_count == 0:
            scene_without_asset.append(a_scene_file_dir)

    meta_table = {'scene_id': scene_id_list,
                   'cloud_cover': cloud_cover_list,
                   'acquisitionDate': acqui_date_list,
                   'downloadTime': modife_time_list,
                   'asset_count': asset_count_list,
                   'image_type': image_type_list,
                   'asset_files': asset_files_list,
                   'geojson': geojson_file_list,
                   'folder': scene_folder_list
                   }

    # # put then in order
    # import collections
    # dict_top1_per = collections.OrderedDict(sorted(dict_top1_per.items()))


    return meta_table, scene_without_asset


def save_planet_images_to_excel(image_dir,save_xlsx):

    if os.path.isfile(image_dir):
        basic.outputlogMessage('Warning, Input %s is a file, expected a folder, skip it'%image_dir)
        return False

    # read save_xlsx if it exist

    # may be not good to exclude these the scene id if we want to update some records.
    # if os.path.isfile(save_xlsx):
    #     df = pd.read_excel(save_xlsx)
    #     scene_id_list.extend(df['scene_id'].to_list())



    # get scene folders (idealy, the number of scene folder are the same to the one of geojson files)
    scene_geojson_folders = io_function.get_file_list_by_pattern(image_dir,'????????_??????_*')     # acquired date_time
    if len(scene_geojson_folders) < 1:
        raise ValueError('There is no scene folder or geojson in %s'%image_dir)

    scene_table, scene_without_asset = get_meta_dict(scene_geojson_folders)
    df = pd.DataFrame(scene_table)

    with pd.ExcelWriter(save_xlsx) as writer:
        df.to_excel(writer)
        basic.outputlogMessage('write records of downloaded scenes to %s'%save_xlsx)

    scene_folder_no_assets_txt = os.path.splitext(save_xlsx)[0] + '_scenes_noAsset.txt'
    with open(scene_folder_no_assets_txt, 'w') as f_obj:
        for scene_dir in scene_without_asset:
            f_obj.writelines(scene_dir + '\n')

    return True

def group_geojson_by_date(geojson_list):
    geojson_group = {}
    for item in geojson_list:
        basename = os.path.basename(item)
        # 20200817_192455_82_1061.geojson
        date = timeTools.str2date(basename[:8],format='%Y%m%d')
        if date in geojson_group.keys():
            geojson_group[date].append(item)
        else:
            geojson_group[date] = [item]

    return geojson_group

def save_planet_images_to_shapefile(geojson_list, save_shp_path, wkt_string, extent_polygon=None, b_group_date=False):
    '''
    get the meta data and extent of download
    :param geojson_list: geojson_list
    :param save_shp_path:
    :param extent_polygon: a extent polygon
    :param b_group_date:
    :return:
    '''

    # remove incomplete scenes
    geojson_list = [ item for item in geojson_list if 'incomplete_scenes' not in item ]
    if len(geojson_list) < 1:
        raise ValueError('No geojson files (exclude incomplete_scenes) the given folder')

    if extent_polygon is not None and len(extent_polygon) > 1:
        raise ValueError('Only support one extent polygon')

    extent = extent_polygon[0]

    if b_group_date is False:
        geojson_group = {'all': geojson_list}
    else:
        geojson_group = group_geojson_by_date(geojson_list)

    for key in geojson_group.keys():
        sub_geojsons = geojson_group[key]
        if len(sub_geojsons) < 1:
            continue

        sel_geojson_list, sel_polygons = get_geojson_list_overlap_a_polygon(extent,sub_geojsons)
        scene_table, scene_without_asset = get_meta_dict(sel_geojson_list)


        if len(scene_table['scene_id']) != len(sel_polygons):
            raise ValueError('The count of scence ID and polygon are different ')

        # to strings, ESRI Shapefile does not support datetime fields
        scene_table['acquisitionDate'] = [ timeTools.datetime2str(item) for item in  scene_table['acquisitionDate']]
        scene_table['downloadTime'] = [ timeTools.datetime2str(item) for item in  scene_table['downloadTime']]

        scene_table['Polygons'] = sel_polygons
        df = pd.DataFrame(scene_table)

        if key=="all":
            save_path = save_shp_path
        else:
            date_str = timeTools.date2str(key)
            save_path = io_function.get_name_by_adding_tail(save_shp_path,date_str)
        vector_gpd.save_polygons_to_files(df,'Polygons', wkt_string, save_path)

    return True

def main(options, args):

    image_dir = args[0]

    # get the file list in folder and save to excel
    if options.save_xlsx_path is not None:
        save_xlsx = options.save_xlsx_path
        save_planet_images_to_excel(image_dir,save_xlsx)
        return True

    cloud_cover_thr = options.cloud_cover  # 0.3
    cloud_cover_thr =  cloud_cover_thr* 100     # in xml, it is percentage

    # try to get all geojson under image_dir
    geojson_list = io_function.get_file_list_by_ext('.geojson',image_dir,bsub_folder=True)
    if len(geojson_list) < 1:
        raise ValueError('There is no geojson files in %s'%image_dir)

    extent_polygons = read_extent_shapefile_epgs4326(options.extent_shp)


    if options.save_meta_shapefile is not None:
        wkt_string = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
        save_planet_images_to_shapefile(geojson_list, options.save_meta_shapefile, wkt_string,
                                        extent_polygon=extent_polygons, b_group_date=options.group_date)
        return True

    for idx, polygon in enumerate(extent_polygons):
        save_list_txt = 'planet_sr_image_poly_%d.txt'%idx
        get_Planet_SR_image_list_overlap_a_polygon(polygon,geojson_list,cloud_cover_thr,save_list_path=save_list_txt)

        pass




    #TODO: filter the image using the acquired time


    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] image_dir "
    parser = OptionParser(usage=usage, version="1.0 2020-07-16")
    parser.description = 'Introduction: get planet image list  '
    # parser.add_option("-s", "--start_date",default='2018-04-30',
    #                   action="store", dest="start_date",
    #                   help="start date for inquiry, with format year-month-day, e.g., 2018-05-23")
    # parser.add_option("-e", "--end_date",default='2018-06-30',
    #                   action="store", dest="end_date",
    #                   help="the end date for inquiry, with format year-month-day, e.g., 2018-05-23")
    parser.add_option("-e", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the path for extent, shapefile (EPSG:4326)")

    parser.add_option("-c", "--cloud_cover",
                      action="store", dest="cloud_cover", type=float,default=0.3,
                      help="the could cover threshold, only accept images with cloud cover less than the threshold")
    parser.add_option("-g", "--group_date",
                      action="store_true", dest="group_date",default=False,
                      help="true to group image if their acquisition date is the same")
    parser.add_option("-x", "--save_xlsx_path",
                      action="store", dest="save_xlsx_path",
                      help="save the sence lists to xlsx file, if this is set, save_meta_shapefile will be ignored")
    parser.add_option("-m", "--save_meta_shapefile",
                      action="store", dest="save_meta_shapefile",
                      help="the path for saving meta and extent of downloaded planet images")



    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    basic.setlogfile('get_planet_image_list_%s.log'%str(datetime.date(datetime.now())))

    main(options, args)


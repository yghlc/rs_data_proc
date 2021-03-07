#!/usr/bin/env python
# Filename: dem_mosaic_crop.py
"""
introduction: mosaic of arctic DEM and crop them to a given extent

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 27 February, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.RSImageProcess as RSImageProcess
import basic_src.timeTools as timeTools
import raster_io

import re
re_stripID='[0-9]{8}_[0-9A-F]{16}_[0-9A-F]{16}'

import multiprocessing
from multiprocessing import Pool
import operator


def subset_image_by_polygon_box(in_img, out_img, polygon,resample_m='bilinear',o_format='GTiff', out_res=None,same_extent=False, thread_num=1):
    if same_extent:
        return RSImageProcess.subset_image_by_polygon_box(out_img,in_img,polygon,resample_m=resample_m, o_format=o_format,
                                                          xres=out_res,yres=out_res,compress='lzw', tiled='yes', bigtiff='if_safer',thread_num=thread_num)
    else:
        # crop to the min extent (polygon or the image)
        return RSImageProcess.subset_image_by_polygon_box_image_min(out_img,in_img,polygon,resample_m=resample_m,o_format=o_format,
                                                            xres=out_res,yres=out_res,compress='lzw', tiled='yes', bigtiff='if_safer',thread_num=thread_num)


def group_demTif_yearmonthDay(demTif_list, diff_days=30):
    '''
    groups DEM tif if the acquisition of their raw images is close (less than 30 days or others)
    :param demTif_list:
    :return:
    '''
    dem_groups = {}
    for tif in demTif_list:
        yeardate =  timeTools.get_yeardate_yyyymmdd(os.path.basename(tif))

        b_assgined = False
        for time in dem_groups.keys():
            if timeTools.diff_yeardate(time,yeardate) <= diff_days:
                dem_groups[time].append(tif)
                b_assgined = True
                break
        if b_assgined is False:
            dem_groups[yeardate] = [tif]

    return dem_groups

def group_demTif_strip_pair_ID(demTif_list):
    '''
    group dem tif based on the same pair ID, such as 20170226_1030010066648800_1030010066CDE700 (did not include satellite name)
    :param demTif_list:
    :return:
    '''

    dem_groups = {}
    for tif in demTif_list:
        strip_ids = re.findall(re_stripID, os.path.basename(tif))
        if len(strip_ids) != 1:
            print(strip_ids)
            raise ValueError('found zero or multiple strip IDs in %s, expect one'%tif)
        strip_id = strip_ids[0]
        if strip_id in dem_groups.keys():
            dem_groups[strip_id].append(tif)
        else:
            dem_groups[strip_id] = [tif]

    return dem_groups

def mosaic_dem_list(key, dem_list, save_tif_dir,resample_method,save_source, o_format, thread_num=1):

    save_mosaic = os.path.join(save_tif_dir, key + '.tif')
    # check file existence
    # if os.path.isfile(save_mosaic):
    b_save_mosaic = io_function.is_file_exist_subfolder(save_tif_dir, key + '.tif')
    if b_save_mosaic is not False:
        basic.outputlogMessage('warning, mosaic file: %s exist, skip' % b_save_mosaic)
        return b_save_mosaic
        # mosaic_list.append(b_save_mosaic)
        # continue
    # save the source file for producing the mosaic
    if save_source:
        save_mosaic_source_txt = os.path.join(save_tif_dir, key + '_src.txt')
        io_function.save_list_to_txt(save_mosaic_source_txt, dem_list)

    # if only one dem, then copy it if it's not VRT format
    if len(dem_list) == 1:
        if raster_io.get_driver_format(dem_list[0]) != 'VRT':
            io_function.copy_file_to_dst(dem_list[0], save_mosaic)
            return save_mosaic

    # create mosaic, can handle only input one file, but is slow
    result = RSImageProcess.mosaic_crop_images_gdalwarp(dem_list, save_mosaic, resampling_method=resample_method,
                                               o_format=o_format,
                                               compress='lzw', tiled='yes', bigtiff='if_safer',thread_num=thread_num)
    if result is False:
        sys.exit(1)
        # return False
    return save_mosaic

def mosaic_dem_same_stripID(demTif_groups,save_tif_dir, resample_method, process_num=1, save_source=False, o_format='GTiff'):
    if os.path.isdir(save_tif_dir):
        io_function.mkdir(save_tif_dir)

    mosaic_list = []
    if process_num == 1:
        for key in demTif_groups.keys():
            save_mosaic = mosaic_dem_list(key, demTif_groups[key], save_tif_dir,resample_method,save_source, o_format,thread_num=process_num)
            mosaic_list.append(save_mosaic)
    elif process_num > 1:
        # change backk to multi process of gdalwarp, when get mosaic, gdalwarp multi-thread cannot fully utlized CPUs
        theadPool = Pool(process_num)  # multi processes

        parameters_list = [(key, demTif_groups[key], save_tif_dir, resample_method, save_source, o_format,1) for key in demTif_groups.keys()]

        results = theadPool.starmap(mosaic_dem_list, parameters_list)  # need python3
        mosaic_list = [ out for out in results if out is not False]
    else:
        raise ValueError('Wrong process_num: %d'%process_num)

    return mosaic_list

def mosaic_dem_date(demTif_date_groups,save_tif_dir, resample_method,process_num=1,save_source=False,o_format='GTiff'):

    # convert the key in demTif_date_groups to string
    date_groups = {}
    for key in demTif_date_groups.keys():
        new_key = key.strftime('%Y%m%d') + '_dem'
        date_groups[new_key] = demTif_date_groups[key]

    # becuase the tifs have been grouped, so we can use mosaic_dem_same_stripID
    return mosaic_dem_same_stripID(date_groups,save_tif_dir,resample_method, process_num=process_num,save_source=save_source,o_format=o_format)

def check_dem_valid_per(dem_tif_list, work_dir, process_num =1, move_dem_threshold = None, area_pixel_num=None):
    '''
    get the valid pixel percentage for each DEM
    :param dem_tif_list:
    :param work_dir:
    :param move_dem_threshold: move a DEM to a sub-folder if its valid percentage small then the threshold
    :return:
    '''

    keep_dem_list = []
    print('start getting valid pixel percent for %d files'%len(dem_tif_list))
    dem_tif_valid_per = {}
    if process_num == 1:
        for idx,tif in enumerate(dem_tif_list):
            # RSImage.get_valid_pixel_count(tif)
            # per = RSImage.get_valid_pixel_percentage(tif,total_pixel_num=area_pixel_num)
            print('(%d/%d) get valid pixel percent for %s'%(idx+1, len(dem_tif_list),tif))
            per = raster_io.get_valid_pixel_percentage(tif, total_pixel_num=area_pixel_num)
            if per is False:
                return False
            dem_tif_valid_per[tif] = per
            keep_dem_list.append(tif)
    elif process_num > 1:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(tif, area_pixel_num, '%d/%d'%(idx+1, len(dem_tif_list)) ) for idx,tif in enumerate(dem_tif_list)]
        results = theadPool.starmap(raster_io.get_valid_pixel_percentage, parameters_list)  # need python3
        for res, tif in zip(results, dem_tif_list):
            if res is False:
                return False
            dem_tif_valid_per[tif] = res
            keep_dem_list.append(tif)
    else:
        raise ValueError("Wrong process_num: %d"%process_num)
    # sort
    dem_tif_valid_per_d = dict(sorted(dem_tif_valid_per.items(), key=operator.itemgetter(1), reverse=True))
    percent_txt = os.path.join(work_dir,'dem_valid_percent.txt')
    with open(percent_txt,'w') as f_obj:
        for key in dem_tif_valid_per_d:
            f_obj.writelines('%s %.4f\n'%(os.path.basename(key),dem_tif_valid_per_d[key]))
        basic.outputlogMessage('save dem valid pixel percentage to %s'%percent_txt)

    # only keep dem with valid pixel greater than a threshold
    if move_dem_threshold is not None:  # int or float
        keep_dem_list = []      # reset the list
        mosaic_dir_rm = os.path.join(work_dir,'dem_valid_lt_%.2f'%move_dem_threshold)
        io_function.mkdir(mosaic_dir_rm)
        for tif in dem_tif_valid_per.keys():
            if dem_tif_valid_per[tif] < move_dem_threshold:
                io_function.movefiletodir(tif,mosaic_dir_rm)
            else:
                keep_dem_list.append(tif)

    return keep_dem_list

def mosaic_crop_dem(dem_tif_list, save_dir, extent_id, extent_poly, b_mosaic_id, b_mosaic_date, process_num,
                         keep_dem_percent, o_res, pre_name, resample_method='average'):

    # crop to the same extent
    crop_tif_dir = os.path.join(save_dir, 'dem_crop_sub_%d' % extent_id)
    if os.path.isdir(crop_tif_dir) is False:
        io_function.mkdir(crop_tif_dir)
    crop_tif_list = []
    for tif in dem_tif_list:
        save_crop_path = os.path.join(crop_tif_dir, os.path.basename(io_function.get_name_by_adding_tail(tif, 'sub_poly_%d' % extent_id)) )
        if os.path.isfile(save_crop_path):
            basic.outputlogMessage('%s exists, skip cropping' % save_crop_path)
            crop_tif_list.append(save_crop_path)
        else:
            crop_tif = subset_image_by_polygon_box(tif, save_crop_path, extent_poly, resample_m='near',
                            o_format='VRT', out_res=o_res,same_extent=True,thread_num=process_num)
            if crop_tif is False:
                raise ValueError('warning, crop %s failed' % tif)
            crop_tif_list.append(crop_tif)
    dem_tif_list = crop_tif_list

    # area pixel count
    area_pixel_count = int(extent_poly.area / (o_res*o_res))
    basic.outputlogMessage('Area pixel count: %d'%area_pixel_count)

    dem_groups = group_demTif_strip_pair_ID(dem_tif_list)

    # create mosaic (dem with the same strip pair ID)
    mosaic_dir = os.path.join(save_dir, 'dem_stripID_mosaic_sub_%d' % extent_id)
    if b_mosaic_id:
        if os.path.isfile(os.path.join(mosaic_dir, 'dem_valid_percent.txt')):
            basic.outputlogMessage('mosaic based on stripID exists, skip mosaicking')
            with open(os.path.join(mosaic_dir, 'dem_valid_percent.txt')) as f_job:
                tif_valid_per_list = [line.strip().split() for line in f_job.readlines()]
                # check keep_dem_percent
                dem_tif_list = [ os.path.join(mosaic_dir,tif) for tif, per in tif_valid_per_list 
                if float(per) >= keep_dem_percent]

        else:
            io_function.mkdir(mosaic_dir)
            # when create mosaic using VRT end some wrong results, so choose to use 'GTiff'
            # for creating a mosaic with VRT format, we should use "gdalbuildvrt"
            mosaic_list = mosaic_dem_same_stripID(dem_groups, mosaic_dir, resample_method, process_num=process_num,
                                                  o_format='GTiff')
            dem_tif_list = mosaic_list

            # get valid pixel percentage (due to large memory used in raster_io, may set process_num as 1)
            dem_tif_list = check_dem_valid_per(dem_tif_list, mosaic_dir, process_num=process_num,
                                               move_dem_threshold=keep_dem_percent, area_pixel_num=area_pixel_count)


    # groups DEM with original images acquired at the same year months
    dem_groups_date = group_demTif_yearmonthDay(dem_tif_list,diff_days=31)
    # sort based on yeardate in accending order : operator.itemgetter(0)
    dem_groups_date = dict(sorted(dem_groups_date.items(), key=operator.itemgetter(0)))
    # save to txt (json format)
    year_date_txt = os.path.join(mosaic_dir, 'year_date_tif.txt')
    io_function.save_dict_to_txt_json(year_date_txt,dem_groups_date)

    # merge DEM with close acquisition date
    mosaic_yeardate_dir = os.path.join(save_dir,'dem_date_mosaic_sub_%d'%extent_id)
    if b_mosaic_date:
        if os.path.isfile(os.path.join(mosaic_yeardate_dir,'dem_valid_percent.txt')):
            basic.outputlogMessage('mosaic based on acquisition date exists, skip mosaicking')
            with open(os.path.join(mosaic_yeardate_dir,'dem_valid_percent.txt')) as f_job:
                tif_valid_per_list = [line.strip().split() for line in f_job.readlines()]
                # check keep_dem_percent
                dem_tif_list = [ os.path.join(mosaic_yeardate_dir,tif) for tif, per in tif_valid_per_list
                if float(per) >= keep_dem_percent]
        else:
            io_function.mkdir(mosaic_yeardate_dir)
            # this is the output of mosaic, save to 'GTiff' format.
            mosaic_list = mosaic_dem_date(dem_groups_date,mosaic_yeardate_dir,resample_method, process_num=process_num, save_source=True, o_format='GTiff')
            dem_tif_list = mosaic_list

            # get valid pixel percentage
            dem_tif_list = check_dem_valid_per(dem_tif_list,mosaic_yeardate_dir,process_num=process_num, move_dem_threshold = keep_dem_percent,area_pixel_num=area_pixel_count)


    return dem_tif_list

def get_dem_tif_ext_polygons(dem_list):
    bounding_box_list = [ raster_io.get_image_bound_box(tif) for tif in dem_list ]
    ext_polys = [ vector_gpd.convert_image_bound_to_shapely_polygon(box) for box in  bounding_box_list]
    return ext_polys


def main(options, args):

    save_dir = options.save_dir
    extent_shp = options.extent_shp
    process_num = options.process_num
    o_res = options.out_res
    b_mosaic_id = options.create_mosaic_id
    b_mosaic_date = options.create_mosaic_date
    keep_dem_percent = options.keep_dem_percent

    dem_dir_or_txt = args[0]
    if os.path.isfile(dem_dir_or_txt):
        dem_list = io_function.read_list_from_txt(dem_dir_or_txt)
    else:
        dem_list = io_function.get_file_list_by_ext('.tif', dem_dir_or_txt, bsub_folder=False)
        dem_list = [ tif for tif in dem_list if 'matchtag' not in tif ] # remove matchtag
    dem_count = len(dem_list)
    if dem_count < 1:
        raise ValueError('No input dem files in %s' % dem_dir_or_txt)

    resample_method= 'average'


    if extent_shp is None:
        # groups DEM based on the same strip ID
        dem_groups = group_demTif_strip_pair_ID(dem_list)
        # mosaic them direclty without consider the extent
        mosaic_dir = os.path.join(save_dir, 'dem_stripID_mosaic' )
        mosaic_dem_same_stripID(dem_groups, mosaic_dir, resample_method, process_num=process_num, save_source=True,
                                o_format='GTiff')
    else:
        extent_shp_base = os.path.splitext(os.path.basename(extent_shp))[0]
        dem_prj = map_projection.get_raster_or_vector_srs_info_epsg(dem_list[0])
        extent_prj = map_projection.get_raster_or_vector_srs_info_epsg(extent_shp)

        # # check projection (time-consuming if there are many tif files)
        # for dem_tif in dem_list:
        #     prj = map_projection.get_raster_or_vector_srs_info_epsg(dem_tif)
        #     if dem_prj != prj:
        #         raise ValueError('The projection inconsistent among dems (%s is different)'%dem_tif)

        dem_ext_polys = get_dem_tif_ext_polygons(dem_list)

        if extent_prj==dem_prj:
            extent_polys = vector_gpd.read_polygons_gpd(extent_shp)
        else:
            extent_polys = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp,dem_prj)

        if len(extent_polys) < 1:
            raise ValueError('No polygons in %s' % extent_shp)
        else:
            basic.outputlogMessage('%d extent polygons in %s' % (len(extent_polys), extent_shp))

        extPolys_ids = vector_gpd.read_attribute_values_list(extent_shp, 'id')
        if extPolys_ids is None or None in extPolys_ids:
            basic.outputlogMessage('Warning, field: id is not in %s, will create default ID for each grid' % extent_shp)
            extPolys_ids = [id + 1 for id in range(len(extent_polys))]

        for idx, ext_poly in zip(extPolys_ids, extent_polys):
            basic.outputlogMessage('mosaic and crop DEM for the %d th extent (%d in total)' % (idx, len(extent_polys)))
            # get subset of DEM
            dem_poly_ids = vector_gpd.get_poly_index_within_extent(dem_ext_polys, ext_poly)
            if len(dem_poly_ids) < 1:
                basic.outputlogMessage('no dem tifs within %d polygons'%idx)
                continue
            dem_list_sub = [dem_list[id] for id in dem_poly_ids]

            mosaic_crop_dem(dem_list_sub, save_dir, idx, ext_poly, b_mosaic_id, b_mosaic_date,
                                 process_num, keep_dem_percent, o_res, extent_shp_base, resample_method='average')



if __name__ == '__main__':
    usage = "usage: %prog [options] dem_tif_dir or dem_list_txt "
    parser = OptionParser(usage=usage, version="1.0 2020-12-27")
    parser.description = 'Introduction: mosaic and crop for dem files '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")

    parser.add_option("-e", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the extent file for cropping")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    parser.add_option("-p", "--keep_dem_percent",
                      action="store", dest="keep_dem_percent",type=float,default=30.0,
                      help="keep dem with valid percentage greater than this value")

    parser.add_option("-o", "--out_res",
                      action="store", dest="out_res",type=float,default=2.0,
                      help="the resolution for final output")

    parser.add_option("-m", "--create_mosaic_id",
                      action="store_true", dest="create_mosaic_id",default=False,
                      help="for a small region, if true, then get a mosaic of dem with the same ID (date_catalogID_catalogID)")

    parser.add_option("-t", "--create_mosaic_date",
                      action="store_true", dest="create_mosaic_date",default=False,
                      help="for a small region, if true, then get a mosaic of dem with close acquisition date ")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
#!/usr/bin/env python
# Filename: dem_difference.py 
"""
introduction: conduct DEM difference

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 27 February, 2021
"""
import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.timeTools as timeTools
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function
import split_image

import numpy as np
from itertools import combinations
import operator


import multiprocessing
from multiprocessing import Pool

def check_dem_valid_per(dem_tif_list, work_dir, process_num =1, move_dem_threshold = None, area_pixel_num=None):
    '''
    get the valid pixel percentage for each DEM
    :param dem_tif_list:
    :param work_dir:
    :param move_dem_threshold: move a DEM to a sub-folder if its valid percentage small then the threshold
    :return:
    '''

    keep_dem_list = []

    dem_tif_valid_per = {}
    if process_num == 1:
        for tif in dem_tif_list:
            # RSImage.get_valid_pixel_count(tif)
            # per = RSImage.get_valid_pixel_percentage(tif,total_pixel_num=area_pixel_num)
            per = raster_io.get_valid_pixel_percentage(tif, total_pixel_num=area_pixel_num)
            if per is False:
                return False
            dem_tif_valid_per[tif] = per
            keep_dem_list.append(tif)
    elif process_num > 1:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(tif, area_pixel_num) for tif in dem_tif_list]
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


def read_date_dem_to_memory(pair_idx, pair, date_pair_list_sorted,dem_data_dict, dem_groups_date, less_memory=False,boundary=None):

    if less_memory is False:
        # read data to memory if need, then store in memory, avoid to read them again.
        # for a large area, because we read all raster to memory, it will cause "out of memory problem"
        if pair[0] not in dem_data_dict.keys():
            data_old, nodata_old = raster_io.read_raster_one_band_np(dem_groups_date[pair[0]][0],boundary=boundary)
            data_old[data_old == nodata_old] = np.nan
            dem_data_dict[pair[0]] = data_old
        else:
            data_old = dem_data_dict[pair[0]]

        # read data to memory if need
        if pair[1] not in dem_data_dict.keys():
            data_new, nodata_new = raster_io.read_raster_one_band_np(dem_groups_date[pair[1]][0],boundary=boundary)
            data_new[data_new == nodata_new] = np.nan
            dem_data_dict[pair[1]] = data_new
        else:
            data_new = dem_data_dict[pair[1]]
    else:
        # if we don't have enough memory, don't store the all DEM data in memory, only read two needed.
        # wil increase reading operation from disk
        data_old, nodata_old = raster_io.read_raster_one_band_np(dem_groups_date[pair[0]][0],boundary=boundary)
        data_new, nodata_new = raster_io.read_raster_one_band_np(dem_groups_date[pair[1]][0],boundary=boundary)

        # replace nodata with nan
        data_old[ data_old == nodata_old ] = np.nan
        data_new[ data_new == nodata_new ] = np.nan

    # release some memory if we can (NO)

    return data_old, data_new


def dem_diff_newest_oldest(dem_tif_list, out_dem_diff, out_date_diff):
    '''
    get DEM difference, for each pixel, newest vaild value - oldest valid value
    :param dem_list:
    :param output:
    :return:
    '''
    if len(dem_tif_list) < 2:
        basic.outputlogMessage('error, the count of DEM is smaller than 2')
        return False


    # groups DEM with original images acquired at the same year months
    dem_groups_date = group_demTif_yearmonthDay(dem_tif_list,diff_days=0)
    # sort based on yeardate in accending order : operator.itemgetter(0)
    dem_groups_date = dict(sorted(dem_groups_date.items(), key=operator.itemgetter(0)))
    txt_save_path = os.path.splitext(out_date_diff)[0]+'.txt'
    io_function.save_dict_to_txt_json(txt_save_path,dem_groups_date)

    date_list = list(dem_groups_date.keys())
    dem_tif_list = [ dem_groups_date[key][0] for key in dem_groups_date.keys()]  # each date, only have one tif
    tif_obj_list = [ raster_io.open_raster_read(tif) for tif in dem_tif_list]


    height, width, _ = raster_io.get_width_heigth_bandnum(tif_obj_list[0])

    # check them have the width and height
    for tif, obj in zip(dem_tif_list[1:],tif_obj_list[1:]):
        h, w, _ = raster_io.get_width_heigth_bandnum(obj)
        if h!=height or w!=width:
            raise ValueError('the height and width of %s is different from others'%tif)

    # divide the image the many small patches, then calcuate one by one, solving memory issues.
    image_patches = split_image.sliding_window(width,height, 1024, 1024,adj_overlay_x=0,adj_overlay_y=0)
    patch_count = len(image_patches)
    tif_obj_list = None

    # read all and their date
    date_pair_list = list(combinations(date_list, 2))
    date_diff_list = [ (item[1] - item[0]).days for item in date_pair_list ]
    # sort based on day difference (from max to min)
    date_pair_list_sorted = [x for _, x in sorted(zip(date_diff_list, date_pair_list),reverse=True)]    # descending


    # get the difference
    date_diff_np = np.zeros((height, width),dtype=np.uint16)
    dem_diff_np = np.empty((height, width),dtype=np.float32)
    dem_diff_np[:] = np.nan

    for idx, patch in enumerate(image_patches):
        print('tile: %d / %d'%(idx+1, patch_count))

        patch_w = patch[2]
        patch_h = patch[3]
        patch_date_diff = np.zeros((patch_h, patch_w),dtype=np.uint16)
        patch_dem_diff = np.empty((patch_h, patch_w),dtype=np.float32)
        patch_dem_diff[:] = np.nan

        # use dict to read data from disk (only need)
        dem_data_dict = {}
        for p_idx, pair in enumerate(date_pair_list_sorted):
            diff_days = (pair[1] - pair[0]).days
            basic.outputlogMessage('Getting DEM difference using the one on %s and %s, total day diff: %d'%
                                (timeTools.date2str(pair[1]), timeTools.date2str(pair[0]),diff_days))
            # print(pair,':',(pair[1] - pair[0]).days)

            data_old, data_new = read_date_dem_to_memory(p_idx, pair, date_pair_list_sorted,dem_data_dict, dem_groups_date,boundary=patch)

            # print('data_old shape:',data_old.shape)
            # print('data_new shape:',data_new.shape)

            diff_two = data_new - data_old
            # print(diff_two)

            # fill the element
            new_ele = np.where(np.logical_and( np.isnan(patch_dem_diff), ~np.isnan(diff_two)))

            patch_dem_diff[new_ele] = diff_two[new_ele]
            patch_date_diff[new_ele] = diff_days

            # check if all have been filled ( nan pixels)
            diff_remain_hole = np.where(np.isnan(patch_dem_diff))
            # basic.outputlogMessage(' remain %.4f percent pixels need to be filled'% (100.0*diff_remain_hole[0].size/patch_dem_diff.size) )
            if diff_remain_hole[0].size < 1:
                break

        # copy to the entire image
        row_s = patch[1]
        row_e = patch[1] + patch[3]
        col_s = patch[0]
        col_e = patch[0] + patch[2]
        dem_diff_np[row_s:row_e, col_s:col_e] = patch_dem_diff
        date_diff_np[row_s:row_e, col_s:col_e] = patch_date_diff

    # save date diff to tif (16 bit)
    raster_io.save_numpy_array_to_rasterfile(date_diff_np,out_date_diff,dem_tif_list[0], nodata=0,compress='lzw',tiled='yes',bigtiff='if_safer')

    # # stretch the DEM difference, save to 8 bit.
    # dem_diff_np_8bit = raster_io.image_numpy_to_8bit(dem_diff_np,10,-10,dst_nodata=0)
    # out_dem_diff_8bit = io_function.get_name_by_adding_tail(out_dem_diff, '8bit')
    # raster_io.save_numpy_array_to_rasterfile(dem_diff_np_8bit, out_dem_diff_8bit, dem_tif_list[0], nodata=0)


    # if possible, save to 16 bit, to save the disk storage.
    # dem_diff_np[0:5,0] = -500
    # dem_diff_np[0,0:5] = 500
    # print(np.nanmin(dem_diff_np))
    # print(np.nanmax(dem_diff_np))
    range = np.iinfo(np.int16)
    # dem_diff_np_cm = dem_diff_np*100
    # if np.nanmin(dem_diff_np_cm) < range.min or np.nanmax(dem_diff_np_cm) > range.max:
    # save dem diff to files (float), meter
    raster_io.save_numpy_array_to_rasterfile(dem_diff_np,out_dem_diff,dem_tif_list[0],nodata=-9999,compress='lzw',tiled='yes',bigtiff='if_safer')
    # else:
    #     # save dem diff to 16bit, centimeter, only handle diff from -327.67 to 327.67 meters
    #     dem_diff_np_cm = dem_diff_np_cm.astype(np.int16)        # save to int16
    #     raster_io.save_numpy_array_to_rasterfile(dem_diff_np_cm, out_dem_diff_cm, dem_tif_list[0],nodata=32767,compress='lzw',tiled='yes',bigtiff='if_safer')


    return True

def check_dem_diff_results(save_dir,pre_name,extent_id):

    save_dem_diff = os.path.join(save_dir, pre_name + '_DEM_diff_sub_%d.tif' % extent_id)
    save_date_diff = os.path.join(save_dir, pre_name + '_date_diff_sub_%d.tif' % extent_id)

    if os.path.isfile(save_dem_diff) and os.path.isfile(save_date_diff):
        basic.outputlogMessage('warning, DEM difference already exist, skipping create new ones')
        return True

    out_dem_diff_cm = io_function.get_name_by_adding_tail(save_dem_diff, 'cm')
    if os.path.isfile(out_dem_diff_cm):
        basic.outputlogMessage('warning, DEM difference already exist, skipping create new ones')
        return True

    return False

def main(options, args):
    extent_shp = args[0]

    dem_list_txt = options.dem_list_txt

    if b_dem_diff:
        # crop each one to the same extent, easy for DEM differnce.
        same_extent = True

    dem_tif_list = io_function.read_list_from_txt(dem_list_txt)
    # check projection
    for dem_tif in dem_tif_list:
        dem_prj = map_projection.get_raster_or_vector_srs_info_epsg(dem_tif)
        if dem_prj != extent_prj:
            raise ValueError('The projection of %s is different from %s' % (dem_prj, extent_prj))

    pass


if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp dem_dir "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: difference for multi-temporal DEM '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")


    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

    pass


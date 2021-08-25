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
import vector_gpd
import basic_src.timeTools as timeTools
import raster_io
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import split_image

import numpy as np
from itertools import combinations
import operator

from dem_mosaic_crop import subset_image_by_polygon_box
from dem_mosaic_crop import group_demTif_yearmonthDay

import multiprocessing
from multiprocessing import Pool

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

def dem_diff_new_old_min_neg_diff_patch(idx, patch, patch_count,date_pair_list_sorted,dem_groups_date):
    # conduct difference by keep the maximum subsidence (negative values)
    print('tile: %d / %d' % (idx + 1, patch_count))

    patch_w = patch[2]
    patch_h = patch[3]
    patch_date_diff = np.zeros((patch_h, patch_w), dtype=np.uint16)
    patch_dem_diff = np.empty((patch_h, patch_w), dtype=np.float32)
    patch_dem_diff[:] = 9999  # positive 9999

    # use dict to read data from disk (only need)
    dem_data_dict = {}
    for p_idx, pair in enumerate(date_pair_list_sorted):
        diff_days = (pair[1] - pair[0]).days

        data_old, data_new = read_date_dem_to_memory(p_idx, pair, date_pair_list_sorted, dem_data_dict, dem_groups_date,
                                                     boundary=patch)
        diff_two = data_new - data_old

        # fill the element
        new_ele = np.where(diff_two < patch_dem_diff )   # keep the negative values

        patch_dem_diff[new_ele] = diff_two[new_ele]
        patch_date_diff[new_ele] = diff_days

    # for locations where valid value
    patch_dem_diff[patch_dem_diff == 9999] = np.nan

    return patch, patch_dem_diff, patch_date_diff


def dem_diff_newest_oldest_a_patch(idx, patch, patch_count,date_pair_list_sorted,dem_groups_date):
    print('tile: %d / %d' % (idx + 1, patch_count))

    patch_w = patch[2]
    patch_h = patch[3]
    patch_date_diff = np.zeros((patch_h, patch_w), dtype=np.uint16)
    patch_old_date_idx = np.empty((patch_h, patch_w), dtype=np.uint8)
    patch_new_date_idx = np.empty((patch_h, patch_w), dtype=np.uint8)
    patch_old_date_idx[:] = 255
    patch_new_date_idx[:] = 255

    patch_dem_diff = np.empty((patch_h, patch_w), dtype=np.float32)
    patch_dem_diff[:] = np.nan

    date_list = [item for item in dem_groups_date.keys()]

    # use dict to read data from disk (only need)
    dem_data_dict = {}
    for p_idx, pair in enumerate(date_pair_list_sorted):
        diff_days = (pair[1] - pair[0]).days
        # basic.outputlogMessage('Getting DEM difference using the one on %s and %s, total day diff: %d' %
        #                        (timeTools.date2str(pair[1]), timeTools.date2str(pair[0]), diff_days))
        # print(pair,':',(pair[1] - pair[0]).days)

        data_old, data_new = read_date_dem_to_memory(p_idx, pair, date_pair_list_sorted, dem_data_dict, dem_groups_date,
                                                     boundary=patch)

        # print('data_old shape:',data_old.shape)
        # print('data_new shape:',data_new.shape)

        diff_two = data_new - data_old
        # print(diff_two)

        # fill the element
        new_ele = np.where(np.logical_and(np.isnan(patch_dem_diff), ~np.isnan(diff_two)))

        patch_dem_diff[new_ele] = diff_two[new_ele]
        patch_date_diff[new_ele] = diff_days
        # output the index of dates
        patch_old_date_idx[new_ele] = date_list.index(pair[0])
        patch_new_date_idx[new_ele] = date_list.index(pair[1])

        # check if all have been filled ( nan pixels)
        diff_remain_hole = np.where(np.isnan(patch_dem_diff))
        # basic.outputlogMessage(' remain %.4f percent pixels need to be filled'% (100.0*diff_remain_hole[0].size/patch_dem_diff.size) )
        if diff_remain_hole[0].size < 1:
            break

    return patch,patch_dem_diff,patch_date_diff, patch_old_date_idx,patch_new_date_idx

def dem_diff_newest_oldest(dem_tif_list, out_dem_diff, out_date_diff, process_num, b_max_subsidence=False,b_save_cm=False):
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

    # change the key to integer number after sorting and save to txt file
    dem_groups_date_sort_idx = {}
    for idx, key in enumerate(dem_groups_date.keys()):
        dem_groups_date_sort_idx[idx] = dem_groups_date[key]
    io_function.save_dict_to_txt_json(txt_save_path,dem_groups_date_sort_idx)

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
    old_date_index = np.zeros((height, width),dtype=np.uint8)
    new_date_index = np.zeros((height, width),dtype=np.uint8)
    dem_diff_np = np.empty((height, width),dtype=np.float32)
    dem_diff_np[:] = np.nan

    if process_num == 1:
        for idx, patch in enumerate(image_patches):
            _,patch_dem_diff,patch_date_diff, patch_old_date_idx,patch_new_date_idx = \
                dem_diff_newest_oldest_a_patch(idx, patch, patch_count,date_pair_list_sorted,dem_groups_date)
            # copy to the entire image
            row_s = patch[1]
            row_e = patch[1] + patch[3]
            col_s = patch[0]
            col_e = patch[0] + patch[2]
            dem_diff_np[row_s:row_e, col_s:col_e] = patch_dem_diff
            date_diff_np[row_s:row_e, col_s:col_e] = patch_date_diff
            old_date_index[row_s:row_e, col_s:col_e] = patch_old_date_idx
            new_date_index[row_s:row_e, col_s:col_e] = patch_new_date_idx
    else:
        theadPool = Pool(process_num)
        parameters_list = [ (idx, patch, patch_count,date_pair_list_sorted,dem_groups_date) for idx, patch in enumerate(image_patches)]
        if b_max_subsidence is False:
            results = theadPool.starmap(dem_diff_newest_oldest_a_patch, parameters_list)
        else:
            results = theadPool.starmap(dem_diff_new_old_min_neg_diff_patch , parameters_list)
        for res in results:
            patch, patch_dem_diff, patch_date_diff,patch_old_date_idx,patch_new_date_idx = res
            # copy to the entire image
            row_s = patch[1]
            row_e = patch[1] + patch[3]
            col_s = patch[0]
            col_e = patch[0] + patch[2]
            dem_diff_np[row_s:row_e, col_s:col_e] = patch_dem_diff
            date_diff_np[row_s:row_e, col_s:col_e] = patch_date_diff
            old_date_index[row_s:row_e, col_s:col_e] = patch_old_date_idx
            new_date_index[row_s:row_e, col_s:col_e] = patch_new_date_idx

    # save date diff to tif (16 bit)
    raster_io.save_numpy_array_to_rasterfile(date_diff_np,out_date_diff,dem_tif_list[0], nodata=0,compress='lzw',tiled='yes',bigtiff='if_safer')
    # save old and new date index to tif (8 bit)
    out_old_date_idx = io_function.get_name_by_adding_tail(out_date_diff,'oldIndex')
    out_new_date_idx = io_function.get_name_by_adding_tail(out_date_diff,'newIndex')
    raster_io.save_numpy_array_to_rasterfile(old_date_index,out_old_date_idx,dem_tif_list[0], nodata=255,compress='lzw',tiled='yes',bigtiff='if_safer')
    raster_io.save_numpy_array_to_rasterfile(new_date_index,out_new_date_idx,dem_tif_list[0], nodata=255,compress='lzw',tiled='yes',bigtiff='if_safer')

    # # stretch the DEM difference, save to 8 bit.
    # dem_diff_np_8bit = raster_io.image_numpy_to_8bit(dem_diff_np,10,-10,dst_nodata=0)
    # out_dem_diff_8bit = io_function.get_name_by_adding_tail(out_dem_diff, '8bit')
    # raster_io.save_numpy_array_to_rasterfile(dem_diff_np_8bit, out_dem_diff_8bit, dem_tif_list[0], nodata=0)


    # if possible, save to 16 bit, to save the disk storage.
    # dem_diff_np[0:5,0] = -500
    # dem_diff_np[0,0:5] = 500
    # print(np.nanmin(dem_diff_np))
    # print(np.nanmax(dem_diff_np))

    # if np.nanmin(dem_diff_np_cm) < range.min or np.nanmax(dem_diff_np_cm) > range.max:
    # save dem diff to files (float), meter
    if b_save_cm is False:
        raster_io.save_numpy_array_to_rasterfile(dem_diff_np,out_dem_diff,dem_tif_list[0],nodata=-9999,compress='lzw',tiled='yes',bigtiff='if_safer')
    else:
        # save dem diff to 16bit, centimeter, only handle diff from -327.67 to 327.67 meters
        bit16_nodata = 32767
        range = np.iinfo(np.int16)
        dem_diff_np_cm = dem_diff_np*100
        dem_diff_np_cm[dem_diff_np_cm < range.min] = range.min
        dem_diff_np_cm[dem_diff_np_cm > range.max] = range.max
        dem_diff_np_cm[np.isnan(dem_diff_np_cm)] = bit16_nodata  # set the nodata for int16
        dem_diff_np_cm = dem_diff_np_cm.astype(np.int16)        # save to int16
        out_dem_diff_cm = out_dem_diff
        basic.outputlogMessage('note, save DEM difference (%s) to centimeter, int16, range: -327.68 to 327.67 m'%os.path.basename(out_dem_diff_cm))
        raster_io.save_numpy_array_to_rasterfile(dem_diff_np_cm, out_dem_diff_cm, dem_tif_list[0],nodata=bit16_nodata,compress='lzw',tiled='yes',bigtiff='if_safer')


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

def crop_to_same_exent_for_diff(dem_tif_list, save_dir, extent_id, extent_poly,process_num):
    # crop to the same extent
    crop_tif_dir = os.path.join(save_dir, 'dem_crop_for_diff_sub_%d' % extent_id)
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
                                                     same_extent=True,thread_num=process_num)
            if crop_tif is False:
                # raise ValueError('warning, crop %s failed' % tif)
                continue
            crop_tif_list.append(crop_tif)
    dem_tif_list = crop_tif_list

    return dem_tif_list

def main(options, args):

    save_dir = options.save_dir
    extent_shp = options.extent_shp
    process_num = options.process_num

    dem_dir_or_txt = args[0]
    if os.path.isfile(dem_dir_or_txt):
        dem_list = io_function.read_list_from_txt(dem_dir_or_txt)
    else:
        dem_list = io_function.get_file_list_by_ext('.tif', dem_dir_or_txt, bsub_folder=False)
    dem_count = len(dem_list)
    if dem_count < 1:
        raise ValueError('No input dem files in %s' % dem_dir_or_txt)

    if extent_shp is not None:
        pre_name = os.path.splitext(os.path.basename(extent_shp))[0]
    else:
        pre_name = os.path.basename(os.path.abspath(save_dir))
    save_dem_diff = os.path.join(save_dir, pre_name + '_DEM_diff.tif')
    save_date_diff = os.path.join(save_dir, pre_name + '_date_diff.tif')
    if os.path.isfile(save_dem_diff) and os.path.isfile(save_date_diff):
        print('%s and %s exists, skip'%(save_dem_diff, save_date_diff))
        return

    if extent_shp is not None:
        # crop the DEM before differencing
        extent_shp_base = os.path.splitext(os.path.basename(extent_shp))[0]
        dem_prj = map_projection.get_raster_or_vector_srs_info_epsg(dem_list[0])
        extent_prj = map_projection.get_raster_or_vector_srs_info_epsg(extent_shp)
        if dem_prj != extent_prj:
            raise ValueError('The projection of extent file (%s) and dem tifs is different'%extent_shp)

        extent_polys = vector_gpd.read_polygons_gpd(extent_shp)
        if len(extent_polys) != 1:
            raise ValueError('Only allow one polygon in %s' % extent_shp)

        extPolys_ids = vector_gpd.read_attribute_values_list(extent_shp, 'id')
        if extPolys_ids is None or None in extPolys_ids:
            basic.outputlogMessage('Warning, field: id is not in %s, will create default ID for each grid' % extent_shp)
            extPolys_ids = [id + 1 for id in range(len(extent_polys))]

        # crop
        for idx, ext_poly in zip(extPolys_ids, extent_polys):
            basic.outputlogMessage('crop and differnce DEM for the %d th extent (%d in total)' % (idx, len(extent_polys)))
            crop_dem_list = crop_to_same_exent_for_diff(dem_list, save_dir, idx, ext_poly, process_num)

            dem_list = crop_dem_list

    dem_diff_newest_oldest(dem_list, save_dem_diff, save_date_diff,process_num, b_max_subsidence=options.max_subsidence)



if __name__ == '__main__':
    usage = "usage: %prog [options] dem_tif_dir or dem_list_txt "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: difference for multi-temporal DEM '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    parser.add_option("-e", "--extent_shp",
                      action="store", dest="extent_shp",
                      help="the extent file for cropping")

    parser.add_option("-m", "--max_subsidence",
                      action="store_true", dest="max_subsidence",default=False,
                      help="for each pixel, keep the maximum elevation reduction values")


    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

    pass


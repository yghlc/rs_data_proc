#!/usr/bin/env python
# Filename: dem_headwall_medial_axis.py 
"""
introduction: try to extract headwall based on skimage.morphology.medial_axis

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 26 August, 2021
"""


import os,sys
from optparse import OptionParser
import time

from skimage import morphology
import numpy as np
import pandas as pd

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import vector_gpd
import raster_io
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection
import basic_src.basic as basic


def slope_tif_to_slope_bin(slope_tif,slope_bin_path,slope_threshold):
    if os.path.isfile(slope_bin_path):
        print('%s exist'%slope_bin_path)
    else:
        slope_data, nodata = raster_io.read_raster_one_band_np(slope_tif)
        bin_slope = np.zeros_like(slope_data,dtype=np.uint8)
        bin_slope[slope_data > slope_threshold] = 1
        bin_slope[slope_data > 88] = 0          # if slope is too large, it may caused by artifacts, so remove them

        # save
        slope_bin = bin_slope*255
        #  set nodata as 0
        if raster_io.save_numpy_array_to_rasterfile(slope_bin,slope_bin_path,slope_tif,nodata=0,compress='lzw',tiled='yes',bigtiff='if_safer') is not True:
            return None

    return slope_bin_path

def slope_bin_to_medial_axis_raster(in_image_path, out_image_path):

    if os.path.isfile(out_image_path):
        print('%s exists, skip slope_bin_to_medial_axis_raster'%out_image_path)
        return out_image_path


    image_np, nodata = raster_io.read_raster_one_band_np(in_image_path)

    out_np = morphology.medial_axis(image_np, mask=None, return_distance=False)
    # out_np,dist = morphology.medial_axis(image_np, mask=None, return_distance=True)
    # bool numpy to uint8
    out_np = out_np.astype(np.uint8)
    # dist = dist.astype(np.float32)

    # save to file
    raster_io.save_numpy_array_to_rasterfile(out_np, out_image_path, in_image_path, compress='lzw',
                                             tiled='yes', bigtiff='if_safer')
    # save distances to file (no need)
    # out_dist_path = io_function.get_name_by_adding_tail(out_image_path,'dist')
    # raster_io.save_numpy_array_to_rasterfile(dist, out_dist_path, in_image_path, compress='lzw',
    #                                          tiled='yes', bigtiff='if_safer')

    return out_image_path

def medial_axis_raster_to_vector(in_medial_axis_tif, out_vector_shp, raster_res=2):

    if vector_gpd.raster2shapefile(in_medial_axis_tif,out_shp=out_vector_shp,connect8=True) is None:
        return None

    vector_shp_buff = io_function.get_name_by_adding_tail(out_vector_shp,'buff')
    if os.path.isfile(vector_shp_buff):
        print('%s exists, skip buffering'%vector_shp_buff)
        return vector_shp_buff


    polys = vector_gpd.read_polygons_gpd(out_vector_shp,b_fix_invalid_polygon=False)

    # calculate the attributes before buffer
    # from the area, we can tell how many pixels in each line (line segment), each pixel have size of 2*2 m^2
    id_list = [item for item in range(len(polys))]
    area_list = [ item.area for item in polys]
    length_list = [ item.length for item in polys]

    pixel_size = raster_res*raster_res
    pixel_count_list = [ item/pixel_size for item in area_list]

    # buffer 0.1 meters, so the width of polygons is around 2.02 meters (2-m ArcticDEM), still have width of one pixel
    polys_buff = [item.buffer(0.01) for item in polys]

    # after buffer, for complex medial axis, it may have holes in the polygons, get hole count
    hole_count_list =[len(list(item.interiors)) for item in polys_buff]

    # save, overwrite out_vector_shp
    wkt = map_projection.get_raster_or_vector_srs_info_proj4(out_vector_shp)
    save_pd = pd.DataFrame({'id':id_list,'area':area_list,'length':length_list,'pixels':pixel_count_list,
                            'holes':hole_count_list,'Polygon':polys_buff})
    vector_gpd.save_polygons_to_files(save_pd, 'Polygon', wkt, vector_shp_buff)

    return vector_shp_buff


def remove_based_on_length_pixel(medial_axis_shp,min_length, max_length,wkt, rm_length_shp):
    polygons, lengths_pixel = vector_gpd.read_polygons_attributes_list(medial_axis_shp,'pixels',b_fix_invalid_polygon=False)

    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, (poly,length) in enumerate(zip(polygons,lengths_pixel)):
        # remove too long or too short ones
        if length > max_length or length < min_length:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on length in pixel, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_length_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_length_shp)

    return rm_length_shp


def remove_based_on_hole(medial_axis_shp, max_hole,wkt, rm_hole_shp):

    polygons, holes_count = vector_gpd.read_polygons_attributes_list(medial_axis_shp,'holes',b_fix_invalid_polygon=False)

    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, (poly,holes) in enumerate(zip(polygons,holes_count)):
        # remove too long or too short ones
        if holes > max_hole:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on holes, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_hole_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_hole_shp)

    return rm_hole_shp

def calculate_remove_based_on_line_segments(medial_axis_shp, max_line_segments,wkt, rm_line_segment_shp):

    # calculate the number of line segments
    polygons = vector_gpd.read_polygons_gpd(medial_axis_shp,b_fix_invalid_polygon=False)
    line_segments_list = []
    for idx, poly in enumerate(polygons):
        out_line = list(poly.exterior)
        in_lines = list(poly.interiors)
        count = len(out_line) + len(in_lines)
        line_segments_list.append(count)
    add_attributes = {'lines':line_segments_list}
    vector_gpd.add_attributes_to_shp(medial_axis_shp,add_attributes)

    # remove based on the number of line segments
    remain_polygons_idx = []
    # remove relative large but narrow ones.
    remove_count = 0
    for idx, lines in enumerate(line_segments_list):
        # remove too long or too short ones
        if lines > max_line_segments:
            remove_count += 1
            continue
        remain_polygons_idx.append(idx)

    basic.outputlogMessage('remove %d polygons based on the count of line segments, remain %d ones saving to %s' %
                           (remove_count, len(remain_polygons_idx), rm_line_segment_shp))

    if len(remain_polygons_idx) < 1:
        return False

    vector_gpd.save_shapefile_subset_as(remain_polygons_idx,medial_axis_shp,rm_line_segment_shp)

    return rm_line_segment_shp



def extract_headwall_based_medial_axis_from_slope(idx, total, slope_tif, work_dir, save_dir,slope_threshold,
                                                  min_length, max_length, max_hole_count,process_num):
    '''
    extract headwall from slope based on medial axis (skimage.morphology.medial_axis)
    :param idx: tif index
    :param total: total slope file count
    :param slope_tif: slope file
    :param work_dir:
    :param save_dir:
    :param slope_threshold:
    :param min_length: min length, the medial axis calculated from skimage.morphology has width of one pixel, the length is based pixel count of line segments
    :param max_length: max length (pixel count)
    :param max_hole_count: some complex line segment may end in holes when forming polygons
    :param process_num:
    :return:
    '''

    headwall_shp = os.path.splitext(os.path.basename(io_function.get_name_by_adding_tail(slope_tif, 'headwall')))[0] + '.shp'
    save_headwall_shp = os.path.join(save_dir, headwall_shp)
    if os.path.isfile(save_headwall_shp):
        print('%s exists, skip' % save_headwall_shp)
        return save_headwall_shp

    print('(%d/%d) extracting headwall from %s' % (idx, total, slope_tif))

    wkt = map_projection.get_raster_or_vector_srs_info_wkt(slope_tif)
    # binary slope
    slope_bin_path = os.path.join(work_dir, os.path.basename(io_function.get_name_by_adding_tail(slope_tif, 'bin')))
    if slope_tif_to_slope_bin(slope_tif, slope_bin_path, slope_threshold) is None:
        return False

    # get medial axis raster
    medial_axis_tif = io_function.get_name_by_adding_tail(slope_bin_path,'medial_axis')
    if slope_bin_to_medial_axis_raster(slope_bin_path, medial_axis_tif) is None:
        return False

    # get madial axis vector (polygons: width of one pixel)
    medial_axis_poly_shp = os.path.join(work_dir, io_function.get_name_no_ext(medial_axis_tif) + '_poly.shp')
    medial_axis_poly_shp_buff = medial_axis_raster_to_vector(medial_axis_tif, medial_axis_poly_shp)
    if medial_axis_poly_shp_buff is None:
        return False

    # only keep not too long or too short line segments
    rm_length_shp = io_function.get_name_by_adding_tail(medial_axis_poly_shp_buff, 'rmLength')
    if os.path.isfile(rm_length_shp):
        print('%s exists, skip removing based on length in pixels' % rm_length_shp)
    else:
        # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
        if remove_based_on_length_pixel(medial_axis_poly_shp_buff, min_length, max_length, wkt, rm_length_shp) is False:
            return False

    # remove based on  hole count, if too many holes, not headwall
    rm_hole_shp = io_function.get_name_by_adding_tail(rm_length_shp, 'rmHole')
    if os.path.isfile(rm_hole_shp):
        print('%s exists, skip removing based holes' % rm_hole_shp)
    else:
        # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
        if remove_based_on_hole(rm_length_shp, max_hole_count, wkt, rm_hole_shp) is False:
            return False

    # get line segments in polygons and remove based on line segments
    max_line_count = 10
    rm_line_shp = io_function.get_name_by_adding_tail(rm_hole_shp, 'rmLine')
    if os.path.isfile(rm_line_shp):
        print('%s exists, skip removing based the count of line segments' % rm_line_shp)
    else:
        # medial_axis_shp,min_length, max_length,wkt, rm_length_shp
        if calculate_remove_based_on_line_segments(rm_hole_shp, max_line_count, wkt, rm_line_shp) is False:
            return False



    # # add some shape info
    # rm_shapeinfo_shp = io_function.get_name_by_adding_tail(slope_bin_shp, 'rmShape')
    # if os.path.isfile(rm_shapeinfo_shp):
    #     print('%s exists, skip removing based on shape' % rm_shapeinfo_shp)
    # else:
    #     if remove_based_on_shapeinfo(rm_area_shp, rm_shapeinfo_shp, max_box_WH) is False:
    #         return False
    #
    #
    # # copy the results.
    # io_function.copy_shape_file(rm_medialAxis_shp, save_headwall_shp)



def test_slope_bin_to_medial_axis():
    data_dir = os.path.expanduser('~/Data/dem_processing/headwall_shp_sub_6174/20080511_dem_slope')
    slope_bin_path = os.path.join(data_dir, '20080511_dem_slope_bin.tif' )
    medial_axis_tif = os.path.join(data_dir, '20080511_dem_slope_bin_medial_axis.tif')

    medial_axis_tif = slope_bin_to_medial_axis_raster(slope_bin_path,medial_axis_tif)

    medial_axis_poly_shp = os.path.join(data_dir, '20080511_dem_slope_bin_medial_axis_poly.shp')
    medial_axis_raster_to_vector(medial_axis_tif, medial_axis_poly_shp)

def test_extract_headwall_based_medial_axis_from_slope():

    data_dir = os.path.expanduser('~/Data/dem_processing/grid_6174_tmp_files/slope_sub_6174')
    slope_tif = os.path.join(data_dir,'20080511_dem_slope.tif')
    work_dir = os.path.expanduser('~/Data/dem_processing')
    save_dir = os.path.expanduser('~/Data/dem_processing/grid_6174_tmp_files')
    slope_threshold = 20
    min_length = 6
    max_length = 500
    max_hole_count = 0
    process_num = 1

    extract_headwall_based_medial_axis_from_slope(0, 1, slope_tif, work_dir, save_dir, slope_threshold,
                                                  min_length, max_length, max_hole_count, process_num)


def main():
    # test_slope_bin_to_medial_axis()
    test_extract_headwall_based_medial_axis_from_slope()
    pass


if __name__ == '__main__':
    main()
    pass
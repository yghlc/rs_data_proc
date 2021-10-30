#!/usr/bin/env python
# Filename: mosaic_images_crop_grid 
"""
introduction: create mosaic of many small images and crop to grid

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 19 July, 2020
"""
import sys,os
from optparse import OptionParser
import re

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import vector_gpd

import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.RSImageProcess as RSImageProcess
import basic_src.timeTools as timeTools
from datetime import datetime

import time

import numpy as np

import operator

import multiprocessing
from multiprocessing import Pool

code_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',)
sys.path.insert(0, os.path.join(code_dir,'planetScripts'))

# sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/ChangeDet_DL/dataTools'))
from  get_planet_image_list import  get_Planet_SR_image_list_overlap_a_polygon

prePlanetImage = os.path.join(code_dir,'planetScripts','prePlanetImage.py')

# some issues for the global var in the multiple processes, so remove the temporal foolder in bash
temporal_dirs = []

def convert_planet_to_rgb_images(tif_path,save_dir='RGB_images', sr_min=0, sr_max=3000, save_org_dir=None, sharpen=True, rgb_nodata=0):

    #if multiple processes try to derive the same rgb images, it may have problem.
    # save output to 'RGB_images' + processID

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    if save_org_dir is not None and os.path.isdir(save_org_dir) is False:
        io_function.mkdir(save_org_dir)

    if save_org_dir is not None:
        copied_org_img_path = os.path.join(save_org_dir,os.path.basename(tif_path))
        io_function.copy_file_to_dst(tif_path,copied_org_img_path)

    # filename_no_ext
    output = os.path.splitext(os.path.basename(tif_path))[0]
    if sharpen:
        fin_output= os.path.join(save_dir, output + '_8bit_rgb_sharpen.tif')
    else:
        fin_output = os.path.join(save_dir, output + '_8bit_rgb.tif')
    if os.path.isfile(fin_output):
        basic.outputlogMessage("Skip, because File %s exists in current folder: %s"%(fin_output,os.getcwd()))
        return fin_output

    # use fix min and max to make the color be consistent to sentinel-images
    src_min=sr_min
    src_max=sr_max
    dst_min=1       # 0 is the nodata, so set as 1
    dst_max=255

    # gdal_translate -ot Byte -scale ${src_min} ${src_max} ${dst_min} ${dst_max} ${image_path} ${output}_8bit.tif
    if 'SR.tif' in tif_path:
        cmd_str = 'gdal_translate -ot Byte -scale %d %d %d %d -of VRT %s %s_8bit.tif'%(src_min,src_max,dst_min,dst_max,tif_path,output)
    else:
        # gdal_contrast_stretch -percentile-range 0.01 0.99 ${output}.tif ${output}_8bit.tif
        cmd_str = 'gdal_contrast_stretch -percentile-range 0.01 0.99 %s %s_8bit.tif' % (tif_path, output)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)

    # the third band is red, second is green, and first is blue
    #gdal_translate -b 3 -b 2 -b 1  ${output}_8bit.tif ${output}_8bit_rgb.tif
    cmd_str = 'gdal_translate -b 3 -b 2 -b 1 -of VRT %s_8bit.tif %s_8bit_rgb.tif'%(output,output)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)

    # python ${code_dir}/planetScripts/prePlanetImage.py ${output}_8bit_rgb.tif ${fin_output}
    if sharpen:
        cmd_str = 'python %s %s_8bit_rgb.tif %s'%(prePlanetImage,output,fin_output)
    else:
        # convert from VRT format to tif format
        cmd_str = 'gdal_translate -of GTiff %s_8bit_rgb.tif %s' % (output, fin_output)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)



    # set nodata
    # gdal_edit.py -a_nodata 0  ${fin_output}
    cmd_str = 'gdal_edit.py -a_nodata %d  %s' % (rgb_nodata, fin_output)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)

    io_function.delete_file_or_dir('%s_8bit.tif'%output)
    io_function.delete_file_or_dir('%s_8bit_rgb.tif'%output)

    return fin_output

def reproject_planet_image(tif_path, new_prj_wkt, new_prj_proj4, save_dir='planet_images_reproj'):
    '''
    reprojection of images
    :param tif_path: image path
    :param new_prj_wkt: new projection in wkt format (more accurate)
    :param new_prj_proj4: new projection in proj format (not accurate, but good for comparision)
    :param save_dir: output save folder.
    :return:
    '''

    # if multiple processes try to derive the same rgb images, it may have problem.
    # save output to 'planet_images_reproj' + processID

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    # filename_no_ext
    output = os.path.splitext(os.path.basename(tif_path))[0]
    fin_output= os.path.join(save_dir, output + '_prj.tif')
    if os.path.isfile(fin_output):
        basic.outputlogMessage("Skip, because File %s exists in current folder: %s"%(fin_output,os.getcwd()))
        return fin_output

    tif_prj4 = map_projection.get_raster_or_vector_srs_info_proj4(tif_path).strip()
    # if they have the same projection, then return False, no need to reproject
    if tif_prj4 == new_prj_proj4:
        return False

    # reproject to the new projection
    # gdalwarp -t_srs EPSG:4326  -overwrite tmp.tif $out
    cmd_str = 'gdalwarp -t_srs %s -of VRT %s %s'%(new_prj_wkt,tif_path,fin_output)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)

    return fin_output

def create_moasic_of_each_grid_polygon(id,polygon, polygon_latlon, out_res, cloud_cover_thr, geojson_list, save_dir,
                                       new_prj_wkt=None,new_prj_proj4=None, sr_min=0, sr_max=3000,to_rgb=True, nodata=0, save_org_dir=None,
                                       resampling_method='min'):
    '''
    create mosaic for Planet images within a grid
    :param polygon:
    :param polygon_latlon:
    :param out_res:
    :param cloud_cover_thr:
    :param geojson_list:
    :param save_dir:
    :param new_prj_wkt:
    :param new_prj_proj4:
    :param sr_min:
    :param sr_max:
    :param to_rgb:
    :param nodata:
    :return:
    '''
    time0 = time.time()
    file_name = os.path.basename(save_dir)
    fin_out = os.path.join(save_dir, file_name + '_sub_%d.tif' % id)
    if os.path.isfile(fin_out):
        basic.outputlogMessage('Warning, skip %s because it already exists, remove it if want to regenerate it'%fin_out)
        return fin_out

    # get image list and cloud cover
    planet_img_list, cloud_covers = get_Planet_SR_image_list_overlap_a_polygon(polygon_latlon,geojson_list,cloud_cover_thr)
    if len(planet_img_list) < 1:
        basic.outputlogMessage('warning, no images within %d grid'%id)
        return False

    io_function.mkdir(save_dir)

    print('images and their cloud cover for %dth grid:'%id)
    for img, cloud_cover in zip(planet_img_list, cloud_covers):
        print(img, cloud_cover)

    proc_id = multiprocessing.current_process().pid

    # convert to RGB images (for Planet)
    rgb_image_list = []
    rgb_dir = 'RGB_images_'+str(proc_id)
    if to_rgb:
        for tif_path in planet_img_list:
            rgb_img = convert_planet_to_rgb_images(tif_path,save_dir=rgb_dir,save_org_dir=save_org_dir, sr_min=sr_min, sr_max=sr_max)
            rgb_image_list.append(rgb_img)
    if len(rgb_image_list) > 0:
        planet_img_list = rgb_image_list

    reproj_img_list = []
    # reproject if necessary
    reproj_dir = 'planet_images_reproj_' + str(proc_id)
    if new_prj_wkt != None and new_prj_proj4 != None:
        for tif_path in planet_img_list:
            prj_out = reproject_planet_image(tif_path, new_prj_wkt, new_prj_proj4, save_dir=reproj_dir)
            # replace the image
            if prj_out is not False and os.path.isfile(prj_out):
                reproj_img_list.append(prj_out)
            else:
                # if not reproject, then append the original image.
                reproj_img_list.append(tif_path)
    if len(reproj_img_list) > 0:
        planet_img_list = reproj_img_list

    # create mosaic using gdal_merge.py
    # because in gdal_merge.py, a later image will replace one, so we put image with largest cloud cover first

    out = os.path.join(save_dir,file_name + '_sub_%d_tmp.tif'%id)
    if os.path.isfile(out):
        io_function.delete_file_or_dir(out)

    # reverse=True to make it in descending order
    img_cloud_list = [(img_path,cloud) for cloud, img_path in sorted(zip(cloud_covers,planet_img_list), key=lambda pair: pair[0],reverse=True)]
    # for checking
    print('Image and its cloud after sorting:')
    for (img_path,cloud)  in img_cloud_list:
        print(img_path,cloud)
    tifs = [img_path for (img_path,cloud)  in img_cloud_list ]
    tifs_str = ' '.join(tifs)

    # cmd_str = 'gdal_merge.py -o %s -n %d -init %d -ps %d %d %s'%(out,nodata,nodata,out_res,out_res,tifs_str)
    cmd_str = 'gdalbuildvrt -resolution user -tr %d %d -srcnodata %d -vrtnodata %d  %s %s'%(out_res,out_res,nodata,nodata,out,tifs_str)
    status, result = basic.exec_command_string(cmd_str)
    if status != 0:
        print(result)
        sys.exit(status)

    # # #  polygon.exterior.coords
    # minx, miny, maxx, maxy =  polygon.bounds    # (minx, miny, maxx, maxy)
    # print(minx, miny, maxx, maxy)
    # results = RSImageProcess.subset_image_projwin(fin_out,out,minx, maxy, maxx, miny, xres=out_res,yres=out_res)
    # print(results)
    results = RSImageProcess.subset_image_by_polygon_box_image_min(fin_out,out,polygon, xres=out_res,yres=out_res,
                                                                   compress='lzw',tiled='yes',bigtiff='if_safer')

    if results is False:
        basic.outputlogMessage('Warning, Crop %s failed, keep the one without cropping'%out)
        io_function.move_file_to_dst(out,fin_out)
    else:
        io_function.delete_file_or_dir(out)

    # ## mosaic and crop at the same time together
    # minx, miny, maxx, maxy =  polygon.bounds    # (minx, miny, maxx, maxy)
    # print(minx, miny, maxx, maxy)
    # results = RSImageProcess.mosaic_crop_images_gdalwarp(tifs,fin_out,src_nodata=nodata,min_x=minx,min_y=miny,max_x=maxx,max_y=maxy,
    #                                                      xres=out_res,yres=out_res,resampling_method=resampling_method)
    #
    # if results is False:
    #     basic.outputlogMessage('Warning, create %s failed' % fin_out)
    #     return False

    # sys.exit(0)
    cost_time_sec = time.time() - time0
    basic.outputlogMessage('finished creating %s cost %.2f seconds (%.2f minutes)' % (fin_out,cost_time_sec,cost_time_sec/60))

    return fin_out

def create_moasic_of_each_grid_polygon_one_proc(parameters):

    return create_moasic_of_each_grid_polygon(parameters[0], parameters[1], parameters[2], parameters[3],
                                       parameters[4], parameters[5], parameters[6],
                                       new_prj_wkt=parameters[7], new_prj_proj4=parameters[8],
                                       sr_min=parameters[9], sr_max=parameters[10],
                                       save_org_dir=parameters[11],resampling_method=parameters[12])

def group_planet_images_date(geojson_list,diff_days=0):
    '''
    group image based on their acquisition date
    :param geojson_list:
    :param diff_days:
    :return: a dict
    '''

    # e.g., 20200830_222236_71_1057.geojson
    img_groups = {}
    for item in geojson_list:
        yeardate =  timeTools.get_yeardate_yyyymmdd(os.path.basename(item))

        b_assgined = False
        for time in img_groups.keys():
            if timeTools.diff_yeardate(time,yeardate) <= diff_days:
                img_groups[time].append(item)
                b_assgined = True
                break
        if b_assgined is False:
            img_groups[yeardate] = [item]


    # convert the key from datetime to string.
    strKey_dict = {}
    for key in img_groups.keys():
        # key.strftime(key,'%Y%m%d')
        strKey_dict[datetime.strftime(key,'%Y%m%d')] = img_groups[key]

    return strKey_dict










def main(options, args):

    time0 = time.time()
    image_dir = args[0]
    geojson_list = io_function.get_file_list_by_ext('.geojson',image_dir,bsub_folder=False)
    # remove some scenes, or maybe we should set bsub_folder=False
    # geojson_list = [item for item in geojson_list if 'incomplete_scenes' not in item ]  # remove those in "incomplete_scenes"
    # geojson_list = [item for item in geojson_list if 'scenes_high_cloud_cover' not in item ]  # remove those in "scenes_high_cloud_cover"

    if len(geojson_list) < 1:
        raise ValueError('There is no geojson files in %s'%image_dir)

    basic.outputlogMessage('Image Dir: %s'%image_dir)
    basic.outputlogMessage("Number of geojson files: %d" % len(geojson_list))

    grid_polygon_shp = args[1]      # the polygon should be in projection Cartesian coordinate system (e.g., UTM )
    basic.outputlogMessage('Image grid polygon shapefile: %s' % grid_polygon_shp)
    process_num = options.process_num
    basic.outputlogMessage('The number of processes for creating the mosaic is: %d' % process_num)

    # read grid polygons
    grid_polygons = vector_gpd.read_polygons_gpd(grid_polygon_shp)
    grid_ids = vector_gpd.read_attribute_values_list(grid_polygon_shp,'id')
    if grid_ids is None:
        basic.outputlogMessage('Warning, field: id is not in %s, will create default ID for each grid'%grid_polygon_shp)
        grid_ids = [ id + 1 for id in range(len(grid_polygons))]

    shp_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_polygon_shp).strip()
    # print(shp_prj)
    grid_polygons_latlon = grid_polygons
    if shp_prj != '+proj=longlat +datum=WGS84 +no_defs':
        # read polygons and reproject to 4326 projection
        grid_polygons_latlon = vector_gpd.read_shape_gpd_to_NewPrj(grid_polygon_shp,'EPSG:4326')
    # else:
    #     raise ValueError(' %s should be in projection of Cartesian coordinate system'%grid_polygon_shp)

    shp_prj_wkt = map_projection.get_raster_or_vector_srs_info_wkt(grid_polygon_shp)

    max_sr = options.max_sr
    min_sr = options.min_sr

    original_img_copy_dir = options.original_img_copy_dir
    b_to_rgb_8bit = options.to_rgb
    basic.outputlogMessage('Convert to 8bit RGB images: %s'%str(b_to_rgb_8bit))

    # group planet image based on acquisition date
    b_group_date = options.group_date
    basic.outputlogMessage('Group Planet image based on acquisition date: %s' % str(b_group_date))
    if b_group_date:
        # diff_days as 0, group images acquired at the same date
        geojson_groups = group_planet_images_date(geojson_list, diff_days=0)

        # sort based on yeardate in accending order : operator.itemgetter(0)
        geojson_groups = dict(sorted(geojson_groups.items(), key=operator.itemgetter(0)))

        save_group_txt = 'geojson_groups_input_folder.txt'
        basic.outputlogMessage('images are divided into %d groups, save to %s' % (len(geojson_groups.keys()),save_group_txt))
        io_function.save_dict_to_txt_json(save_group_txt, geojson_groups)
    else:
        geojson_groups = {'all':geojson_list}

    # create mosaic of each grid
    cloud_cover_thr = options.cloud_cover
    cloud_cover_thr = cloud_cover_thr * 100         # for Planet image, it is percentage
    out_res = options.out_res
    cur_dir = os.getcwd()
    resampling_method = options.merged_method

    for key in geojson_groups.keys():

        # # test
        # if key != '20200701':
        #     continue

        geojson_list = geojson_groups[key]
        save_dir = os.path.basename(cur_dir) + '_mosaic_' + str(out_res) + '_' + key
        # print(save_dir)
        if process_num == 1:
            for id, polygon, poly_latlon in zip(grid_ids,grid_polygons,grid_polygons_latlon):
                # if id != 34:
                #     continue
                create_moasic_of_each_grid_polygon(id, polygon, poly_latlon, out_res,
                                                   cloud_cover_thr, geojson_list,save_dir,
                                                   new_prj_wkt=shp_prj_wkt, new_prj_proj4=shp_prj,
                                                   sr_min=min_sr, sr_max=max_sr,
                                                   to_rgb = b_to_rgb_8bit,
                                                   save_org_dir=original_img_copy_dir,
                                                   resampling_method=resampling_method)
        elif process_num > 1:
            theadPool = Pool(process_num)  # multi processes

            parameters_list = [
                (id, polygon, poly_latlon, out_res,cloud_cover_thr, geojson_list,save_dir,shp_prj_wkt,shp_prj,min_sr,max_sr,b_to_rgb_8bit,0,original_img_copy_dir) for
                id, polygon, poly_latlon in zip(grid_ids,grid_polygons,grid_polygons_latlon)]
            results = theadPool.starmap(create_moasic_of_each_grid_polygon, parameters_list)  # need python3
            theadPool.close()
        else:
            raise ValueError('incorrect process number: %d'% process_num)



    cost_time_sec = time.time() - time0
    basic.outputlogMessage('Done, total time cost %.2f seconds (%.2f minutes or %.2f hours)' % (cost_time_sec,cost_time_sec/60,cost_time_sec/3600))


    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] image_dir polygon_shp "
    parser = OptionParser(usage=usage, version="1.0 2020-07-19")
    parser.description = 'Introduction: create mosaic of Planet images and save to files '
    # parser.add_option("-s", "--start_date",default='2018-04-30',
    #                   action="store", dest="start_date",
    #                   help="start date for inquiry, with format year-month-day, e.g., 2018-05-23")
    # parser.add_option("-e", "--end_date",default='2018-06-30',
    #                   action="store", dest="end_date",
    #                   help="the end date for inquiry, with format year-month-day, e.g., 2018-05-23")
    parser.add_option("-u", "--upper_sr",
                      action="store", dest="max_sr", type=float,default=3000,
                      help="the upper limit of surface reflectance (maximum of surface reflectance)")
    parser.add_option("-l", "--lower_sr",
                      action="store", dest="min_sr", type=float,default=0,
                      help="the lower limit of surface reflectance (minimum of surface reflectance) ")

    parser.add_option("-c", "--cloud_cover",
                      action="store", dest="cloud_cover", type=float,default=0.3,
                      help="the could cover threshold, only accept images with cloud cover less than the threshold")
    parser.add_option("-r", "--out_res",
                      action="store", dest="out_res", type=float,default=30,
                      help="the output resolution of mosaic")
    parser.add_option("-o", "--original_img_copy_dir",
                      action="store", dest="original_img_copy_dir",
                      help="the folder to copy and save original images")

    parser.add_option("-m", "--merged_method",
                      action="store", dest="merged_method", default='min',
                      help="the method to merge pixels at the same location, such as min, max, or med, etc")

    parser.add_option("-p", "--process_num",
                      action="store", dest="process_num",type=int,default=1,
                      help="number of processes to create the mosaic")

    parser.add_option("-t", "--to_rgb",
                      action="store_true", dest="to_rgb",default=False,
                      help="true to convert all images to 8 bit rgb")

    parser.add_option("-g", "--group_date",
                      action="store_true", dest="group_date",default=False,
                      help="true to group image if their acquisition date is the same")

    # parser.add_option("-i", "--item_types",
    #                   action="store", dest="item_types",default='PSScene4Band',
    #                   help="the item types, e.g., PSScene4Band,PSOrthoTile")
    # parser.add_option("-a", "--planet_account",
    #                   action="store", dest="planet_account",default='huanglingcao@link.cuhk.edu.hk',
    #                   help="planet email account, e.g., huanglingcao@link.cuhk.edu.hk")



    (options, args) = parser.parse_args()
    # print(options.to_rgb)
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    basic.setlogfile('mosaic_images_crop_grid_%s.log'%str(datetime.date(datetime.now())))

    main(options, args)
#!/usr/bin/env python
# Filename: create_vrt_mosaic_for_gee_Img.py
"""
introduction: create virtual mosaic for images downloaded from GEE

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 July, 2025
"""

import sys,os
from optparse import OptionParser

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic

from delete_bad_images import get_image_valid_percent_entropy, get_img_grid_id_from_path

def get_image_list_with_max_entropy(image_list, valid_percent_list, img_entropy_list):

    # as we already set valid pixel percent threshold during download images,
    # so only use the entropy, keep the image with maximum entropy if multiple one exist
    gird_id_img_path = {}
    grid_id_max_entropy = {}
    for img, v_per, entropy in zip(image_list, valid_percent_list, img_entropy_list):
        grid_id = get_img_grid_id_from_path(img)
        if grid_id in gird_id_img_path.keys():
            if grid_id_max_entropy[grid_id] < entropy:
                gird_id_img_path[grid_id] = img
                grid_id_max_entropy[grid_id] = entropy
        else:
            gird_id_img_path[grid_id] = img
            grid_id_max_entropy[grid_id] = entropy

    return list(gird_id_img_path.values())


def create_virtual_mosaic(image_list, save_path, b_overview=True):
    # dir=ARTS-v2_1_0_4bands_S2_SR_HARMONIZED_20240701_2024830_images
    # find ${dir} | grep .tif > infile_list.txt
    # save_vrt=ARTS-v2_1_0_4bands_S2_SR_HARMONIZED_20240701_2024830_images.vrt
    # gdalbuildvrt -resolution average -r nearest -input_file_list infile_list.txt ${save_vrt}
    # #rm infile_list.txt

    img_list_txt = 'mosaic_infile_list.txt'
    io_function.save_list_to_txt(img_list_txt,image_list)

    cmd_str = f"gdalbuildvrt -resolution average -r bilinear -input_file_list {img_list_txt} {save_path}"
    basic.os_system_exit_code(cmd_str)

    # generated overview
    if b_overview:
        cmd_str = f"gdaladdo  -ro {save_path} 4 8 16 32 64 128"
        basic.os_system_exit_code(cmd_str)


def main(options, args):

    file_pattern = options.file_pattern
    image_dir = args[0]

    image_list = io_function.get_file_list_by_pattern(image_dir, file_pattern)
    valid_percent_list, img_entropy_list = get_image_valid_percent_entropy(image_list)

    img_list_max_entropy = get_image_list_with_max_entropy(image_list,valid_percent_list,img_entropy_list)

    vrt_mosaic_path = image_dir + '.vrt'
    create_virtual_mosaic(img_list_max_entropy,vrt_mosaic_path, b_overview=True)



if __name__ == '__main__':
    usage = "usage: %prog [options] image_folder "
    parser = OptionParser(usage=usage, version="1.0 2025-07-18")
    parser.description = 'Introduction: create mosaic for images downloaded from Google Earth Engine '


    parser.add_option("", "--file_pattern",
                      action="store", dest="file_pattern",default='*.tif',
                      help="the pattern to get raster list in a folder")

    parser.add_option("", "--nodata",
                      action="store", dest="nodata",type=int,default='0',
                      help="the user defined Nodata")

    # parser.add_option("-i", "--id_column",
    #                   action="store", dest="id_column",default='id',
    #                   help="the name of unique ID column")


    (options, args) = parser.parse_args()
    # print(options.planet_geojson)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

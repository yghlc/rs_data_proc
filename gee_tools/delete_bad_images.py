#!/usr/bin/env python
# Filename: delete_bad_images.py 
"""
introduction: check the quality of downloaded images from GEE.

1. remove valid pixel percent less than 80.
2. remove images with low entropy (less than 0.5?), we didn't set a threshold when downloading GEE
 because it may remove snow (or water etc) images
3. remove *.submit files is not images left after deleting for an image ID.
4. increase "max_count" then run the downloaing processing again, hope to download more at least one image.



authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 16 July, 2025
"""
import sys,os
from optparse import OptionParser
import time
import re
from tqdm import tqdm

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import vector_gpd
import raster_io
import basic_src.io_function as io_function


def get_img_grid_id_from_path(item):
    # return str
    return re.findall('img\d+', os.path.basename(item))[0][3:]


def get_image_valid_percent_entropy(image_list, nodata_user=0):
    valid_percent_list = []
    img_entropy_list = []
    for idx, img in enumerate(tqdm(image_list, desc="Calculating image valid and entropy")):
        valid_per, entropy = raster_io.get_valid_percent_shannon_entropy(img,nodata_input=nodata_user)
        valid_percent_list.append(valid_per)
        img_entropy_list.append(entropy)
    return valid_percent_list, img_entropy_list


def get_image_lower_than_thresholds(image_list, valid_percent, entropy_list, valid_threshold=80, entropy_thr=0.5):
    save_image_values_dict = {}
    for img, valid_p, entropy_v in zip(image_list, valid_percent, entropy_list):
        if valid_p < valid_threshold or entropy_v < entropy_thr:
            # print(f"Image: {img}, Value: {value:.2f} (below threshold: {threshold})")
            save_image_values_dict[img] = [valid_p, entropy_v]

    io_function.save_dict_to_txt_json(f'images_below_{valid_threshold}_or_{entropy_thr}.json', save_image_values_dict)

    return  list(save_image_values_dict.keys())

def delete_bad_images(image_list, bak_dir='./deleted_files'):
    if os.path.isdir(bak_dir):
        io_function.mkdir(bak_dir)
    for img in image_list:
        io_function.movefiletodir(img,bak_dir,overwrite=True)

def get_ids_without_images(all_image_id, org_image_list, delete_img_list):

    org_set = set(org_image_list)
    remove_set = set(delete_img_list)
    remain_image_list = list(org_set - remove_set)

    remain_image_ids = [get_img_grid_id_from_path(item) for item in remain_image_list]
    id_without_images = list(set(all_image_id) - set(remain_image_ids))

    io_function.save_list_to_txt('ids_without_download_image.txt',id_without_images)

    return id_without_images


def delete_submit_file_if_no_image_for_the_id(id_without_images, image_dir, bak_dir='./deleted_files'):
    if os.path.isdir(bak_dir):
        io_function.mkdir(bak_dir)

    for id_no_img in id_without_images:
        submit_file = os.path.join(image_dir,f'img{id_no_img}.submit')
        if os.path.isfile(submit_file):
            io_function.movefiletodir(submit_file,bak_dir,overwrite=True)

def main(options, args):

    extent_shp = args[0]
    id_column_name = options.id_column
    file_pattern = options.file_pattern
    image_dir = args[1]
    back_up_dir = './deleted_files'

    extent_polygons = vector_gpd.read_shape_gpd_to_NewPrj(extent_shp, 'EPSG:4326')
    extent_ids = vector_gpd.read_attribute_values_list(extent_shp, id_column_name)

    if extent_ids is None:
        extent_ids = [ str(item) for item in range(len(extent_polygons))]
    else:
        extent_ids = [item for item in extent_ids]

    image_list = io_function.get_file_list_by_pattern(image_dir,file_pattern)

    valid_percent_list, img_entropy_list = get_image_valid_percent_entropy(image_list)

    deleted_images = get_image_lower_than_thresholds(image_list, valid_percent_list, img_entropy_list,
                                    valid_threshold=80, entropy_thr=0.5)


    id_without_images = get_ids_without_images(extent_ids,image_list,deleted_images)

    delete_bad_images(deleted_images,bak_dir=back_up_dir)
    delete_submit_file_if_no_image_for_the_id(id_without_images, image_dir, bak_dir=back_up_dir)




if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp image_folder "
    parser = OptionParser(usage=usage, version="1.0 2025-07-16")
    parser.description = 'Introduction: delete bad images download from Google Earth Engine '


    parser.add_option("", "--file_pattern",
                      action="store", dest="file_pattern",default='*.tif',
                      help="the pattern to get raster list in a folder")

    parser.add_option("", "--nodata",
                      action="store", dest="nodata",type=int,default='0',
                      help="the user defined Nodata")

    parser.add_option("-i", "--id_column",
                      action="store", dest="id_column",default='id',
                      help="the name of unique ID column")


    (options, args) = parser.parse_args()
    # print(options.planet_geojson)

    if len(sys.argv) < 2 or len(args) < 2:
        parser.print_help()
        sys.exit(2)

    main(options, args)


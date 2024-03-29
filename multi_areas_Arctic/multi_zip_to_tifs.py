#!/usr/bin/env python
# Filename: multi_zip_to_tifs.py
"""
introduction: unzip and organize Rapid and PlanetScope imagery

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 24 June, 2023
"""

import os, sys
import time
from optparse import OptionParser
from datetime import datetime

deeplabforRS = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io

data_dir = os.path.expanduser('~/Data/Arctic/canada_arctic')
work_dir = os.getcwd()
# regions = ['Willow_River', 'Banks_east', 'Ellesmere_Island']
region_dirs = ['Willow_River','Banks_east','Ellesmere_Island','TukCoast']

py8bit = os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py')

def is_file_exist(file_path):
    if os.path.isfile(file_path):
        print('warning, %s exists, skip' % file_path)
        return True
    return False

def create_tif_mosaic(tif_list, save_path):
    if is_file_exist(save_path):
        return save_path
    cmd_str = 'gdal_merge.py -o %s '%save_path
    nodata = raster_io.get_nodata(tif_list[0])
    if nodata is not None:
        cmd_str += '-n %s -a_nodata %s'%(str(nodata), str(nodata))

    for tif in tif_list:
        cmd_str += ' %s'%tif
    basic.os_system_exit_code(cmd_str)

def convert_to_8bit(input_tif,save_dir):
    save_path = os.path.join(save_dir, os.path.basename(io_function.get_name_by_adding_tail(input_tif,'8bit')))
    if is_file_exist(save_path):
        return save_path
    cmd_str = '%s -u 0.99 -l 0.01 %s  %s' % (py8bit,input_tif,save_path)
    basic.os_system_exit_code(cmd_str)
    return save_path


def extract_three_bands(input_tif, save_dir, bands, name_tail='rgb'):
    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)
    save_path = os.path.join(save_dir, os.path.basename(io_function.get_name_by_adding_tail(input_tif, name_tail)))
    if is_file_exist(save_path):
        return save_path
    cmd_str = 'gdal_translate -b %d -b %d -b %d '%(bands[0], bands[1], bands[2])
    cmd_str += '-of GTiff -a_nodata 0 -co compress=lzw -co tiled=yes -co bigtiff=if_safer '
    cmd_str += ' %s %s'%(input_tif, save_path)
    basic.os_system_exit_code(cmd_str)

    raster_io.remove_nodata_from_raster_metadata(save_path)

def extract_PlanetScope_rgb_bands(input_tif,save_dir):
    # https://developers.planet.com/docs/apis/data/sensors/
    extract_three_bands(input_tif,save_dir,[3,2,1],name_tail='rgb') # PlanetScope: RGB is band 3, 2, 1

def extract_PlanetScope_nirGB_bands(input_tif, save_dir):
    extract_three_bands(input_tif, save_dir, [4, 2, 1],name_tail='nirGB')  # PlanetScope: Nir GB is band 4, 2, 1

def extract_RapidEye_rgb_bands(input_tif,save_dir):
    # https://earth.esa.int/eogateway/missions/rapideye
    extract_three_bands(input_tif,save_dir,[3,2,1],name_tail='rgb') # RapidEye: RGB is band 3, 2, 1

def extract_RapidEye_nirGB_bands(input_tif,save_dir):
    extract_three_bands(input_tif,save_dir,[5,2,1],name_tail='nirGB') # RapidEye: Nir GB is band 3, 2, 1

def one_zip_to_images(zip_path, save_dir):

    zip_folder = io_function.unzip_file(zip_path,work_dir)
    zip_name_lower = os.path.basename(zip_path).lower()

    composite_tif = os.path.join(zip_folder,'composite.tif')
    mosaic_tif = os.path.join(zip_folder, os.path.basename(zip_folder) + '_composite.tif')
    if os.path.isfile(composite_tif):
        # rename
        io_function.move_file_to_dst(composite_tif,mosaic_tif)
    elif os.path.isfile(mosaic_tif):
        pass
    else:
        # get a mosaic
        tif_list = io_function.get_file_list_by_pattern(zip_folder,'*/*SR*.tif')
        if len(tif_list) < 1:
            raise IOError('Not tif in %s'%zip_folder)
        mosaic_tif = os.path.join(zip_folder, os.path.basename(zip_folder) + '_mosaic.tif')
        create_tif_mosaic(tif_list,mosaic_tif)

    # co-registration (sometime, there is an offset between RapidEye and PlanetScope images)
    # I built a software (ImageMatchsiftGPU) based on SiftGPU, a legacy software depending on openGL or CUDA,
    # I built a virtual machine (Ubuntu 14.04) to run the co-registraion,
    # need to copy data to the Virtual machine, then run the program, then copy results back,
    # see ../coregistration_siftGPU/

    # use the mosaic after co-registration if it exists
    mosaic_tif_coreg = io_function.get_name_by_adding_tail(mosaic_tif,'new_warp')
    if os.path.isfile(mosaic_tif_coreg):
        mosaic_tif = mosaic_tif_coreg

    # convert to 8bit
    mosaic_tif_8bit = convert_to_8bit(mosaic_tif,zip_folder)

    # to RGB and nirGB
    f_base_name = os.path.basename(save_dir)
    save_rgb_dir = os.path.join(save_dir,f_base_name+'_8bit_RGB')
    save_nirGB_dir = os.path.join(save_dir,f_base_name+'_8bit_nirRGB')
    if 'rapideye'.lower() in zip_name_lower:
        extract_RapidEye_rgb_bands(mosaic_tif_8bit,save_rgb_dir)
        extract_RapidEye_nirGB_bands(mosaic_tif_8bit,save_nirGB_dir)
    elif 'planetscope'.lower() in zip_name_lower:
        extract_PlanetScope_rgb_bands(mosaic_tif_8bit,save_rgb_dir)
        extract_PlanetScope_nirGB_bands(mosaic_tif_8bit,save_nirGB_dir)
    else:
        raise ValueError('does not find RapidEye or PlanetScope in the zip file name: %s'%zip_path)

    return zip_folder

def main(options, args):
    global data_dir
    if len(args)==1:
        data_dir = args[0]
        # region_dirs = options.region_name
        zip_list = io_function.get_file_list_by_pattern(data_dir, '*.zip')
        if len(zip_list)< 1:
            print('No download zip files in %s/*download_images'%data_dir)
            return
        [print(item) for item in zip_list]
        out_folders = []
        for idx, a_zip in enumerate(zip_list):
            o_folder = one_zip_to_images(a_zip, data_dir)
            out_folders.append(o_folder)

    else:
        for reg in region_dirs:
            print(datetime.now(), 'working on the region: %s'%reg)
            reg_dir = os.path.join(data_dir, reg)
            zip_list = io_function.get_file_list_by_pattern(reg_dir,'*download_images/*.zip')
            if len(zip_list)< 1:
                print('No download zip files in %s/*download_images'%reg_dir)
                continue
            else:
                print('zip (%d) files:'%len(zip_list))
            [print(item) for item in zip_list]
            # print(zip_list)
            out_folders = []
            for idx, a_zip in enumerate(zip_list):
                o_folder = one_zip_to_images(a_zip, reg_dir)
                out_folders.append(o_folder)
                # for test
                # if idx > 3:
                #     break

        # copy and process 2020 images


        # # remove and clear
        # for out_dir in out_folders:
        #     print('removing %s'%out_dir)
        #     io_function.delete_file_or_dir(out_dir)


if __name__ == '__main__':
    usage = "usage: %prog [options] data_dir "
    parser = OptionParser(usage=usage, version="1.0 2023-6-24")
    parser.description = 'Introduction: unzip and organize Rapid and PlanetScope imagery  '

    parser.add_option('-n',"--region_name",
                      action="store", dest="region_name",
                      help="region name")

    # parser.add_option("-p", "--process_num",
    #                   action="store", dest="process_num", type=int, default=8,
    #                   help="number of processes to checking invalid tifs")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    # if len(sys.argv) < 2 or len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)



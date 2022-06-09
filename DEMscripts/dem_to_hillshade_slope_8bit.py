#!/usr/bin/env python
# Filename: dem_to_hillshade_slope_8bit 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 20 March, 2021
"""
import os,sys
from optparse import OptionParser
import time
from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic

machine_name = os.uname()[1]

# some folder paths
from dem_common import ArcticDEM_tmp_dir
from dem_common import check_create_lock, release_lock

py8bit= os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py')

from multiprocessing import Pool

def slope_to_8bit(input, output):

    dst_nodat = 255
    hist_max_percent=0.98
    hist_min_percent=0.02
    MIN_MAX_VALUE = '0 70'   # slop range from 0 to 70

    command_str = py8bit + ' ' + input + ' ' + output
    command_str += ' -n ' + str(dst_nodat)
    command_str += ' -u ' + str(hist_max_percent) + ' -l ' + str(hist_min_percent)
    command_str += ' -m ' + MIN_MAX_VALUE

    # print(command_str)
    basic.os_system_exit_code(command_str)
    return True

def dem_to_slope(input,output,slope_file_bak):

    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return output
    if os.path.isfile(slope_file_bak):
        basic.outputlogMessage('%s exists, skip' % slope_file_bak)
        return slope_file_bak

    if os.path.isfile(input) is False:
        basic.outputlogMessage('Waring, %s does not exist'%input)
        return False

    # # use the default setting in QGIS
    command_str = 'gdaldem slope %s %s -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -b 1 -s 1.0'%(input,output)
    basic.os_system_exit_code(command_str)
    return output

def tpi_to_8bit(input,output):
    dst_nodat = 255
    hist_max_percent=0.98
    hist_min_percent=0.02
    MIN_MAX_VALUE = '-1 1'   # tpi range from -1 to 1

    command_str = py8bit + ' ' + input + ' ' + output
    command_str += ' -n ' + str(dst_nodat)
    command_str += ' -u ' + str(hist_max_percent) + ' -l ' + str(hist_min_percent)
    command_str += ' -m ' + MIN_MAX_VALUE

    # print(command_str)
    basic.os_system_exit_code(command_str)
    return True

def dem_to_tpi_save_8bit(input,output):
    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    if os.path.isfile(input) is False:
        basic.outputlogMessage('Waring, %s does not exist'%input)
        return False

    # Topographic Position Index
    tpi_file = os.path.basename(io_function.get_name_by_adding_tail(input, 'tpi'))
    command_str = 'gdaldem TPI %s %s -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -b 1 '%(input,tpi_file)
    basic.os_system_exit_code(command_str)

    # to 8bit
    if tpi_to_8bit(tpi_file,output) is True:
        io_function.delete_file_or_dir(tpi_file)
    return True


def dem_to_hillshade(input,output):

    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True
    if os.path.isfile(input) is False:
        basic.outputlogMessage('Waring, %s does not exist'%input)
        return False

    # use the default setting in QGIS
    # gdaldem hillshade ${dem} ${hillshade} -of GTiff -b 1 -z 1.0 -s 1.0 -az 315.0 -alt 45.0
    command_str = 'gdaldem hillshade %s  %s -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -b 1 -z 1.0 -s 1.0 -az 315.0 -alt 45.0'%(input,output)
    basic.os_system_exit_code(command_str)
    return True

def process_one_dem(idx, count, tif,product_list, arcticDEM_slope_dir,arcticDEM_slope_8bit_dir, arcticDEM_hillshade_dir,arcticDEM_tpi_8bit_dir):
    print('%d/%d convert %s to slope (8bit) and hillshade' % (idx + 1, count, tif))

    try:
        slope_file = os.path.basename(io_function.get_name_by_adding_tail(tif, 'slope'))
        slope_file_bak = os.path.join(arcticDEM_slope_dir, os.path.basename(slope_file))
        if 'slope' in product_list or 'slope_8bit' in product_list:
            slope_out = dem_to_slope(tif,slope_file,slope_file_bak)
            if slope_out is False:
                pass
            else:
                if 'slope_8bit' in product_list:
                    slope_8bit = io_function.get_name_by_adding_tail(tif, 'slope8bit')
                    slope_8bit = os.path.join(arcticDEM_slope_8bit_dir, os.path.basename(slope_8bit))
                    slope_to_8bit(slope_file, slope_8bit)

                # delete or move the slope file
                if 'slope' in product_list:
                    io_function.move_file_to_dst(slope_file,slope_file_bak)
                else:
                    io_function.delete_file_or_dir(slope_file)


        if 'hillshade' in product_list:
            hillshapde = io_function.get_name_by_adding_tail(tif, 'hillshade')
            hillshapde = os.path.join(arcticDEM_hillshade_dir, os.path.basename(hillshapde))
            dem_to_hillshade(tif, hillshapde)

        if 'tpi' in product_list:
            tip_8bit = io_function.get_name_by_adding_tail(tif, 'TPI8bit')
            tip_8bit = os.path.join(arcticDEM_tpi_8bit_dir, os.path.basename(tip_8bit))
            dem_to_tpi_save_8bit(tif, tip_8bit)

        return True
    except:
        print('failed in process %s'%tif)
        return tif

def main(options, args):
    b_mosaic_ArcticDEM = options.b_mosaic_ArcticDEM
    process_num = options.process_num
    product_list = args # subset of [slope, slope_8bit, hillshade, tpi]
    print('Will produce products includes: %s'%(product_list))

    if b_mosaic_ArcticDEM:
        print('Input is the mosaic version of AricticDEM')
        from dem_common import arcticDEM_tile_reg_tif_dir,arcticDEM_tile_hillshade_dir,arcticDEM_tile_slope_8bit_dir,\
            arcticDEM_tile_slope_dir,arcticDEM_tile_tpi_8bit_dir
        arcticDEM_reg_tif_dir = arcticDEM_tile_reg_tif_dir
        arcticDEM_hillshade_dir = arcticDEM_tile_hillshade_dir
        arcticDEM_slope_8bit_dir = arcticDEM_tile_slope_8bit_dir
        arcticDEM_slope_dir = arcticDEM_tile_slope_dir
        arcticDEM_tpi_8bit_dir = arcticDEM_tile_tpi_8bit_dir
        dem_pattern = '*reg_dem.tif'
    else:
        from dem_common import arcticDEM_reg_tif_dir,dem_hillshade_dir,dem_slope_dir,dem_slope_8bit_dir,dem_tpi_8bit_dir
        arcticDEM_reg_tif_dir = arcticDEM_reg_tif_dir
        arcticDEM_hillshade_dir = dem_hillshade_dir
        arcticDEM_slope_8bit_dir = dem_slope_8bit_dir
        arcticDEM_slope_dir = dem_slope_dir
        arcticDEM_tpi_8bit_dir = dem_tpi_8bit_dir
        dem_pattern = '*dem_reg.tif'

    basic.setlogfile('log_dem_to_slope8bit_hillshade.txt')

    if os.path.isdir(arcticDEM_slope_8bit_dir) is False:
        io_function.mkdir(arcticDEM_slope_8bit_dir)
    if os.path.isdir(arcticDEM_hillshade_dir) is False:
        io_function.mkdir(arcticDEM_hillshade_dir)
    if os.path.isdir(arcticDEM_tpi_8bit_dir) is False:
        io_function.mkdir(arcticDEM_tpi_8bit_dir)
    if os.path.isdir(arcticDEM_slope_dir) is False:
        io_function.mkdir(arcticDEM_slope_dir)

    failed_tifs = []

    # create a lock file (make sure only one workstation is working on producing slope)
    arcticDEM_slope_lock = os.path.join(arcticDEM_slope_dir,'arcticDEM_slope_lock.txt')
    check_create_lock(arcticDEM_slope_lock,'because other program is creating slope files into %s'%arcticDEM_slope_dir)

    dem_reg_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,dem_pattern)
    count = len(dem_reg_list)
    print('Find %d DEM in %s, with pattern: %s'%(count,arcticDEM_reg_tif_dir,dem_pattern))
    if process_num == 1:
        for idx, tif in enumerate(dem_reg_list):
            res = process_one_dem(idx, count, tif, product_list, arcticDEM_slope_dir,arcticDEM_slope_8bit_dir, arcticDEM_hillshade_dir,arcticDEM_tpi_8bit_dir)
            if res is not True:
                failed_tifs.append(res)
    else:
        parameters_list = [(idx, count, tif, product_list, arcticDEM_slope_dir,arcticDEM_slope_8bit_dir, arcticDEM_hillshade_dir,arcticDEM_tpi_8bit_dir) for
                           idx, tif in enumerate(dem_reg_list)]
        theadPool = Pool(process_num)
        results = theadPool.starmap(process_one_dem, parameters_list)
        for res in results:
            if res is not True:
                failed_tifs.append(res)
        theadPool.close()

    with open('to_hillshade_slope8bit_failed_cases.txt','w') as f_obj:
        for item in failed_tifs:
            f_obj.writelines(item + '\n')

    # delete the lock file
    release_lock(arcticDEM_slope_lock)

if __name__ == '__main__':
    usage = "usage: %prog [options] products (slope, slope_8bit, hillshade, tpi) "
    parser = OptionParser(usage=usage, version="1.0 2021-3-20")
    parser.description = 'Introduction: producing DEM hillshade, slope, and TPI (Topographic Position Index), and save to 8bit '


    parser.add_option("-m", "--b_mosaic_ArcticDEM",
                      action="store_true", dest="b_mosaic_ArcticDEM",default=False,
                      help="whether indicate the input is ArcticDEM mosaic version")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")


    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
    pass
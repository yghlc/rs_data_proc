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

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic

machine_name = os.uname()[1]

# some folder paths
if machine_name == 'uist':
    ArcticDEM_tmp_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
elif machine_name == 'ubuntu':  # tesia
    ArcticDEM_tmp_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir'
else:
    ArcticDEM_tmp_dir = './'

py8bit= os.path.expanduser('~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py')

arcticDEM_reg_tif_dir = os.path.join(ArcticDEM_tmp_dir,'registration_tifs')
arcticDEM_hillshade_dir = os.path.join(ArcticDEM_tmp_dir,'dem_hillshade')
arcticDEM_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_slope_8bit')


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

def dem_to_slope_save_8bit(input,output):

    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True
    slope_file = os.path.basename(io_function.get_name_by_adding_tail(input,'slope'))

    # # use the default setting in QGIS
    command_str = 'gdaldem slope %s %s -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -b 1 -s 1.0'%(input,slope_file)
    basic.os_system_exit_code(command_str)

    # to 8bit
    if slope_to_8bit(slope_file,output) is True:
        io_function.delete_file_or_dir(slope_file)
    return True


def dem_to_hillshade(input,output):

    if os.path.isfile(output):
        basic.outputlogMessage('%s exists, skip'%output)
        return True

    # use the default setting in QGIS
    # gdaldem hillshade ${dem} ${hillshade} -of GTiff -b 1 -z 1.0 -s 1.0 -az 315.0 -alt 45.0
    command_str = 'gdaldem hillshade %s  %s -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer -b 1 -z 1.0 -s 1.0 -az 315.0 -alt 45.0'%(input,output)
    basic.os_system_exit_code(command_str)
    return True


def main():

    basic.setlogfile('log_dem_to_slope8bit_hillshade.txt')

    if os.path.isdir(arcticDEM_slope_8bit_dir) is False:
        io_function.mkdir(arcticDEM_slope_8bit_dir)
    if os.path.isdir(arcticDEM_hillshade_dir) is False:
        io_function.mkdir(arcticDEM_hillshade_dir)

    dem_reg_list = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir,'*dem_reg.tif')
    count = len(dem_reg_list)
    for idx, tif in enumerate(dem_reg_list):
        print('%d/%d convert %s to slope (8bit) and hillshade'%(idx+1, count, tif))

        slope_8bit = io_function.get_name_by_adding_tail(tif,'slope8bit')
        slope_8bit = os.path.join(arcticDEM_slope_8bit_dir, os.path.basename(slope_8bit))
        dem_to_slope_save_8bit(tif, slope_8bit)

        hillshapde = io_function.get_name_by_adding_tail(tif,'hillshade')
        hillshapde = os.path.join(arcticDEM_hillshade_dir, os.path.basename(hillshapde))
        dem_to_hillshade(tif,hillshapde)

if __name__ == '__main__':
    main()
    pass
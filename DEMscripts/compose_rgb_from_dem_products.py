#!/usr/bin/env python
# Filename: dem_to_hillshade_slope_8bit 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 23 April, 2021
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

from multiprocessing import Pool

def find_only_one_file(folder, pattern):
    tif_list = io_function.get_file_list_by_pattern(folder, pattern)
    if len(tif_list) == 1:
        return tif_list[0]
    return False

def combine_hillshade_slope_tip(hillshade, slope, tpi,save_path):
    if os.path.isfile(save_path):
        print('%s exists, skip'%save_path)
        return

    # combine them.
    command_str = 'gdal_merge.py -o %s -separate -a_nodata 255 -of GTiff -co compress=lzw -co tiled=yes -co bigtiff=if_safer %s %s %s'%(save_path, hillshade, slope, tpi)
    basic.os_system_exit_code(command_str)


def compose_rgb_from_hillshade_slope_tpi(idx, count, hillshade, arcticDEM_slope_8bit_dir, arcticDEM_tpi_8bit_dir,arcticDEM_topo_rgb_dir):
    print('%d/%d get RGB from hillshade (%s), slope, and tpi' % (idx + 1, count, hillshade))

    try:
        base_name = os.path.basename(hillshade)[:-14] # remove '_hillshade.tif'
        slope_file = find_only_one_file(arcticDEM_slope_8bit_dir, base_name + '*.tif')
        if slope_file is False:
            print('Cannot find corresponding slope file')
            return hillshade
        tip_file = find_only_one_file(arcticDEM_tpi_8bit_dir, base_name + '*.tif')
        if tip_file is False:
            print('Cannot find corresponding TPI file')
            return hillshade

        save_path = os.path.join(arcticDEM_topo_rgb_dir, base_name + '_topoRGB.tif')
        combine_hillshade_slope_tip(hillshade,slope_file,tip_file,save_path)

        return True
    except:
        return hillshade

def main(options, args):
    b_mosaic_ArcticDEM = options.b_mosaic_ArcticDEM
    process_num = options.process_num

    if b_mosaic_ArcticDEM:
        print('Input is the mosaic version of AricticDEM')
        arcticDEM_hillshade_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_hillshade')
        arcticDEM_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'arcticdem_mosaic_slope_8bit')
        arcticDEM_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir,'arcticdem_mosaic_tpi_8bit')
        file_pattern = '*.tif'
        arcticDEM_topo_rgb_dir = os.path.join(ArcticDEM_tmp_dir,'arcticdem_mosaic_topoRGB_8bit')
    else:
        arcticDEM_hillshade_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_hillshade')
        arcticDEM_slope_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_slope_8bit')
        arcticDEM_tpi_8bit_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_tpi_8bit')
        file_pattern = '*.tif'
        arcticDEM_topo_rgb_dir = os.path.join(ArcticDEM_tmp_dir, 'dem_topoRGB_8bit')

    basic.setlogfile('log_compose_dem_products.txt')

    if os.path.isdir(arcticDEM_topo_rgb_dir) is False:
        io_function.mkdir(arcticDEM_topo_rgb_dir)

    failed_tifs = []

    hillshade_list = io_function.get_file_list_by_pattern(arcticDEM_hillshade_dir,file_pattern)
    count = len(hillshade_list)
    print('Find %d Hillshade in %s, with pattern: %s'%(count,arcticDEM_hillshade_dir,file_pattern))
    if process_num == 1:
        for idx, tif in enumerate(hillshade_list):
            res = compose_rgb_from_hillshade_slope_tpi(idx, count, tif, arcticDEM_slope_8bit_dir,
                                                       arcticDEM_tpi_8bit_dir,arcticDEM_topo_rgb_dir)
            if res is not True:
                failed_tifs.append(res)
    else:
        parameters_list = [(idx, count, tif, arcticDEM_slope_8bit_dir,arcticDEM_tpi_8bit_dir,arcticDEM_topo_rgb_dir) for
                           idx, tif in enumerate(hillshade_list)]
        theadPool = Pool(process_num)
        results = theadPool.starmap(compose_rgb_from_hillshade_slope_tpi, parameters_list)
        for res in results:
            if res is not True:
                failed_tifs.append(res)
        theadPool.close()

    with open('to_rgb_from_hillshade_slope8bit_failed_cases.txt','w') as f_obj:
        for item in failed_tifs:
            f_obj.writelines(item + '\n')

if __name__ == '__main__':
    usage = "usage: %prog [options]  "
    parser = OptionParser(usage=usage, version="1.0 2021-4-23")
    parser.description = 'Introduction: '


    parser.add_option("-m", "--b_mosaic_ArcticDEM",
                      action="store_true", dest="b_mosaic_ArcticDEM",default=False,
                      help="whether indicate the input is ArcticDEM mosaic version")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")


    (options, args) = parser.parse_args()

    # if len(sys.argv) < 2 or len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)
    pass
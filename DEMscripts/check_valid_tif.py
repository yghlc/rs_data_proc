#!/usr/bin/env python
# Filename: check_valid_tif.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 09 March, 2021
"""

import os,sys
from optparse import OptionParser
import time

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic
import raster_io

machine_name = os.uname()[1]

# some folder paths
from dem_common import arcticDEM_reg_tif_dir, grid_dem_diffs_dir

import multiprocessing
from multiprocessing import Pool

def check_one_tif(idx,total,tif_path, good_tif_list):
    if os.path.basename(tif_path) in good_tif_list:
        return True
    try:
        print('checking %d/%d' % (idx, total))
        # src = raster_io.open_raster_read(tif_path)
        data, nodata = raster_io.read_raster_all_bands_np(tif_path)
    except:
        print('invalid tif: %s' % tif_path)
        return False

    return True


def main(options, args):

    # process_num = multiprocessing.cpu_count()
    process_num = options.process_num

    tifs = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir, '*.tif') # _dem_reg check all tifs
    save_invalid_txt_path = os.path.basename(arcticDEM_reg_tif_dir) + '_invalid_list.txt'
    save_good_txt_path = os.path.basename(arcticDEM_reg_tif_dir) + '_good_list.txt'
    tif_count = len(tifs)

    good_tifs = []
    if os.path.isfile(save_good_txt_path):
        good_tifs.extend(io_function.read_list_from_txt(save_good_txt_path))
    invalid_tif = []

    # remove good one for the list
    if len(good_tifs)>0:
        tifs = [item for item in tifs if os.path.basename(item) not in good_tifs]

    if process_num == 1:
        # tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir, bsub_folder=False)
        for idx,tif in enumerate(tifs):
            if check_one_tif(idx,tif_count,tif,good_tifs):
                good_tifs.append(os.path.basename(tif))
            else:
                invalid_tif.append(os.path.basename(tif))
    else:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(idx,tif_count,tif,good_tifs) for idx,tif in enumerate(tifs)]
        results = theadPool.starmap(check_one_tif, parameters_list)  # need python3
        for tif, res in zip(tifs,results):
            if res:
                good_tifs.append(os.path.basename(tif))
            else:
                invalid_tif.append(os.path.basename(tif))

    io_function.save_list_to_txt(save_invalid_txt_path, invalid_tif)
    io_function.save_list_to_txt(save_good_txt_path,good_tifs)
        

    # tifs = io_function.get_file_list_by_ext('.tif', grid_dem_diff_dir, bsub_folder=False)
    # invalid_tif = check_tif(tifs)
    # for tif in invalid_tif:
    #     print('removing %s'%tif)
        # io_function.delete_file_or_dir(tif)


if __name__ == '__main__':
    usage = "usage: %prog [options]  "
    parser = OptionParser(usage=usage, version="1.0 2021-9-6")
    parser.description = 'Introduction: find out invalid tifs  '

    parser.add_option("-p", "--process_num",
                      action="store", dest="process_num", type=int, default=8,
                      help="number of processes to checking invalid tifs")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    # if len(sys.argv) < 2 or len(args) < 1:
    #     parser.print_help()
    #     sys.exit(2)

    main(options, args)
    pass
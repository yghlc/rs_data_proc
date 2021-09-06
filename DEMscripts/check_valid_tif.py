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

def check_one_tif(idx,total,tif_path):
    try:
        print('checking %d/%d' % (idx, total))
        # src = raster_io.open_raster_read(tif_path)
        data, nodata = raster_io.read_raster_all_bands_np(tif_path)
    except:
        print('invalid tif: %s' % tif_path)
        return tif_path

    return None


def main(options, args):

    # process_num = multiprocessing.cpu_count()
    process_num = options.process_num

    tifs = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir, '*_dem_reg.tif')
    save_invalid_txt_path = os.path.basename(arcticDEM_reg_tif_dir) + '_invalid_list.txt'
    # save_good_txt_path = os.path.basename(arcticDEM_reg_tif_dir) + '_good_list.txt'
    tif_count = len(tifs)

    invalid_tif = []

    if process_num == 1:
        # tifs = io_function.get_file_list_by_ext('.tif',arcticDEM_reg_tif_dir, bsub_folder=False)
        for idx,tif in enumerate(tifs):
            tif_name = os.path.basename(tif)
            result = check_one_tif(idx,tif_count,tif)
            if result is None:
                pass
            else:
                invalid_tif.append(tif_name)
    else:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(idx,tif_count,tif) for idx,tif in enumerate(tifs)]
        results = theadPool.starmap(check_one_tif, parameters_list)  # need python3
        invalid_tif = [ os.path.basename(out) for out in results if out is not None]

    io_function.save_list_to_txt(save_invalid_txt_path, invalid_tif)
        

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
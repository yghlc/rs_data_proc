#!/usr/bin/env python
# Filename: dem_coregistration 
"""
introduction: dem co-registration

using demcoreg, install it and its dependencies:

    git clone https://github.com/dshean/imview.git
    pip install -e imview

   git clone https://github.com/dshean/pygeotools.git
   pip install -e pygeotools

   git clone https://github.com/dshean/demcoreg.git
   pip install -e demcoreg

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 March, 2021
"""

import os, sys
from optparse import OptionParser
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function


def choose_reference_dem(dem_list, dem_valid_per_txt):
    if dem_valid_per_txt is None:
        raise ValueError('NO information of valid percentage of DEMs, cannot choose a reference DEM')

    with open(dem_valid_per_txt, 'r') as f_obj:
        tif_valid_per_list = [line.strip().split() for line in f_obj.readlines()]
        tif_higest_per = None
        per_max = 0
        for tif, per in tif_valid_per_list:
            if float(per) > per_max:
                per_max = float(per)
                tif_higest_per = tif

        for dem_tif in dem_list:
            if tif_higest_per in dem_tif:
                return dem_tif

    return None

def co_registration_parallel(ref_dem, dem_list, save_dir, process_num):

    # parallel --progress --delay 10 -j 14 dem_align.py ${ref} {} ::: $(ls *_dem.tif | grep -v 2012)
    commond_str = 'parallel --progress --delay 5 -j %d %s'%(process_num, ref_dem)
    commond_str += ' \{} ::: ' + ' '.join(dem_list)
    os.system(commond_str)

def main(options, args):

    save_dir = options.save_dir
    dem_dir_or_txt = args[0]
    ref_dem = options.ref_dem
    dem_valid_per_txt = options.dem_valid_per_txt
    process_num = options.process_num

    if os.path.isfile(dem_dir_or_txt):
        dem_list = io_function.read_list_from_txt(dem_dir_or_txt)
    else:
        dem_list = io_function.get_file_list_by_ext('.tif', dem_dir_or_txt, bsub_folder=False)
        if dem_valid_per_txt is None:
            dem_valid_per_txt = os.path.join(dem_dir_or_txt,'dem_valid_percent.txt')
    dem_count = len(dem_list)
    if dem_count < 1:
        raise ValueError('No input dem files in %s' % dem_dir_or_txt)

    if ref_dem is None:
        ref_dem = choose_reference_dem(dem_list, dem_valid_per_txt)
        if ref_dem is None:
            raise ValueError('Cannot find a reference DEM')

    if ref_dem in dem_list:
        dem_list.remove(ref_dem)
    co_registration_parallel(ref_dem,dem_list,save_dir,process_num)


if __name__ == '__main__':
    usage = "usage: %prog [options] dem_tif_dir or dem_list_txt "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: co-registration for multi-temporal DEMs '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir", default='./',
                      help="the folder to save results")

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes")

    parser.add_option("-r", "--ref_dem",
                      action="store", dest="ref_dem",
                      help="the reference DEM, if not set, it will select with the most coverage")

    parser.add_option("-p", "--dem_valid_per_txt",
                      action="store", dest="dem_valid_per_txt",
                      help="a txt file storing the valid percentage of all the DEM")

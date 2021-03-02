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
import basic_src.basic as basic
import raster_io

from multiprocessing import Process

dem_dem_align = os.path.expanduser('~/codes/github_public_repositories/demcoreg/demcoreg/dem_align.py')

# use psutil to check available memory.
import psutil
import time

# check the resource used by python
# import resource
# resource.getrusage(resource.RUSAGE_SELF)

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

def check_coreg_results(dem_tif, save_dir):
    dem_align = os.path.join(save_dir, 'dem_coreg', os.path.basename(io_function.get_name_by_adding_tail(dem_tif, 'coreg')))
    if os.path.isfile(dem_align):
        return True
    return False

def check_align_folder(dem_tif):
    # by default, dem_align.py save the results to where dem_tif is
    res_dir = os.path.dirname(dem_tif)
    align_folder = os.path.splitext(os.path.basename(dem_tif))[0] + '_dem_align'
    align_dir = os.path.join(res_dir,align_folder)
    # after dem_align.py usually have 9 files
    align_outputs = io_function.get_file_list_by_pattern(align_dir,'*')
    # print(align_outputs)
    return align_outputs


def move_align_results(ref_dem, dem_tif, save_dir):

    coreg_save_dir = os.path.join(save_dir, 'dem_coreg')
    if os.path.isdir(coreg_save_dir) is False:
        io_function.mkdir(coreg_save_dir)

    align_outputs = check_align_folder(dem_tif)
    if len(align_outputs) < 9:
        raise ValueError('the output of dem_align.py is less than 9 files')

    dem_align = os.path.join(coreg_save_dir, os.path.basename(io_function.get_name_by_adding_tail(dem_tif, 'coreg')))
    # align DEM and a filt version, which one should I use? what filter they apply?
    # visually check one results (Banks east) , a the same location, align DEM and a filt one have exact values,
    # but the filt version have more nodata.  Let's use the filt version.
    # the nodata pixels usually are water pixels, but also some inside the thaw slumps
    align_filt = [out for out in align_outputs if out.endswith('align_filt.tif')][0]
    io_function.move_file_to_dst(align_filt,dem_align, overwrite=True)

    # copy reference dem if necessary
    ref_dem_copy = os.path.join(coreg_save_dir, os.path.basename(ref_dem))
    if os.path.isfile(ref_dem_copy) is False:
        io_function.copy_file_to_dst(ref_dem,ref_dem_copy)

    # move the elevation difference?
    ele_diff_folder = os.path.join(save_dir,'dem_diff_from_demcoreg')
    if os.path.isdir(ele_diff_folder) is False:
        io_function.mkdir(ele_diff_folder)
    dem_diff_filt = [out for out in align_outputs if out.endswith('align_diff_filt.tif')][0]
    io_function.movefiletodir(dem_diff_filt,ele_diff_folder, overwrite=True)

    return True


# def co_registration_parallel(ref_dem, dem_list, save_dir, process_num):
#     print('ref_dem', ref_dem)
#     print('source dem:')
#     for dem_tif in dem_list:
#         print(dem_tif)
#
#     # dem_align.py requires large memory, for example, a region of 50000 by 50000 pixels, may requires more than 110 GB memory.
#
#     # parallel --progress --delay 10 -j 14 dem_align.py ${ref} {} ::: $(ls *_dem.tif | grep -v 2012)
#     commond_str = 'parallel --progress --delay 5 -j %d %s %s'%(process_num, dem_dem_align ,ref_dem)
#     commond_str += ' {} ::: ' + ' '.join(dem_list)
#     print(commond_str)
#     os.system(commond_str)

    # move results to another folder

def co_registration_one_dem(ref_dem, dem_tif, save_dir):
    if check_coreg_results(dem_tif,save_dir):
        return 0

    align_outputs = check_align_folder(dem_tif)
    if len(align_outputs) >= 9:
        print('%s has been co-registered, skip'%dem_tif)
    else:
        commond_str = dem_dem_align + ' ' + ref_dem + ' ' + dem_tif
        print(commond_str)
        res = os.system(commond_str)
        if res != 0:
            sys.exit(1)

    return move_align_results(ref_dem, dem_tif, save_dir)


def co_registration_multi_process(ref_dem, dem_list, save_dir, process_num):
    print('ref_dem', ref_dem)
    print('source dem:')
    for dem_tif in dem_list:
        print(dem_tif)

    sub_tasks = []
    for dem_tif in dem_list:
        height, width, band_num, daet_type = raster_io.get_height_width_bandnum_dtype(dem_tif)
        # estimate memory need
        need_memory = height*width*4*12 # usually, each pixel need 4 Bytes (float), dem_align.py need more than memory 12 times of file size
        avai_memory = psutil.virtual_memory().available

        while need_memory > avai_memory:
            print('waiting more available memory, need: %.4f GB, available: %.4f GB'%(need_memory/(1000*1000*1000), avai_memory/(1000*1000*1000)))
            time.sleep(10)
            avai_memory = psutil.virtual_memory().available

        while basic.alive_process_count(sub_tasks) >= process_num:
            time.sleep(10)

        sub_process = Process(target=co_registration_one_dem, args=(ref_dem, dem_tif,save_dir))
        sub_process.start()
        sub_tasks.append(sub_process)
        time.sleep(10)   # wait 10 seconds
        if sub_process.exitcode is not None and sub_process.exitcode != 0:
            sys.exit(1)



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
    # co_registration_parallel(ref_dem,dem_list,save_dir,process_num)
    co_registration_multi_process(ref_dem, dem_list, save_dir, process_num)


if __name__ == '__main__':
    usage = "usage: %prog [options] dem_tif_dir or dem_list_txt "
    parser = OptionParser(usage=usage, version="1.0 2020-3-1")
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

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
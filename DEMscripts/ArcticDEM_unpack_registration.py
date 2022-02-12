#!/usr/bin/env python
# Filename: ArcticDEM_unpack_registration 
"""
introduction: unpackage the (strip) ArcticDEM tarball and also apply registration.

# using /usr/local/bin/parallel or something similar to run parallel,
# multiprocessing module (python) may end with deadlock and hang there.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 27 February, 2021
"""

import os,sys
from optparse import OptionParser
import time


deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic

import re

reg_py=os.path.expanduser('~/codes/github_public_repositories/pgcdemtools/apply_setsm_registration.py')

machine_name = os.uname()[1]
# some folder paths
if machine_name == 'uist':
    arcticDEM_reg_tif_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
elif machine_name == 'ubuntu':  # tesia
    arcticDEM_reg_tif_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
else:
    arcticDEM_reg_tif_dir= ''

arcticDEM_mosaic_reg_tif_dir = arcticDEM_reg_tif_dir.replace('registration_tifs','arcticdem_mosaic_reg_tifs')

def is_ArcticDEM_tiles(tar_list):
    '''
    check whether it is mosaic version (tiles) of DEM
    :param tar_list:
    :return:
    '''
    tile_pattern = '^\d{2}_\d{2}_'
    for tar in tar_list:
        tar_base = os.path.basename(tar)
        tiles = re.findall(tile_pattern,tar_base)
        if len(tiles) == 1:
            pass
        else:
            basic.outputlogMessage('%s is not a tile of ArcticDEM'%tar)
            return False

    return True

def get_dem_path_in_unpack_tarball(out_dir, pre_name=None):
    file_end = ['_dem.tif','_reg_dem.tif']   # Arctic strip and tile (mosaic) version
    if pre_name is None:
        pre_name = os.path.basename(out_dir)
    for end in file_end:
        dem_tif = os.path.join(out_dir, pre_name + end)
        if os.path.isfile(dem_tif):
            return dem_tif
    return False

def check_files_existence(dir, pre_name):
    file_pattern = ['*dem_reg.tif', '*reg_dem.tif'] # Arctic strip and tile (mosaic) version
    for pattern in file_pattern:
        file_list = io_function.get_file_list_by_pattern(dir,pre_name + pattern)
        if len(file_list) > 0:
            return True
        else:
            # check if in the archived dir
            if os.path.isdir(arcticDEM_reg_tif_dir):
                file_list_archived = io_function.get_file_list_by_pattern(arcticDEM_reg_tif_dir, pre_name + pattern)
                if len(file_list_archived) > 0:
                    return True

            if os.path.isdir(arcticDEM_mosaic_reg_tif_dir):
                file_list_archived = io_function.get_file_list_by_pattern(arcticDEM_mosaic_reg_tif_dir, pre_name + pattern)
                if len(file_list_archived) > 0:
                    return True

    return False

def arcticDEM_strip_registration(strip_dir):
    command_str = 'python '+ reg_py + ' ' +strip_dir
    basic.os_system_exit_code(command_str)
    end = '_dem_reg.tif'
    reg_tif = os.path.join(strip_dir, os.path.basename(strip_dir) + end)
    if os.path.isfile(reg_tif):
        return reg_tif
    else:
        # with open('no_registration_strips.txt','a') as f_obj:
        #     f_obj.writelines('%s\n'%strip_dir)
        f_obj=open('no_registration_strips.txt','a')
        f_obj.writelines('%s\n'%strip_dir)
        f_obj.close()   # close file manually, to avoid dead lock when using multiple threads
        return None

def process_dem_one_tarball(targz,work_dir,apply_registration):

    out_tif = None
    # unpack, it can check whether it has been unpacked
    out_dir = io_function.unpack_tar_gz_file(targz, work_dir)
    if out_dir is not False:
        out_tif = get_dem_path_in_unpack_tarball(out_dir)
        if os.path.isfile(out_tif):
            # registration for each DEM using dx, dy, dz in *reg.txt file
            if apply_registration:
                reg_tif = arcticDEM_strip_registration(out_dir)
                if reg_tif is None:
                    return None, None
                else:
                    out_tif = reg_tif
        else:
            basic.outputlogMessage('warning, no *_dem.tif in %s' % out_dir)

        # tar_folder_list.append(out_dir)
    else:
        basic.outputlogMessage('warning, unpack %s faild' % targz)

    return out_tif, out_dir


def process_dem_tarball(tar_list, work_dir, save_dir, remove_inter_data=False, rm_tarball=False, apply_registration=False):
    '''
    process dem tarball, unpack, apply registration
    :param tar_list:
    :param work_dir:
    :param save_dir:
    :param remove_inter_data:
    :param apply_registration:
    :return:
    '''

    if os.path.isdir(save_dir) is False:
        io_function.mkdir(save_dir)

    if os.path.isfile('no_registration_strips.txt'):
        no_registration_strips = io_function.read_list_from_txt('no_registration_strips.txt')
    else:
        no_registration_strips = []

    out_dir_list = []
    out_reg_tifs = []
    for idx, targz in enumerate(tar_list):
        tar_base = os.path.basename(targz)[:-7]
        # check if no registraion information for this tarball
        if './'+tar_base in no_registration_strips:
            continue

        if check_files_existence(save_dir,tar_base):
            print("registration result of %s already exists, skip"%targz)
            continue

        # check free disk space
        free_GB = io_function.get_free_disk_space_GB(work_dir)
        total_wait_time = 0
        while free_GB < 50 and total_wait_time < 60*60*12:
            basic.outputlogMessage(' The free disk space (%.4f) is less than 50 GB, wait 60 seconds'%free_GB)
            time.sleep(60)
            total_wait_time += 60
            free_GB = io_function.get_free_disk_space_GB(work_dir)

        out_tif, out_dir = process_dem_one_tarball(targz,work_dir,apply_registration)
        if out_tif is None:
            continue
        out_dir_list.append(out_dir)

        # move file to save_dir
        io_function.movefiletodir(out_tif,save_dir)
        dem_log = os.path.join(out_dir,tar_base+'_dem.log')
        if os.path.isfile(dem_log):
            io_function.movefiletodir(dem_log,save_dir)
        matchtag_tif = os.path.join(out_dir,tar_base+'_matchtag_reg.tif')
        if os.path.isfile(matchtag_tif):
            io_function.movefiletodir(matchtag_tif,save_dir)
        matchtag_tif_log = os.path.join(out_dir, tar_base + '_matchtag.log')
        if os.path.isfile(matchtag_tif_log):
            io_function.movefiletodir(matchtag_tif_log, save_dir)

        out_reg_tifs.append(os.path.join(save_dir, os.path.basename(out_tif)))
        # remove folder
        if remove_inter_data:
            io_function.delete_file_or_dir(out_dir)
        if rm_tarball:
            io_function.delete_file_or_dir(targz)

    # remove folder (in case failed in the previous step)
    if remove_inter_data:
        for dir in out_dir_list:
            if os.path.isdir(dir):
                io_function.delete_file_or_dir(dir)

    return out_reg_tifs

def main(options, args):

    save_dir = options.save_dir
    b_rm_inter = options.remove_inter_data
    b_rm_tarball = options.remove_tarball

    tar_dir = args[0]
    if os.path.isfile(tar_dir):
        tar_list = [tar_dir]
    else:
        tar_list = io_function.get_file_list_by_ext('.gz', tar_dir, bsub_folder=False)
        tar_count = len(tar_list)
        if tar_count < 1:
            raise ValueError('No input tar.gz files in %s' % tar_dir)

    if is_ArcticDEM_tiles(tar_list):
        apply_registration = False
    else:
        apply_registration = True

    work_dir = './'
    process_dem_tarball(tar_list, work_dir, save_dir, remove_inter_data=b_rm_inter, rm_tarball=b_rm_tarball, apply_registration=apply_registration)



if __name__ == '__main__':
    usage = "usage: %prog [options] tarball_dir or a_tarball "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: untar ArcticDEM tarball and apply registration '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")

    parser.add_option("-r", "--remove_inter_data",
                      action="store_true", dest="remove_inter_data",default=False,
                      help="True to remove intermediate data")

    parser.add_option("-t", "--remove_tarball",
                      action="store_true", dest="remove_tarball",default=False,
                      help="if true, to remove tarball after unpacking")

    (options, args) = parser.parse_args()

    # print('remove_inter_data',options.remove_inter_data)
    # print('remove_tarball',options.remove_tarball)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
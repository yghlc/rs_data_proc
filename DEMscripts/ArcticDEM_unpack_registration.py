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



deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function
import basic_src.basic as basic

reg_py=os.path.expanduser('~/codes/github_public_repositories/pgcdemtools/apply_setsm_registration.py')

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
    file_list = io_function.get_file_list_by_pattern(dir,pre_name + '*')
    if len(file_list) > 1:
        return True
    else:
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


def process_dem_tarball(tar_list, work_dir, save_dir, remove_inter_data=False, apply_registration=False):
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

    out_dir_list = []
    out_reg_tifs = []
    for idx, targz in enumerate(tar_list):
        tar_base = os.path.basename(targz)[:-7]

        if check_files_existence(save_dir,tar_base):
            print("registration result of %s already exists, skip"%targz)
            continue

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


    # remove folder (in case failed in the previous step)
    if remove_inter_data:
        for dir in out_dir_list:
            if os.path.isdir(dir):
                io_function.delete_file_or_dir(dir)

    return out_reg_tifs

def main(options, args):

    save_dir = options.save_dir
    b_rm_inter = options.remove_inter_data

    tar_dir = args[0]
    if os.path.isfile(tar_dir):
        tar_list = [tar_dir]
    else:
        tar_list = io_function.get_file_list_by_ext('.gz', tar_dir, bsub_folder=False)
        tar_count = len(tar_list)
        if tar_count < 1:
            raise ValueError('No input tar.gz files in %s' % tar_dir)

    work_dir = './'
    process_dem_tarball(tar_list, work_dir, save_dir, remove_inter_data=b_rm_inter, apply_registration=True)



if __name__ == '__main__':
    usage = "usage: %prog [options] tarball_dir or a_tarball "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: untar ArcticDEM tarball and apply registration '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")

    parser.add_option("-r", "--remove_inter_data",
                      action="store_true", dest="remove_inter_data",default=False,
                      help="True to keep intermediate data")

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
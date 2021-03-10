#!/usr/bin/env python
# Filename: move_old_files_folders.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 08 March, 2021
"""
import os,sys
# import difflib
import time
from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic


machine_name = os.uname()[1]
# some folder paths
if machine_name == 'uist':
    arcticDEM_reg_tif_dir = '/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
elif machine_name == 'ubuntu':  # tesia
    arcticDEM_reg_tif_dir = '/home/lihu9680/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/registration_tifs'
else:
    arcticDEM_reg_tif_dir= ''

time_hour_thr = 1

def check_file_or_dir_is_old(file_folder, time_hour_thr):
    # if not exists, then return False
    if os.path.isfile(file_folder) is False and os.path.isdir(file_folder) is False:
        return False
    now = datetime.now()
    m_time = datetime.fromtimestamp(os.path.getmtime(file_folder))
    print('%s modified time: %s'%(file_folder,str(m_time)))
    diff_time = now - m_time
    diff_time_hour = diff_time.total_seconds()/3600
    if diff_time_hour > time_hour_thr:
        return True
    else:
        return False


def main():

    reg_tif_dir = 'arcticdem_registration_tifs'
    while True:
        print(str(datetime.now()),'start moving or removing files or folders\n\n')
        reg_files = io_function.get_file_list_by_pattern(reg_tif_dir, '*')
        print('registration file count: %d in %s'%(len(reg_files),reg_tif_dir ))
        for file in reg_files:
            if check_file_or_dir_is_old(file,time_hour_thr):
                print('%s is older than %f hours, will be moved to archieved dir'%(file, time_hour_thr))
                io_function.movefiletodir(file, arcticDEM_reg_tif_dir,overwrite=True)

        SETSM_dir = io_function.get_file_list_by_pattern('./', 'SETSM_*2m_v3.0')
        print('SETSM folder count: %d in %s'%(len(SETSM_dir),'./' ))
        for folder in SETSM_dir:
            if check_file_or_dir_is_old(folder,time_hour_thr):
                print('%s is older than %f hours, will be removed'%(folder, time_hour_thr))
                io_function.delete_file_or_dir(folder)

        grid_tmp_dir = io_function.get_file_list_by_pattern('./', 'grid*files')
        print('Grid tmp folder count: %d in %s'%(len(grid_tmp_dir),'./'))
        for folder in grid_tmp_dir:
            if check_file_or_dir_is_old(folder, time_hour_thr):
                print('%s is older than %f hours, will be removed' % (folder, time_hour_thr))
                io_function.delete_file_or_dir(folder)

        open_files = basic.get_all_processes_openfiles('python')
        print('open_files count by process with python name: %d'%len(open_files))
        for item in open_files:
            print(item)

        time.sleep(60)  # wait

    pass


if __name__ == '__main__':
    main()
    pass
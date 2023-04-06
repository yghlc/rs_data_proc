#!/usr/bin/env python
# Filename: group_s1_orbit_files.py 
"""
introduction:  1. in a folder with S1*.zip, run  eof (https://github.com/scottstanie/sentineleof)
               2. after download EOF (orbit) file, run this script in the same folder
               # this script will organize the orbit file into /home/lingcao/.snap/auxdata/Orbits/Sentinel-1/POEORB for SNAP

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 05 April, 2023
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic
from datetime import datetime

orbit_folder = os.path.expanduser('~/.snap/auxdata/Orbits/Sentinel-1')

def move_a_eof_file(file_path):
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    #e.g. S1B_OPER_AUX_POEORB_OPOD_20210312T044626_V20180216T225942_20180218T005942.EOF
    info_list = file_name.split('_')
    year = info_list[7][:4]
    month = info_list[7][4:6]
    dst_dir = os.path.join(orbit_folder,info_list[3], info_list[0], year,month) # POEORB, S1B,
    if os.path.isdir(dst_dir) is False:
        io_function.mkdir(dst_dir)

    # io_function.movefiletodir(file_path,dst_dir)
    # zip it, instead of moving it
    zip_path = os.path.join(dst_dir,os.path.basename(file_path)+'.zip')
    if os.path.isfile(zip_path):
        print('%s already exist, skip'%zip_path)
        return
    cmd_str = 'zip %s %s'%(zip_path,file_path)
    basic.os_system_exit_code(cmd_str)
    print(datetime.now(),'compress a EOF file into %s'%zip_path)



def main():
    eof_list = io_function.get_file_list_by_ext('.EOF','./',bsub_folder=False)
    for eof in eof_list:
        move_a_eof_file(eof)

    print('organize %d orbit files'%len(eof_list))

    pass

if __name__ == '__main__':
    main()
    pass
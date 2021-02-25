#!/usr/bin/env python
# Filename: backup_curc_scripts.py
"""
introduction: backup scripts

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 25 February, 2021
"""

import os,sys
# import difflib

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function

def main():
    org_dir='/scratch/summit/lihu9680'
    bak_dir='/home/lihu9680/scripts_para_bak'

    # get bash file (*.sh) list
    print('find *.sh files')
    sh_list = io_function.get_file_list_by_ext('.sh',org_dir,bsub_folder=True)
    changed_list = []
    for sh in sh_list:
        new_path = bak_dir + sh.replace(org_dir,'')
        # print(bak_path)
        mo_time = os.path.getmtime(sh)
        if os.path.isfile(new_path):
            bak_mo_time = os.path.getmtime(new_path)
            if mo_time == bak_mo_time:
                continue

        new_dir = os.path.dirname(new_path)
        if os.path.isdir(new_dir) is False:
            io_function.mkdir(new_dir)

        res = os.system('cp -p %s %s'%(sh, new_path))
        if res!=0:
            sys.exit(1)
        changed_list.append(sh)

    if len(changed_list) < 1 :
        print("no new or modified files")
    else:
        print("backup files:")
        for sh in changed_list:
            print(sh)



if __name__ == '__main__':
    main()
    pass



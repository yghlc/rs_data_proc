#!/usr/bin/env python
# Filename: back_curc_scripts.py 
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
    sh_list = io_function.get_file_list_by_ext('.sh',org_dir,bsub_folder=True)
    changed_list = []
    for sh in sh_list:
        bak_path = bak_dir + sh.replace(org_dir,'')
        mo_time = os.path.getmtime(sh)
        if os.path.isfile(bak_path):
            bak_mo_time = os.path.getmtime(bak_path)
            if mo_time == bak_mo_time:
                continue

        bak_dir = os.path.dirname(bak_path)
        if os.path.isdir(bak_dir) is False:
            io_function.mkdir(bak_dir)

        res = os.system('cp -p %s %s'%(sh, bak_path))
        if res!=0:
            sys.exit(1)
        changed_list.append(sh)

    if len(changed_list) != 0:
        print("no new or modified files")
    else:
        print("backup files:")
        for sh in changed_list:
            print(sh)



if __name__ == '__main__':
    main()
    pass



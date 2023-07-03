#!/usr/bin/env python
# Filename: coregistration_siftGPU.py
"""
introduction: co-register remote sensing imagery using SiftGPU

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 28 June, 2023
"""

import os,sys
from optparse import OptionParser
from datetime import datetime
import time
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.basic as basic
from basic_src.xml_rw import OffsetMetaDataClass
import basic_src.io_function as io_function
import basic_src.RSImageProcess as RSImageProcess

def registration_two_images(ref_image, warp_image):
    io_function.is_file_exist(ref_image)
    io_function.is_file_exist(warp_image)
    ref_name = io_function.get_name_no_ext(ref_image)
    warp_name = io_function.get_name_no_ext(warp_image)
    meta_xml = 'registration' + '_' + ref_name + '_' + warp_name + '.xml'
    meta_obj = OffsetMetaDataClass(meta_xml)
    bkeepmidfile = True
    RSImageProcess.coregistration_siftGPU(ref_image,warp_image,bkeepmidfile, meta_obj)

def main(options, args):

    ref_img_path = args[0]
    warp = args[1]
    if warp.endswith('.txt'):
        warp_img_list = io_function.read_list_from_txt(warp)
    else:
        warp_img_list = [warp]

    for idx, warp_img in enumerate(warp_img_list):
        if os.path.basename(warp_img) == os.path.basename(ref_img_path):
            basic.outputlogMessage('the reference and warp image is the same: %s, skip'%os.path.basename(warp_img))
            continue
        basic.outputlogMessage(' (%d/%d) registration two images: %s vs %s'
              %(idx+1, len(warp_img_list), os.path.basename(ref_img_path), os.path.basename(warp_img)))
        curr_dir = os.getcwd()
        # working a sub-directory
        work_dir = os.path.join(curr_dir,'coreg_%d_%s'%((idx+1),io_function.get_name_no_ext(warp_img)))
        io_function.mkdir(work_dir)
        os.chdir(work_dir)
        ref_img_path = os.path.abspath(ref_img_path)
        warp_img = os.path.abspath(warp_img)
        registration_two_images(ref_img_path, warp_img)
        os.chdir(curr_dir)

        # copy file
        out_warp_file = os.path.join(work_dir,os.path.basename(io_function.get_name_by_adding_tail(warp_img,'new_warp')))
        copy_dest = io_function.get_name_by_adding_tail(warp_img,'coreg')
        io_function.copy_file_to_dst(out_warp_file,copy_dest)



    pass

if __name__ == '__main__':
    usage = "usage: %prog [options] reference warp_image( or warp_image_list.txt) "
    parser = OptionParser(usage=usage, version="1.0 2023-6-28")
    parser.description = 'Introduction: co-registration of remote sensing imagery   '

    (options, args) = parser.parse_args()

    if len(sys.argv) < 2 or len(args) < 2:
        parser.print_help()
        sys.exit(2)

    main(options, args)

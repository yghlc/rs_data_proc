#!/usr/bin/env python
# Filename: crop_image_min 
"""
introduction: crop a image based on the extent of another images

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 29 October, 2021
"""

import os,sys
from optparse import OptionParser

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import vector_gpd

import basic_src.RSImageProcess as RSImageProcess
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function


def subset_image_by_ref_image(in_img, out_img, base_image,resample_m='bilinear',o_format='GTiff'):
    return RSImageProcess.subset_image_baseimage(out_img,in_img,base_image,same_res=True,resample_m=resample_m)


def main(options, args):

    ref_raster = args[0]
    img_path = args[1]
    save_dir = options.save_dir
    save_path = options.save_path

    #check projection
    extent_prj = map_projection.get_raster_or_vector_srs_info_proj4(ref_raster)
    img_prj = map_projection.get_raster_or_vector_srs_info_proj4(img_path)
    if img_prj != extent_prj:
        raise ValueError('Project of %s and %s is different'%(ref_raster, img_path))

    if save_path is None:
        out_img = io_function.get_name_by_adding_tail(img_path,'sub')
        out_img = os.path.join(save_dir, os.path.basename(out_img))
    else:
        out_img = save_path

    subset_image_by_ref_image(img_path,out_img,ref_raster,resample_m='near')


    pass

if __name__ == "__main__":

    usage = "usage: %prog [options] ref_raster image_path"
    parser = OptionParser(usage=usage, version="1.0 2021-10-29")
    parser.description = 'Introduction: crop image to the extent of another raster   '

    parser.add_option("-d", "--save_dir",
                      action="store", dest="save_dir",default='./',
                      help="the folder to save pre-processed results")

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="the output image path")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if  len(args) < 2:
        parser.print_help()
        sys.exit(2)

    main(options, args)




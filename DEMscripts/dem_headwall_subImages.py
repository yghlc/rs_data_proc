#!/usr/bin/env python
# Filename: dem_headwall_subImages 
"""
introduction: Base on the shapefile of headwall, get subImages of Hillshade or other data

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 30 June, 2021
"""

import os,sys
deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)

import basic_src.io_function as io_function

# path of headwall shapefiles
from dem_common import dem_headwall_shp_dir
# from dem_common import grid_dem_headwall_shp_dir

from dem_common import dem_hillshade_dir,dem_hillshade_subImages_headwall

extract_py = os.path.expanduser('~/codes/PycharmProjects/ChangeDet_DL/dataTools/extract_subimage_timeSeries.py')

def set_image_dir_patter_description(tif_dir,shapefile):

    # for ones like: SETSM_WV03_20170813_104001003237BC00_1040010031CAA600_seg4_2m_v3.0_dem_reg_slope_headwall.shp
    basename = os.path.basename(shapefile).split('_slope_')[0]
    tifs = io_function.get_file_list_by_pattern(tif_dir,basename+'*')
    if len(tifs) == 1:
        img_dir = os.path.dirname(tifs[0])
        img_name = os.path.basename(tifs[0])
        with open('image_dir_patter_description.txt','w') as f_obj:
            f_obj.writelines(img_dir+':'+img_name+':'+'hillshade\n')
        return True
    elif len(tifs) < 1:
        print('Error, No tif file for %s in %s'%(shapefile,tif_dir))
        return False
    else:
        with open('image_dir_patter_description.txt','w') as f_obj:
            f_obj.writelines(tif_dir+':'+basename+':'+'hillshade\n')
        return True

def main():

    # get shapefile list
    headwall_shp_list = io_function.get_file_list_by_ext('.shp',dem_headwall_shp_dir,bsub_folder=False)
    if len(headwall_shp_list) < 1:
        raise ValueError('NO shapefile in %s'%dem_headwall_shp_dir)

    failed_shp = []
    out_dir = dem_hillshade_subImages_headwall
    if len(sys.argv) == 2:
        # change the output dir
        out_dir = sys.argv[1]

    for idx, shp in enumerate(headwall_shp_list):
        print('(%d/%d) extract sub images for %s'%(idx,len(headwall_shp_list), shp))

        if set_image_dir_patter_description(dem_hillshade_dir,shp) is False:
            continue
        
        save_dir = os.path.join(out_dir,os.path.splitext(os.path.basename(shp))[0])
        if os.path.isdir(save_dir):
            print('Warning, skip due to subImages for %s may exist'%shp)
            continue
        io_function.mkdir(save_dir)

        res = os.system(extract_py + ' -p para_file_subImage.ini -o %s '%save_dir + shp)
        if res !=0:
            failed_shp.append(shp)


    if len(failed_shp) > 0:
        io_function.save_list_to_txt('failed_shp.txt',failed_shp)

        



if __name__ == '__main__':
    print(sys.argv[0],sys.argv[1])
    print(len(sys.argv))
    # main()
    pass
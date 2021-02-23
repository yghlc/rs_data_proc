#!/usr/bin/env python
# Filename: bash_orthorectified.py 
"""
introduction: ortho-rectified using gdalwarp

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 23 February, 2021
"""

import os, sys

# path of DeeplabforRS
codes_dir2 = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, codes_dir2)
import basic_src.basic as basic
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection


dir='/BhaltosMount/Bhaltos/GREENLAND/2021_NNA_PROJECT/PERMAFROST/ITTOQQOTROORMIIT'

# -te xmin ymin xmax ymax
# extent from "ogrinfo extent_ITTO.shp extent_ITTO"  ~20 by 20 km
# extent="451984.972763 7810372.206539 476010.392359, 7830488.990841"

#xmin ymin xmax ymax from "ogrinfo ITTOQQOTROORMIIT_ext.shp ITTOQQOTROORMIIT_ext"   1.5 by 0.7 km
extent="463513.839669 7819968.221432 464884.176408 7820678.186371"

out_res=0.5
thread_num=8

def ortho_rectified_gdalwarp(input, output, dem_tif):

    #     gdalwarp -co compress=lzw -co tiled=yes -co bigtiff=if_safer -tr ${res} ${res} -t_srs ${prj} \
#           -multi -wo NUM_THREADS=8  -r cubic -et 0.01 -to "RPC_DEM=${dem}" ${img}  ${out}
    CommandString = 'gdalwarp -co compress=lzw -co tiled=yes -co bigtiff=if_safer -tr '+ str(out_res) + ' '+ str(out_res)

    dem_prj = map_projection.get_raster_or_vector_srs_info_epsg(dem_tif)
    CommandString += ' -t_srs ' + dem_prj
    CommandString += ' -te ' + extent
    CommandString += ' -multi -wo NUM_THREADS=' + str(thread_num)
    CommandString += ' -r cubic -et 0.01 -to "RPC_DEM=%s" '%dem_tif
    CommandString += ' -dstnodata 0 '
    CommandString += ' %s %s  '%(input,output)

    return basic.exec_command_string_one_file(CommandString, output)


def main():
    ntf_list = io_function.get_file_list_by_ext('.ntf',os.path.join(dir,'DATA'), bsub_folder=True)
    io_function.save_list_to_txt('ntf_list.txt',ntf_list)
    dem_list = io_function.get_file_list_by_ext('.tif',os.path.join(dir,'PRODUCTS'), bsub_folder=True)
    dem_list = [ item for item in dem_list if item.endswith('_dem.tif') and 'strips' in item]
    io_function.save_list_to_txt('dem_list.txt', dem_list)

    for idx, ntf in enumerate(ntf_list):
        print(' (%d/%d) working on '%(idx+1, len(ntf_list)),ntf)
        name = os.path.basename(ntf)
        scene_id = name.split('_')[2]
        print('scene_id:', scene_id)

        dem_path = None
        for dem_tif in dem_list:
            if scene_id in os.path.basename(dem_tif):
                dem_path = dem_tif
                break
        if dem_path is None:
            raise ValueError('Cannot find the corresponding DEM')

        output = os.path.splitext(name)[0] + '_ortho_sub.tif'
        ortho_rectified_gdalwarp(ntf, output, dem_path)
        # break

    pass


if __name__ == '__main__':
    main()
    pass
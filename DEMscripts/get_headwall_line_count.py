#!/usr/bin/env python
# Filename: get_headwall_line_count.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 October, 2022
"""

import os,sys
from optparse import OptionParser
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))

import basic_src.io_function as io_function
import basic_src.basic as basic

# import geopandas as gpd

def main(options, args):
    dir = args[0]
    shp_list = io_function.get_file_list_by_pattern(dir,'dem_headwall_shp_grid/*/*.shp')
    shp_list = [item for item in shp_list if 'rippleSel' not in os.path.basename(item)]
    print('find %s headwall lines file in %s'%(len(shp_list), dir))

    total_num = 0
    for shp in shp_list:
        # it look like geopandas will read all the data from the disk
        # shapefile = gpd.read_file(shp)
        # count = len(shapefile)

        layer_name = os.path.splitext(os.path.basename(shp))[0] 
        command_str = 'ogrinfo -so %s %s | grep Feature'%(shp,layer_name)
        status, result = basic.exec_command_string(command_str)
        #print('status, result',status, result)
        if status != 0:
            print(result)
            sys.exit(status)
        # Feature Count: 28992
        _, num_str = result.split(':')
        count = int(num_str.strip())
        print('line count in %s is %d'%(shp,count))
        total_num += count

    # save the total count to file
    with open('headwall-line-count.txt', 'a') as f_obj:
        f_obj.writelines('%s, Feature count: %d\n'%(dir,total_num))


if __name__ == '__main__':
    usage = "usage: %prog [options] dir  "
    parser = OptionParser(usage=usage, version="1.0 2022-10-1")
    parser.description = 'Introduction: get the features counts  '

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
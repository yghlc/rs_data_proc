#!/usr/bin/env python
# Filename: ArcticDEM_proc_grid.py
"""
introduction: produce elevation differences from  multi temporal ArcticDEM

# I divide the coverage of ArcticDEM to 20 km by 20 km grids:
~/Data/Arctic/ArcticDEM/grid_shp/ArcticDEM_grid_20km.shp, each grid have a ID, in projection fo EPSG:3413, Polar Stereographic North.

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 06 March, 2021
"""

import os,sys
from optparse import OptionParser


deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)




def main(options, args):
    extent_shp = args[0]

    # get grid ids based on input extent

    # download ArcticDEM and applying registration

    # mosaic and crop

    # dem co-registration (cancel, the result in not good with the default setting)

    # dem differencing

    # reprojection if necessary


if __name__ == '__main__':
    usage = "usage: %prog [options] extent_shp or grid_id_list.txt "
    parser = OptionParser(usage=usage, version="1.0 2020-12-26")
    parser.description = 'Introduction: produce DEM difference from multiple temporal ArcticDEM  '


    pass
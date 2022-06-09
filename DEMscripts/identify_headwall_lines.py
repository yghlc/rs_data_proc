#!/usr/bin/env python
# Filename: dem_headwall_pattern.py 
"""
introduction: base on the "headwall" lines extracted from slope files, try to identiy some lines
that are really headwall lines (most are not).
if we can find the real headwall lines, potential, we can group them and identify the locations of thaw slumps

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 08 June, 2022
"""

import os,sys
from optparse import OptionParser
import time
from datetime import datetime
machine_name = os.uname()[1]

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import raster_io
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

from multiprocessing import Pool
from shapely.strtree import STRtree

# object.parallel_offset

def calculate_one_hausdorff_dis_closest(center_obj, geometry_list, max_extent=None):
    t0 = time.time()
    object_list = [ item for item in geometry_list if item!=center_obj] # remove itself, otherwise, all dis will be 0
    # tree = STRtree(object_list)
    if max_extent is not None:
        # find objects within extent
        # consider using clip_by_rect() to clip the data
        # object_list = vector_gpd.get_poly_index_within_extent() # not sure how to do that
        pass

    # nearest_one = tree.nearest(center_obj)
    # dis_list = [ center_obj.hausdorff_distance(nearest_one) ]

    dis_list = [center_obj.hausdorff_distance(item) for item in object_list ]
    print(datetime.now(), 'time cost for calculate_one_hausdorff_dis_closest: %f seconds' % (time.time() - t0))
    return min(dis_list)

def calculate_hausdorff_dis(line_list, max_extent = 100, process_num = 1):
    '''
    calculate the hausdorff distance between lines,  only calculate the one between other within "max_extent" meters
    :param line_list: a list of lines (it's ok if the input is polygons)
    :param max_extent: max extent, only calculate them within this extent
    :return:
    '''
    if len(line_list) < 1:
        basic.outputlogMessage('Warning, No geometry')
        return False

    ########## Not test yet ##########

    t0=time.time()

    if process_num==1:
        hausdorff_dis_nearest = []
        for idx, a_line in enumerate(line_list):
            dis = calculate_one_hausdorff_dis_closest(a_line,line_list,max_extent=max_extent)
            hausdorff_dis_nearest.append(dis)
    elif process_num > 1:
        hausdorff_dis_nearest = []
        pass
    else:
        raise ValueError('Wrong process number: %s'%str(process_num))

    print(datetime.now(),'time cost for calculate_hausdorff_dis: %f seconds'%(time.time()-t0))
    print(hausdorff_dis_nearest)


def test_calculate_hausdorff_dis():
    data_dir = os.path.expanduser('~/Data/dem_processing/Alaska_grid10741_results')
    lines_shp = os.path.join(data_dir,'dem_headwall_shp_grid/headwall_shps_grid10741/headwall_shp_multiDates_10741_subset.shp')
    # lines_shp = os.path.join(data_dir,'dem_headwall_shp_grid/headwall_shps_grid10741/headwall_shp_multiDates_10741.shp')

    line_list, dem_year_list = vector_gpd.read_lines_attributes_list(lines_shp, 'dem_year')
    calculate_hausdorff_dis(line_list, max_extent=100, process_num=1)
    pass


def main(options, args):
    test_calculate_hausdorff_dis()

    # lines_shp = args[0]
    # print(lines_shp)

    # read the vector files
    # line_list, dem_year_list = vector_gpd.read_lines_attributes_list(lines_shp,'dem_year')







if __name__ == '__main__':
    usage = "usage: %prog [options] lines_multiTemporal.shp "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: identify real headwall lines and thaw slumps  '

    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to create the mosaic")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
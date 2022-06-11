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
import pandas as pd
import geopandas as gpd
import raster_io
import basic_src.basic as basic
import basic_src.timeTools as timeTools
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

from multiprocessing import Pool
from shapely.strtree import STRtree
import numpy as np

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

def one_line_ripple(id, line, lTime, dataframe, delta=2, total_steps=50, max_extent=100, sim_range=[0.5, 2]):

    # clip
    # max_extent = 50
    mask = line.buffer(max_extent + 1)

    # got an error: IllegalArgumentException: Argument must be Polygonal or LinearRing
    # uninstall pygeos on my Mac fix this problem
    # clip cut some lines and create new lines, no good.
    # clip_dataframe = gpd.clip(dataframe,mask,keep_geom_type=True)   # clip(gdf, mask, keep_geom_type=False)

    intersects = dataframe.intersects(mask)# (dataframe,mask,keep_geom_type=True)
    select_dataframe = dataframe[intersects]

    # remove the "line" itself
    select_dataframe.set_index('id',inplace=True)
    select_dataframe = select_dataframe.drop(id) # ,inplace=True

    # save results to check.
    # # clip_dataframe.to_file('clip_results_%d.shp'%id, driver='ESRI Shapefile')  # save to checking
    # select_dataframe.to_file('intersect_results_%d.shp' % id, driver='ESRI Shapefile')  # save to checking
    # m_pd = pd.DataFrame({'ext_polygon':[mask]})
    # vector_gpd.save_lines_to_files(m_pd,'ext_polygon',select_dataframe.crs,'masks_%d.shp'%id)


    # plot figures to check
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(figsize=(12, 8))
    # clip_dataframe.plot(ax=ax, color="purple")
    # ax.set_title(" Clipped", fontsize=20)
    # plt.show()
    # sys.exit(1)



    # buffer operations
    buffer_dis = 0.1
    # only count these lines 1: have similar length with the center line, 2.not in the same year has been recorded.
    line_count_per_step = np.zeros((total_steps,),dtype=np.uint16)
    previous_inter = pd.Series([False]*len(select_dataframe), index=select_dataframe.index,dtype=bool)
    min_length = line.length * sim_range[0]
    max_length = line.length * sim_range[1]
    recorded_Times = [lTime]

    for step in range(total_steps):
        dis = buffer_dis + delta*step
        ripple_poly = line.buffer(dis)

        b_inters = select_dataframe.intersects(ripple_poly) # pandas.Series type, True or False
        # b_new = b_inters.drop_duplicates(previous)
        # b_new = b_inters.compare(previous)
        # b_new = b_inters.ne(previous)
        b_inters[previous_inter] = False    # remove previous ones
        # print(b_inters)
        if not b_inters.any():  # "is False" not working
            continue

        new_dataframe = select_dataframe[b_inters]   # find new intersected lines
        new_line_count = 0
        # 1. have similar length with the center line,
        # 2.not in the same year has been recorded.
        for idx, row in new_dataframe.iterrows():
            if min_length <= row['length_m'] <= max_length and row['dem_year'] not in recorded_Times:
                new_line_count += 1
                recorded_Times.append(row['dem_year'])
        line_count_per_step[step] = new_line_count

        # update previous iteration
        previous_inter[b_inters] = True
        if previous_inter.all():  # if all of them have been checked, quit
            break

    print(line_count_per_step)





def line_ripple_statistics(lines_multiTemporal_path, delta=2, total_steps=50, max_extent=100, sim_range=[0.5, 2],process_num=1):
    '''
    statistic lines information soemthing like a ripple using buffer operation
    :param lines_multiTemporal_path: input path
    :param delta: the increase of each buffer operation
    :param total_steps: run buffer operation for  total_steps times
    :param max_extent: max extent of a ripple
    :param sim_range: if a line is within a range of sim_range[0]*length to sim_range[1]*length, then consider them as similar
    :param process_num:
    :return:
    '''

    # read the data frame (lines and all attributes)
    io_function.is_file_exist(lines_multiTemporal_path)
    line_dataframe = gpd.read_file(lines_multiTemporal_path)

    crop_ext = max(delta*total_steps, max_extent)

    # initial
    for ri, row in line_dataframe.iterrows():
        line = row['geometry']
        if row['id'] != 709:
            continue
        one_line_ripple(row['id'],line, row['dem_year'], line_dataframe,delta=delta, total_steps=total_steps, max_extent=crop_ext,sim_range=sim_range)





def test_line_ripple_statistics():
    data_dir = os.path.expanduser('~/Data/dem_processing/Alaska_grid10741_results')
    lines_shp = os.path.join(data_dir,'dem_headwall_shp_grid/headwall_shps_grid10741/headwall_shp_multiDates_10741_subset.shp')

    line_ripple_statistics(lines_shp, delta=2, total_steps=50, max_extent=100, process_num=1)


def main(options, args):
    # test_calculate_hausdorff_dis()
    test_line_ripple_statistics()

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
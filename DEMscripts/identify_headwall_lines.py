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
    # previous_inter = pd.Series([False]*len(select_dataframe), index=select_dataframe.index,dtype=bool)
    min_length = line.length * sim_range[0]
    max_length = line.length * sim_range[1]
    recorded_Times = [lTime]

    for step in range(total_steps):
        dis = buffer_dis + delta*step
        ripple_poly = line.buffer(dis)

        b_inters = select_dataframe.intersects(ripple_poly) # pandas.Series type, True or False
        if not b_inters.any():  # "is False" not working
            continue

        new_dataframe = select_dataframe[b_inters]   # find new intersected lines
        # re-calculate the length within "ripple_poly"
        # to avoid the problem that some portion of headwall don't retreat but other parts retreat over years.
        inters_lines = new_dataframe.intersection(ripple_poly)
        # got SettingWithCopyWarning when adding a new column to this
        # new_dataframe.loc[:,'recal_len'] = inters_lines.length      # added length within the buffer zone

        new_line_count = 0
        rm_index = []
        # 1. have similar length (within buffer) with the center line,
        # 2.not in the same year has been recorded.
        for idx, row in new_dataframe.iterrows():
            if row['dem_year'] in recorded_Times:
                rm_index.append(idx)
                continue
            if min_length <= inters_lines.length[idx] <= max_length:
                new_line_count += 1
                recorded_Times.append(row['dem_year'])
                rm_index.append(idx)
            elif inters_lines.length[idx] >= row['length_m']:
                rm_index.append(idx)
            else:
                pass
        line_count_per_step[step] = new_line_count

        # remove the lines have been checked
        select_dataframe = select_dataframe.drop(rm_index)
        if len(select_dataframe) < 1:  # if all of them have been checked, quit
            break

    # print(line_count_per_step)
    # print(recorded_Times,pd.Series(recorded_Times).is_monotonic_increasing, pd.Series(recorded_Times).is_monotonic_decreasing)  # be monotonically increasing

    # return np.sum(line_count_per_step) + 1       # +1 itself back
    # return len(recorded_Times)                     # same to the one above

    # calculate ripple attributes
    line_count_per_step[0] = 1  # add the started line back
    # # min, max, and average distance of each ripple lines
    non_zero_idx = np.where(line_count_per_step > 0)[0]
    if len(non_zero_idx) <= 1:
        min_ripple_delta = None
        max_ripple_delta = None
        avg_ripple_delta = None
    else:
        ripple_distance = non_zero_idx - np.roll(non_zero_idx,1)    # roll offset one element
        ripple_distance = ripple_distance[1:]
        if line_count_per_step.max() > 1:       # if two line overlap
            min_ripple_delta = 0
        else:
            min_ripple_delta = ripple_distance.min()
        max_ripple_delta = ripple_distance.max()
        avg_ripple_delta = ripple_distance.mean()


    #  be monotonically increasing or decreasing
    time_series = pd.Series(recorded_Times)
    b_mono_increase = time_series.is_monotonic_increasing   # likely an oldest headwall line
    b_mono_decrease = time_series.is_monotonic_decreasing   # likely the most recent headwall line

    return len(recorded_Times),b_mono_increase,b_mono_decrease,min_ripple_delta,max_ripple_delta,avg_ripple_delta


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
    total_count = len(line_dataframe)

    ripple_count_list = []
    b_mono_increase_list = []
    b_mono_decrease_list = []
    min_ripple_delta_list = []
    max_ripple_delta_list = []
    avg_ripple_delta_list = []

    if process_num==1:
        for ri, row in line_dataframe.iterrows():
            ripple_count, b_mono_increase, b_mono_decrease, min_ripple_delta, max_ripple_delta, avg_ripple_delta \
                = one_line_ripple(row['id'],row['geometry'], row['dem_year'], line_dataframe,delta=delta,
                                  total_steps=total_steps, max_extent=crop_ext,sim_range=sim_range)
            ripple_count_list.append(ripple_count)
            b_mono_increase_list.append(b_mono_increase)
            b_mono_decrease_list.append(b_mono_decrease)
            min_ripple_delta_list.append(min_ripple_delta)
            max_ripple_delta_list.append(max_ripple_delta)
            avg_ripple_delta_list.append(avg_ripple_delta)

            print('%d/%d'%(ri,total_count),ripple_count, b_mono_increase, b_mono_decrease, min_ripple_delta, max_ripple_delta, avg_ripple_delta)
    elif process_num > 1:
        threadpool = Pool(process_num)
        para_list = [(row['id'],row['geometry'], row['dem_year'],line_dataframe,delta,total_steps,crop_ext,sim_range) for ri, row in line_dataframe.iterrows()]
        stats_res_list = threadpool.starmap(one_line_ripple, para_list)
        threadpool.close()
        for res in stats_res_list:
            ripple_count_list.append(res[0])
            b_mono_increase_list.append(res[1])
            b_mono_decrease_list.append(res[2])
            min_ripple_delta_list.append(res[3])
            max_ripple_delta_list.append(res[4])
            avg_ripple_delta_list.append(res[5])
    else:
        raise ValueError('uknown process_num %s'%str(process_num))

    # save attributes into shapefile
    add_attributes = {'ri_count':ripple_count_list,
                     'mono_incre':b_mono_increase_list,
                     'mono_decre':b_mono_decrease_list,
                     'minRdelta':min_ripple_delta_list,
                     'maxRdelta':max_ripple_delta_list,
                     'avgRdelta':avg_ripple_delta_list}
    vector_gpd.add_attributes_to_shp(lines_multiTemporal_path,add_attributes)
    basic.outputlogMessage('Save ripple attributes into %s'%lines_multiTemporal_path)



def test_line_ripple_statistics():
    data_dir = os.path.expanduser('~/Data/dem_processing/Alaska_grid10741_results')
    lines_shp = os.path.join(data_dir,'dem_headwall_shp_grid/headwall_shps_grid10741/headwall_shp_multiDates_10741_subset.shp')

    line_ripple_statistics(lines_shp, delta=2, total_steps=50, max_extent=100, process_num=1)


def main(options, args):
    # test_calculate_hausdorff_dis()
    # test_line_ripple_statistics()
    t0 = time.time()

    lines_shp = args[0]
    process_num = options.process_num
    buffer_delta = options.buffer_delta
    total_steps = options.total_steps
    max_extent = options.max_extent
    lower_similarity = options.lower_similarity
    upper_similarity = options.upper_similarity

    sim_range = [lower_similarity,upper_similarity]

    # print(lines_shp)

    # read the vector files
    line_ripple_statistics(lines_shp, delta=buffer_delta, total_steps=total_steps, max_extent=max_extent,
                           sim_range=sim_range, process_num=process_num)

    print('total time cost of identify_headwall_lines.py', time.time() - t0, 'seconds')







if __name__ == '__main__':
    usage = "usage: %prog [options] lines_multiTemporal.shp "
    parser = OptionParser(usage=usage, version="1.0 2021-3-6")
    parser.description = 'Introduction: identify real headwall lines and thaw slumps  '

    parser.add_option("-d", "--buffer_delta",
                      action="store", dest="buffer_delta", type=float, default=2.0,
                      help="the delta of each step of buffer operation")

    parser.add_option("-t", "--total_steps",
                      action="store", dest="total_steps", type=int, default=50,
                      help="the total steps of buffer operation")

    parser.add_option("-e", "--max_extent",
                      action="store", dest="max_extent", type=float, default=100.0,
                      help="the maximum extent to check the surrounding area of a line")

    parser.add_option("-l", "--lower_similarity",
                      action="store", dest="lower_similarity", type=float, default=0.5,
                      help="the lower range to check if the length of a line is similar")

    parser.add_option("-u", "--upper_similarity",
                      action="store", dest="upper_similarity", type=float, default=2.0,
                      help="the upper range to check if the length of a line is similar")


    parser.add_option("", "--process_num",
                      action="store", dest="process_num", type=int, default=4,
                      help="number of processes to run")

    (options, args) = parser.parse_args()
    # print(options.create_mosaic)

    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
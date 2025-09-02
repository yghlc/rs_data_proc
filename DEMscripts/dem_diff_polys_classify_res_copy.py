#!/usr/bin/env python
# Filename: dem_diff_polys_classify_res_copy.py 
"""
introduction: dem_diff_polys_classify.py saved the results into shapefile but only keep the centroid of original polygons,

This script to get the polygons with class of 1 (or non 0) and merge to new file for each ext

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 02 September, 2025
"""

import os,sys
from optparse import OptionParser
import time

from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import basic_src.io_function as io_function
import basic_src.basic as basic

from dem_common import get_grid_id_from_path
# ArcticDEM results on ygAlpha
ArcticDEM_results_dir = os.path.expanduser('~/data1/ArcticDEM_results')

import pandas as pd
import geopandas as gpd

def merge_two_shapefile(sam_shp, classify_res_shp, keep_colums=['preClassID'], values=[1]):
    # t0=time.time()
    sam_res = gpd.read_file(sam_shp)
    class_res = gpd.read_file(classify_res_shp)
    if sam_res.crs != class_res.crs:
        raise ValueError(f'The projection of {sam_shp} and {classify_res_shp} is different')

    if len(sam_res) != len(class_res):
        raise ValueError(f'The row count in  {sam_shp} and {classify_res_shp} is different')

    if len(keep_colums) != len(values):
        raise ValueError(f'The count in keep_colums ({keep_colums}) and values ({values}) should be the same')

    # t1 = time.time()
    # merge sam_res and class_res,and only keep the columns in keep_colums and values in [1]
    merged = sam_res.join(class_res.drop(columns='geometry'), lsuffix='_sam', rsuffix='_class')
    # print(merged)

    # t2 = time.time()

    # Only keep the columns in keep_colums (if they exist in merged)
    columns_to_keep = [col for col in keep_colums if col in merged.columns]
    # print(columns_to_keep)
    if not columns_to_keep:
        raise ValueError(f'None of the columns in {keep_colums} found in merged shapefile.')

    # Always keep geometry column
    columns_to_keep = columns_to_keep + ['geometry']

    filtered = merged[columns_to_keep]
    # print(filtered)
    # t3 = time.time()
    # Filter rows: only keep those where all keep_colums have values in `values`
    mask = filtered[keep_colums].apply(lambda row: all(val in values for val in row), axis=1)
    filtered = filtered[mask]
    t4 = time.time()

    # print('time cost:', t1-t0, t2-t1, t3-t2, t4-t3)

    return filtered


def test_merge_two_shapefile():
    sam_shp = os.path.expanduser('~/Downloads/tmp/DEM_diff_grid16508_sam_seg.shp')
    class_shp = os.path.expanduser('~/Downloads/tmp/grid16508_DEMdiffColor_2008-17-predicted_classID.shp')

    save_path = os.path.expanduser('~/Downloads/tmp/grid16508_sam_seg_class_1.shp')
    merge_gpd = merge_two_shapefile(sam_shp, class_shp)
    merge_gpd.to_file(save_path)


def copy_classify_result_to_sam_result(ext_str, dem_diff_poly_dir,dem_diff_poly_classify_res_dir):

    done_indicator = f'{ext_str}.done'
    if os.path.isfile(done_indicator):
        print(f'{ext_str} has beem processed')
        return

    sam_resut_list = io_function.get_file_list_by_ext('.shp',dem_diff_poly_dir,bsub_folder=False)
    classify_res_list = io_function.get_file_list_by_ext('.shp',dem_diff_poly_classify_res_dir,bsub_folder=False)
    if len(sam_resut_list) != len(classify_res_list):
        raise IOError(f"The count for sam results {len(sam_resut_list)} and classifed results {len(classify_res_list)} are different, for {ext_str}")

    # pair each results
    sam_classfiy_res = {}
    for shp in sam_resut_list:
        grid_id = get_grid_id_from_path(shp)
        sam_classfiy_res[grid_id] = [shp]
    for shp in classify_res_list:
        grid_id = get_grid_id_from_path(shp)
        sam_classfiy_res[grid_id].append(shp)

    # combine results
    all_gpds = []
    for idx, g_id in enumerate(sam_classfiy_res.keys()):
        print(datetime.now(), f'({idx+1}/{len(sam_classfiy_res.keys())}) Working on {sam_classfiy_res[g_id][0]}')
        filtered= merge_two_shapefile(sam_classfiy_res[g_id][0], sam_classfiy_res[g_id][1])
        all_gpds.append(filtered)
        # if idx > 10: 
        #     break

    # Merge (concatenate) them:
    merged_gdf = gpd.GeoDataFrame(pd.concat(all_gpds, ignore_index=True))
    save_path = f'{ext_str}_sam_polys_class_1.gpkg'
    print(f'saving to {save_path}')
    merged_gdf.to_file(save_path,driver='GPKG')

    io_function.save_text_to_file(done_indicator,f'completed at {datetime.now()}')


def main():

    org_dir = os.getcwd()
    dem_diff_dir_list = io_function.get_file_list_by_pattern(ArcticDEM_results_dir, 'ext??_*/grid_dem_diffs')

    for dem_diff_dir in dem_diff_dir_list:
        # os.chdir(org_dir)  # change to original dir, because inside sam_segment_a_big_region, they change to other folder

        ext_dir = os.path.dirname(dem_diff_dir)
        work_folder = os.path.basename(ext_dir)
        basic.outputlogMessage(f'processing {ext_dir}')
        dem_diff_poly_dir = os.path.join(ext_dir, 'grid_dem_diffs_sam_results')
        dem_diff_poly_classify_res_dir = os.path.join(ext_dir, 'grid_dem_diffs_reduction_poly_classify_res')

        # select region
        work_folder_begin_str = work_folder[:5]
        copy_classify_result_to_sam_result(work_folder_begin_str, dem_diff_poly_dir, dem_diff_poly_classify_res_dir)


if __name__ == '__main__':
    # test_merge_two_shapefile()

    main()

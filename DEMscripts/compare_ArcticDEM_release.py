#!/usr/bin/env python
# Filename: compare_ArcticDEM_release.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 13 March, 2024
"""

import os,sys
from datetime import datetime
import geopandas as gpd
import pandas as pd

sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS/'))
import vector_gpd
import basic_src.io_function as io_function

def save_ArcticDEM_release_years(release_path,data_dir):

    # save 2019 or newer "acqdate1", "acqdate2"
    gdf = gpd.read_file(release_path)

    # Convert the "acqdate1" field to a pandas DateTime object
    gdf['acqdate1'] = pd.to_datetime(gdf['acqdate1'])

    # Filter the GeoDataFrame based on the condition for "acqdate1"
    filtered_gdf = gdf[gdf['acqdate1'].dt.year >= 2019]
    print('find %d records'%len(filtered_gdf))

    # change to str before saving
    filtered_gdf['acqdate1'] = filtered_gdf['acqdate1'].astype(str)

    # Save the filtered GeoDataFrame to a new shapefile
    output_shapefile_path = io_function.get_name_by_adding_tail(os.path.basename(release_path),'2019toNEW')
    output_shapefile_path = os.path.join(data_dir, output_shapefile_path)
    filtered_gdf.to_file(output_shapefile_path)
    print('saved to %s'%output_shapefile_path)


def compare_ArcticDEM_release(rel_old, rel_new):

    print(datetime.now(), 'comparing two ArcticDEM release:')
    print(rel_old)
    print(rel_new)

    #  pairname is not unique

    # read attribute
    old_pair_names = vector_gpd.read_attribute_values_list(rel_old,'pairname')
    new_pair_names = vector_gpd.read_attribute_values_list(rel_new,'pairname')
    print('count in old_pair_names:', len(old_pair_names), len(set(old_pair_names)))
    print('count in new_pair_names:', len(new_pair_names), len(set(new_pair_names)))

    diff_new_to_old = list(set(new_pair_names) - set(old_pair_names))
    diff_old_to_new = list(set(old_pair_names) - set(new_pair_names))

    print('count in diff_new_to_old:', len(diff_new_to_old))
    print('count in diff_old_to_new:', len(diff_old_to_new))


    # save those in only in new release
    new_gpd = gpd.read_file(rel_new)
    filtered_new_gpd = new_gpd[new_gpd['pairname'].isin(diff_new_to_old)]
    filtered_new_gpd.to_file("diff_new_to_old.shp")

    # save those in only in old release
    old_gpd = gpd.read_file(rel_old)
    filtered_old_gpd = old_gpd[ old_gpd['pairname'].isin(diff_old_to_new)]
    filtered_old_gpd.to_file('diff_old_to_new.shp')


def main():

    data_dir = os.path.expanduser("~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes")

    rel7 = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_Rel7/ArcticDEM_Strip_Index_Rel7.shp')
    rel2022 = os.path.expanduser('~/Data/Arctic/ArcticDEM/BROWSE_SERVER/indexes/ArcticDEM_Strip_Index_s2s041_shp/ArcticDEM_Strip_Index_s2s041.shp')
    # compare_ArcticDEM_release(rel7, rel2022)

    save_ArcticDEM_release_years(rel2022,data_dir)

if __name__ == '__main__':
    main()

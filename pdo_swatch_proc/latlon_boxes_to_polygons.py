#!/usr/bin/env python
# Filename: latlon_boxes_to_polygons.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 29 April, 2022
"""

import os,sys
import csv

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd

import pandas as pd

def latlon_2_polygons(one_nc_info):
    print(one_nc_info)
    # .nc4 filename, top (max) latitude, bottom (min) latitude, left (min) longitude, right (max) longitude
    nc_filename = one_nc_info[0].strip()
    top_lat = float(one_nc_info[1].strip())
    bottom_lat = float(one_nc_info[2].strip())
    left_lon = float(one_nc_info[3].strip())
    right_lon = float(one_nc_info[4].strip())
    bounding_box = (left_lon,bottom_lat,right_lon,top_lat)       # box: (left, bottom, right, top)
    print(nc_filename,bounding_box)
    polygon = vector_gpd.convert_image_bound_to_shapely_polygon(bounding_box)
    return nc_filename, polygon


def main():
    latlon_csv = os.path.expanduser('~/Data/PDO/PDO_statistics_swatchs/LatLonPDOv3.csv')
    nc_name_latlon = []
    with open(latlon_csv) as f_obj:
        reader = csv.reader(f_obj)
        for row in reader:
            if len(row) >=4 and 'PDO' in row[0]:
                nc_name_latlon.append(row)
                print(row)
    print('total %d nc files'%len(nc_name_latlon))

    polygon_latlon = []
    nc_file_name = []
    for idx, nc_info in enumerate(nc_name_latlon):
        fname, polygon = latlon_2_polygons(nc_info)
        polygon_latlon.append(polygon)
        nc_file_name.append(fname)

    # save polygons
    latlon_csv = os.path.expanduser('~/Data/PDO/PDO_statistics_swatchs/swatch_bounding_boxes.shp')
    save_polyons_attributes = {'Polygons':polygon_latlon, 'nc_file':nc_file_name}
    polygon_df = pd.DataFrame(save_polyons_attributes)

    vector_gpd.save_polygons_to_files(polygon_df,'Polygons',{'init' :'epsg:4326'},latlon_csv)


    pass

if __name__ == '__main__':
    main()
    pass
#!/usr/bin/env python
# Filename: exclude_girds.py
"""
introduction: for some *_grid_ids.txt, exclude some grids according to extents (such as Greenland Icesheet)

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 05 August, 2022
"""

import os,sys
import time
sys.path.insert(0, os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS'))
import basic_src.io_function as io_function
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import pandas as pd
import vector_gpd

from dem_common import grid_ids_txt_dir,get_extent_grid_id_txt_done_files
from dem_common import grid_20_shp
from produce_DEM_diff_ArcticDEM import get_grid_20


def exclude_grids_extent(grid_ids_txt, ext_shp):
    # backup grid_ids_txt is not exist
    bak_grid_txt = grid_ids_txt + '_notExcluded'
    if os.path.isfile(bak_grid_txt) is False:
        io_function.copy_file_to_dst(grid_ids_txt,bak_grid_txt)

    exclude_polys = vector_gpd.read_polygons_gpd(ext_shp)

    # read grids and ids
    all_grid_polys, all_ids = vector_gpd.read_polygons_attributes_list(grid_20_shp, 'id')

    #get grid within the polygon, this will create local_grid_id_txt
    grid_polys, grid_ids = get_grid_20(ext_shp, all_grid_polys, all_ids)
    exclude_ids = []

    for exc_poly in exclude_polys:
        for grid_p, grid_i in zip(grid_polys, grid_ids):
            if exc_poly.contains(grid_p):       # only the grid is fully inside the polygon, not touch edge
                exclude_ids.append(grid_i)

    exclude_ids_str = [str(item) for item in exclude_ids]
    org_grid_ids = io_function.read_list_from_txt(grid_ids_txt)

    modified_grid_ids = [ item for item in org_grid_ids if item not in exclude_ids_str ]

    # over write the original file
    io_function.save_list_to_txt(grid_ids_txt, modified_grid_ids)

    # copy to local folder if does not exist
    local_grid_ids_txt = os.path.basename(grid_ids_txt)
    io_function.copy_file_to_dst(grid_ids_txt,local_grid_ids_txt,overwrite=True)

    keep_grids = [grid for grid, id in zip(grid_polys, grid_ids) if str(id) in modified_grid_ids]

    # save the kept grids to file for checking
    dataframe = pd.DataFrame({'Polygon':keep_grids})
    wkt = map_projection.get_raster_or_vector_srs_info_proj4(grid_20_shp)
    save_path = os.path.splitext(os.path.basename(grid_ids_txt))[0] + '_keep.shp'
    vector_gpd.save_polygons_to_files(dataframe,'Polygon',wkt,save_path)



def main():
    grid_txt = os.path.expanduser('~/Bhaltos2/lingcaoHuang/ArcticDEM_tmp_dir/grid_ids_txt/ext06_for_ArcticDEM_proc_grid_ids.txt')
    exclude_ext = os.path.expanduser('~/Data/Arctic/pan_Arctic/extent/green_icesheet_extent/greenland_icesheet_mask_buff100.shp')

    exclude_grids_extent(grid_txt,exclude_ext)


if __name__ == '__main__':
    main()
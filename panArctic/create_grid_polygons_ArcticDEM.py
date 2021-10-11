#!/usr/bin/env python
# Filename: create_grid_polygons_perma 
"""
introduction: based on coverage of ArcticDEM, create 20 by 20 km grid

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 03 March, 2021
"""

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

import pandas as pd

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp')

def main():

    # grid polygons
    grid_20km = os.path.join(shp_dir,'grid_20km.shp')

    ArcticDEM_coverage = os.path.join(shp_dir,'tiles.shp')
    # qtp_grid_50km and qtb_main_perma_area_simp have the same projection
    grid_prj = map_projection.get_raster_or_vector_srs_info_proj4(grid_20km)
    perma_area_prj = map_projection.get_raster_or_vector_srs_info_proj4(ArcticDEM_coverage)
    if grid_prj != perma_area_prj:
        raise ValueError('%s and %s do not have the same projection'%(grid_prj,perma_area_prj))

    grids  = vector_gpd.read_polygons_gpd(grid_20km)
    DEM_areas = vector_gpd.read_polygons_gpd(ArcticDEM_coverage)

    keep_grids = []

    keep_grid_ids = []
    id = 0

    for idx, grid in enumerate(grids):
        print(' processing %dth grid'%idx)


        for dem_area in DEM_areas:
            inte_res = dem_area.intersection(grid)
            if inte_res.is_empty is False:
                if inte_res.area < 100*100: # if it's too small, ignore it
                    continue
                keep_grids.append(grid)
                keep_grid_ids.append(id)
                id += 1

    # save
    save_path = os.path.join(shp_dir,'ArcticDEM_grid_20km.shp')
    save_polyons_attributes = {'id':keep_grid_ids, "Polygons":keep_grids}

    # wkt_string = map_projection.get_raster_or_vector_srs_info_wkt(qtp_main_perma_area_simp)
    wkt_string = map_projection.get_raster_or_vector_srs_info_proj4(grid_20km)
    polygon_df = pd.DataFrame(save_polyons_attributes)
    vector_gpd.save_polygons_to_files(polygon_df, 'Polygons', wkt_string, save_path)


    pass

if __name__ == "__main__":
    main()
#!/usr/bin/env python
# Filename: create_grid_polygons_perma 
"""
introduction: based on coverage of an input extent, create grids that overlap with another group of vectors

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 13 July, 2025
"""

import os,sys
from optparse import OptionParser
from datetime import datetime

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

# import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union
import geopandas as gpd

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp')
overall_coverage = os.path.join(shp_dir,'tiles.shp')


def create_grids_for_overlap_vectors(coverage,input_vector,grid_size_x,grid_size_y,save_path, exclude_id_txt=None):

    # check the projection
    cover_prj = map_projection.get_raster_or_vector_srs_info_proj4(coverage)
    input_prj = map_projection.get_raster_or_vector_srs_info_proj4(input_vector)
    if cover_prj != input_prj:
        raise ValueError('%s and %s do not have the same projection' % (coverage, input_vector))

    print(datetime.now(), f'Start creating grids for overlapping vectors: {input_vector} with coverage: {coverage}')

    # Get the bounds of the study area (minx, miny, maxx, maxy)
    out_minx, out_miny, out_maxx, out_maxy = vector_gpd.get_vector_file_bounding_box(coverage)
    # geometries = vector_gpd.read_polygons_gpd(input_vector)
    input_gdf = gpd.read_file(input_vector)

    # Step 1: Create a regular grid over the bounding box
    grid_cells = []
    grid_ids = []  # Store the unique IDs for each grid cell

    row = 0
    x = out_minx
    while x < out_maxx:
        col = 0
        y = out_miny
        while y < out_maxy:
            # Create a grid cell as a Shapely box
            grid_cell = box(x, y, x + grid_size_x, y + grid_size_y)

            # Assign a unique ID: 5-digit row + 5-digit column
            grid_id = f"{row:05d}{col:05d}"
            grid_cells.append(grid_cell)
            grid_ids.append(grid_id)

            y += grid_size_y
            col += 1
        x += grid_size_x
        row += 1


    print(datetime.now(), f'obtained {len(grid_cells)} grid cells for the coverage area')    

    # Step 2: Convert grid cells to a GeoDataFrame
    grid_gdf = gpd.GeoDataFrame({'geometry': grid_cells, 'RowCol_id': grid_ids, 'grid_id': grid_ids}, crs=cover_prj)


    # Step 3: Perform spatial join to find intersecting grid cells
    intersecting_gdf = gpd.sjoin(grid_gdf, input_gdf, how='inner', predicate='intersects')
    print(datetime.now(), f'after spatial join, {len(intersecting_gdf)} grid cells intersect with the input vector geometries')

    # Remove duplicate grid cells (if any)
    intersecting_gdf = intersecting_gdf.drop_duplicates(subset='RowCol_id')
    print(datetime.now(), f'after removing duplicates, {len(intersecting_gdf)} grid cells remain')

    # # Step 3: Intersect the grid with the input vector geometries
    # input_geometries_union = unary_union(geometries)  # Union of all input geometries
    # grid_gdf['geometry'] = grid_gdf['geometry'].intersection(input_geometries_union)

    # print(datetime.now(), f'after intersection, {len(grid_gdf)} grid cells remain')

    # Remove empty geometries (resulting from non-overlapping grid cells)
    # grid_gdf = grid_gdf[~grid_gdf['geometry'].is_empty]

    # print(datetime.now(), f'{len(intersecting_gdf)} grid cells remain')

    # to exlcude some ids that not
    if exclude_id_txt is not None:
        exclude_ids = io_function.read_list_from_txt(exclude_id_txt)
        # Exclude rows where 'RowCol_id' is in exclude_ids
        intersecting_gdf = intersecting_gdf[~intersecting_gdf['RowCol_id'].isin(exclude_ids)]
        print(datetime.now(), f'after excluding specified grid ids, {len(intersecting_gdf)} grid cells remain')


    # Step 4: Save the resulting GeoDataFrame to the specified path
    save_path = save_path.replace('.gpkg','.shp')
    intersecting_gdf[['geometry', 'RowCol_id']].to_file(save_path, driver='ESRI Shapefile')
    save_grid_id_path = os.path.splitext(io_function.get_name_by_adding_tail(save_path,'grid_id'))[0] + '.txt'
    io_function.save_list_to_txt(save_grid_id_path,intersecting_gdf['RowCol_id'].to_list())
    print(datetime.now(), f'Grids saved to {save_path}')




def main(options, args):

    input_vector = args[0]
    coverage = options.coverage if options.coverage else overall_coverage
    grid_size_x = options.grid_size_x
    grid_size_y = options.grid_size_y

    save_path = options.save_path if options.save_path else io_function.get_name_by_adding_tail(input_vector,f'{grid_size_x}m_grids')

    exclude_grid_ids = options.exclude_grid_ids

    create_grids_for_overlap_vectors(coverage, input_vector, grid_size_x, grid_size_y, save_path, exclude_id_txt=exclude_grid_ids)

if __name__ == "__main__":
    usage = "usage: %prog [options] vector_file "
    parser = OptionParser(usage=usage, version="1.0 2025-07-13")
    parser.description = 'Introduction: generating grids covering input vectors '

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="the path to save the overlap grids")

    parser.add_option("-x", "--grid_size_x",
                  action="store", dest="grid_size_x",type=int, default=20000,
                  help="the default grid size in x direction")

    parser.add_option("-y", "--grid_size_y",
                  action="store", dest="grid_size_y", type=int, default=20000,
                  help="the default grid size in y direction")

    parser.add_option("-c", "--coverage",
                      action="store", dest="coverage",
                      help="the overall coverage of the entire study areas, default is the ArcticDEM coverage")

    parser.add_option("-e", "--exclude_grid_ids",
                      action="store", dest="exclude_grid_ids",
                      help="the grid id lists saved in a txt file to be excluded")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)
#!/usr/bin/env python
# Filename: create_h3_cell_polygons.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 01 September, 2025
"""
import os,sys
from optparse import OptionParser
from datetime import datetime



deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.basic as basic
import basic_src.map_projection as map_projection
import basic_src.io_function as io_function

# import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union
import pandas as pd
import geopandas as gpd

import h3
import shapely
from shapely.geometry import Polygon

from datetime import datetime

import multiprocessing
from multiprocessing import Pool

shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp')
overall_coverage = os.path.join(shp_dir,'tiles.shp')


# Convert Shapely polygon to H3 LatLngPoly
def xy_to_latlng(coords):
    return [(y, x) for x, y in coords]  # shapely: (x=lng, y=lat) -> H3: (lat, lng)

def latlng_to_xy(coords):
    return [(y, x) for x, y in coords]  # H3: (lat, lng) ->  shapely: (x=lng, y=lat)

def obtain_h3_cells_for_a_polygon(polygon, resolution, finest_res=None ):
    """
    Obtain H3 cells (IDs and Shapely polygon boundaries) covering a Shapely polygon in EPSG:4326.

    Parameters:
        polygon: Shapely Polygon in WGS84 (lon/lat). Reproject before calling if needed.
        resolution: H3 resolution (0-15).

    Returns:
        List h3_id and shapely_polygon
    """
    if polygon is None or polygon.is_empty:
        return []

    # Convert exterior ring to [lat, lng] pairs
    exterior = list(polygon.exterior.coords)
    # print('exterior',exterior)
    if len(exterior) >= 2 and exterior[0] == exterior[-1]:
        exterior = exterior[:-1]
    outer_ll = xy_to_latlng(exterior)  # Ensure this function returns [lat, lng] pairs

    # Convert holes (interior rings)
    holes_ll = []
    for ring in polygon.interiors:
        coords = list(ring.coords)
        if len(coords) >= 2 and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) >= 3:  # Only keep valid rings
            holes_ll.append(xy_to_latlng(coords))

    try:
        h3_poly = h3.LatLngPoly(outer_ll, holes_ll)
    except ValueError:
        print('outer_ll:', outer_ll)
        print('holes_ll:', holes_ll)
        raise  # Propagate exception for better debugging

    # to get cells that overlap, not just their centroid was contained
    if finest_res is not None and resolution < finest_res:
        tmp_ids = h3.polygon_to_cells(h3_poly, finest_res)
        cell_ids = [h3.cell_to_parent(item, res=resolution) for item in tmp_ids]
        cell_ids = list(set(cell_ids))  # Remove duplicates
    else:
        cell_ids = h3.polygon_to_cells(h3_poly, resolution)

    if not cell_ids:
        return [], []

    final_c_ids = []
    final_polys = []

    for cid in cell_ids:
        boundary = h3.cell_to_boundary(cid)  # [(lat, lng), ...]
        if not boundary or len(boundary) < 3:
            continue

        boundary = latlng_to_xy(boundary)

        # Ensure closed ring for Shapely
        if boundary[0] != boundary[-1]:
            boundary = boundary + [boundary[0]]


        poly = Polygon(boundary)
        if poly.is_valid and not poly.is_empty:
            final_c_ids.append(cid)
            final_polys.append(poly)
        else:
            print(f'Warning: The cell polygon is empty or invalid ({poly}), with id: {cid}')

    return final_c_ids, final_polys


def process_one_poly(poly,resolution,poly_to_cell_res):
    poly_bound = vector_gpd.convert_bounds_to_polygon(vector_gpd.get_polygon_bounding_box(poly))
    cell_ids, cell_polys = obtain_h3_cells_for_a_polygon(poly_bound, resolution, finest_res=poly_to_cell_res)
    return cell_ids, cell_polys

def obtain_h3_cells_for_overlap_vectors(input_vector, resolution, save_path, exclude_id_txt= None,
                                        poly_to_cell_res = None,buffer_m=None, process_num=1):

    if os.path.isfile(save_path):
        print(datetime.now(), f'{save_path} alreasy exists, skip')
        return

    original_gpd = gpd.read_file(input_vector)
    # print(in_gpd)
    original_crs = original_gpd.crs
    in_gpd = original_gpd.copy()
    if original_gpd.crs != "EPSG:4326":
        if buffer_m is not None:
            in_gpd['geometry'] = in_gpd['geometry'].buffer(buffer_m)
        in_gpd = in_gpd.to_crs("EPSG:4326")
    else:
        in_gpd = original_gpd
    # print(in_gpd)

    all_cell_ids = []
    all_cell_polys = []
    print(datetime.now(),f'Loaded {len(in_gpd)} polygons')
    if process_num == 1:
        for idx, poly in enumerate(in_gpd.geometry.values):
            # poly = poly.buffer(0.000001)
            # print(idx, poly.is_valid)
            cell_ids, cell_polys  = process_one_poly(poly,resolution,poly_to_cell_res)
            all_cell_ids.extend(cell_ids)
            all_cell_polys.extend(cell_polys)
    else:
        theadPool = Pool(process_num)  # multi processes
        parameters_list = [(poly,resolution,poly_to_cell_res) for idx, poly in enumerate(in_gpd.geometry.values)]
        results = theadPool.starmap(process_one_poly, parameters_list)  # need python3
        for res in  results:
            cell_ids, cell_polys = res
            all_cell_ids.extend(cell_ids)
            all_cell_polys.extend(cell_polys)
        theadPool.close()

    print(datetime.now(), f'Obtained {len(all_cell_ids)} cell at res: {resolution}')
    id_colum_name = f"h3_id_{resolution}"

    df = pd.DataFrame({id_colum_name: all_cell_ids, "geometry": all_cell_polys})
    df = df.drop_duplicates(subset=id_colum_name)
    print(datetime.now(), f'After removing duplicates,  {len(df)} cells remains')

    exclude_ids = io_function.read_list_from_txt(exclude_id_txt) if exclude_id_txt is not None else None
    # remove exclude_ids
    if exclude_ids is not None:
        start_len = len(df)
        df = df[~df[id_colum_name].isin(exclude_ids)]
        print(datetime.now(), f'Removed {start_len - len(df)} cells based on exclude list; {len(df)} remain.')

    save_gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
    # print(save_gdf)
    if original_crs != "EPSG:4326":
        save_gdf = save_gdf.to_crs(original_crs)

    # Perform spatial join to find overlapping or touching geometries
    overlap_touch = gpd.sjoin(save_gdf, original_gpd, how='inner', predicate='intersects')
    # Only keep columns from save_gdf
    overlap_touch = overlap_touch[save_gdf.columns]
    # Remove duplicate geometries in the result
    overlap_touch = overlap_touch.drop_duplicates(subset=['geometry'])

    print(datetime.now(), f'After overlap checking with original input vector, kept {len(overlap_touch)} cells')
    save_gdf = overlap_touch

    # print(save_gdf)
    save_gdf.to_file(save_path)
    print(datetime.now(), f'Saved {len(save_gdf)} H3 cells to {save_path}')

    # save cell ids
    cell_ids_txt = os.path.splitext(save_path)[0] + '_cell_ids.txt'
    save_cell_ids = save_gdf[id_colum_name].to_list()
    io_function.save_list_to_txt(cell_ids_txt, save_cell_ids)





def test_obtain_h3_cells_for_overlap_vectors():
    input_vector = os.path.expanduser('~/Data/published_data/'
    'Yili_Yang_etal_ARTS_The-Arctic-Retrogressive-Thaw-Slumps-Data-Set/ARTS_main_dataset_v.3.1.0_ClassInt_org_c1.gpkg')

    input_vector = os.path.expanduser('~/Data/published_data/Dai_etal_2025_largeRTS_ArcticDEM/ArcticRTS_epsg3413_ClassInt.shp')

    resolution = 7
    save_path = f'h3_res{resolution}_cells_set01.shp'
    obtain_h3_cells_for_overlap_vectors(input_vector, resolution, save_path, exclude_id_txt=None,poly_to_cell_res=13)


def main(options, args):

    input_vector = args[0]
    save_path = options.save_path
    h3_resolution = options.h3_resolution
    exclude_grid_ids_txt = options.exclude_grid_ids
    poly_to_cell_res = options.poly_to_cell_res
    buffer_meters = options.buffer_meters
    process_num = options.process_num
    if save_path is None:
        save_path = f'h3_cells_res{h3_resolution}_{io_function.get_name_no_ext(input_vector)}.shp'

    obtain_h3_cells_for_overlap_vectors(input_vector, h3_resolution, save_path, exclude_id_txt=exclude_grid_ids_txt,
                                        poly_to_cell_res=poly_to_cell_res, buffer_m=buffer_meters, process_num=process_num)


if __name__ == "__main__":
    usage = "usage: %prog [options] vector_file "
    parser = OptionParser(usage=usage, version="1.0 2025-07-13")
    parser.description = 'Introduction: generating h3 cells covering input vectors '

    # test_obtain_h3_cells_for_overlap_vectors()
    # sys.exit(0)

    parser.add_option("-s", "--save_path",
                      action="store", dest="save_path",
                      help="the path to save the overlap h3 cells")

    parser.add_option("-r", "--h3_resolution",
                  action="store", dest="h3_resolution",type=int, default=8,
                  help="the resolution of the h3 geo index")

    parser.add_option("-f", "--poly_to_cell_res",
                  action="store", dest="poly_to_cell_res",type=int,
                  help="the resolution for polygon_to_cells, can be larger than h3_resolution. As "
                       "polygon_to_cells only check centroid of cell contained in a polygons, "
                       "set this to draw cells in finer resolution ")

    parser.add_option("-b", "--buffer_meters",
                  action="store", dest="buffer_meters",type=float,
                  help="buffer the input vectors before getting cells, makint it also works for points and lines")

    parser.add_option("-p", "--process_num",
                      action="store", dest="process_num", type=int, default=16,
                      help="number of processes to get h3 cells")


    parser.add_option("-e", "--exclude_grid_ids",
                      action="store", dest="exclude_grid_ids",
                      help="the cell id lists saved in a txt file to be excluded")


    (options, args) = parser.parse_args()
    if len(sys.argv) < 2 or len(args) < 1:
        parser.print_help()
        sys.exit(2)

    main(options, args)

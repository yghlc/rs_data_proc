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
from shapely.geometry import shape, mapping, Polygon, MultiPolygon
from shapely.geometry.polygon import orient
import h3pandas  # registers the .h3 accessor
from typing import Optional, Union, List, Tuple, Iterable, Set


shp_dir = os.path.expanduser('~/Data/Arctic/ArcticDEM/grid_shp')
overall_coverage = os.path.join(shp_dir,'tiles.shp')


# Convert Shapely polygon to H3 LatLngPoly
def xy_to_latlng(coords):
    return [(y, x) for x, y in coords]  # shapely: (x=lng, y=lat) -> H3: (lat, lng)

def obtain_h3_cells_for_a_polygon(polygon, resolution):
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

    finest_res = 14
    # to get cells that overlap, not just their centroid was contained
    if resolution < finest_res:
        tmp_ids = h3.polygon_to_cells(h3_poly, finest_res)
        cell_ids = [h3.cell_to_parent(item, res=resolution) for item in tmp_ids]
        cell_ids = list(set(cell_ids))  # Remove duplicates
    else:
        cell_ids = h3.polygon_to_cells(h3_poly, resolution)

    if not cell_ids:
        return []

    final_c_ids = []
    final_polys = []

    for cid in cell_ids:
        boundary = h3.cell_to_boundary(cid)  # [[lng, lat], ...]
        if not boundary or len(boundary) < 3:
            continue
        # Ensure closed ring for Shapely
        if boundary[0] != boundary[-1]:
            boundary = boundary + tuple([boundary[0]])
        poly = Polygon(boundary)
        if poly.is_valid and not poly.is_empty:
            final_c_ids.append(cid)
            final_polys.append(poly)
        else:
            print(f'Warning: The cell polygon is empty or invalid ({poly}), with id: {cid}')

    return final_c_ids, final_polys



def obtain_h3_cells_for_overlap_vectors(input_vector, resolution, save_path, exclude_id_txt= None):

    in_gpd = gpd.read_file(input_vector)
    # print(in_gpd)
    original_crs = in_gpd.crs
    if in_gpd.crs != "EPSG:4326":
        in_gpd = in_gpd.to_crs("EPSG:4326")
    # print(in_gpd)

    all_cell_ids = []
    all_cell_polys = []
    print(f'Loaded {len(in_gpd)} polygons')
    for idx, poly in enumerate(in_gpd.geometry.values):
        # poly = poly.buffer(0.000001)
        # print(idx, poly.is_valid)

        poly_bound = vector_gpd.convert_bounds_to_polygon(vector_gpd.get_polygon_bounding_box(poly))
        cell_ids, cell_polys  = obtain_h3_cells_for_a_polygon(poly_bound, resolution)
        all_cell_ids.extend(cell_ids)
        all_cell_polys.extend(cell_polys)

    print(f'Obtained {len(all_cell_ids)} cell at res: {resolution}')
    id_colum_name = f"h3_id_{resolution}"

    df = pd.DataFrame({id_colum_name: all_cell_ids, "geometry": all_cell_polys})
    df = df.drop_duplicates(subset=id_colum_name)
    print(f'After removing duplicates,  {len(df)} cells remains')

    exclude_ids = io_function.read_list_from_txt(exclude_id_txt) if exclude_id_txt is not None else None
    # remove exclude_ids
    if exclude_ids is not None:
        pass


    save_gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
    if original_crs != "EPSG:4326":
        save_gdf = save_gdf.to_crs(original_crs)
    save_gdf.to_file(save_path)




def test_obtain_h3_cells_for_overlap_vectors():
    input_vector = os.path.expanduser('~/Data/published_data/'
    'Yili_Yang_etal_ARTS_The-Arctic-Retrogressive-Thaw-Slumps-Data-Set/ARTS_main_dataset_v.3.1.0_ClassInt_org_c1.gpkg')

    input_vector = os.path.expanduser('~/Data/published_data/Dai_etal_2025_largeRTS_ArcticDEM/ArcticRTS_epsg3413_ClassInt.shp')

    resolution = 8
    save_path = f'h3_res{resolution}_cells_set01.shp'
    obtain_h3_cells_for_overlap_vectors(input_vector, resolution, save_path, exclude_id_txt=None)


def main(options, args):

    pass


if __name__ == "__main__":
    usage = "usage: %prog [options] vector_file "
    parser = OptionParser(usage=usage, version="1.0 2025-07-13")
    parser.description = 'Introduction: generating grids covering input vectors '

    test_obtain_h3_cells_for_overlap_vectors()
    sys.exit(0)

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

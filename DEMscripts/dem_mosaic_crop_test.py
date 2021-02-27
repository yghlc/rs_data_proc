#!/usr/bin/env python

import os,sys

deeplabforRS =  os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, deeplabforRS)
import vector_gpd
import basic_src.io_function as io_function
import basic_src.map_projection as map_projection

import dem_mosaic_crop
import pandas as pd

def test_get_dem_tif_ext_polygons():
    work_dir = os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff/dem_tifs')
    os.chdir(work_dir)

    tifs = io_function.get_file_list_by_ext('.tif', work_dir,bsub_folder=False)
    polygons = dem_mosaic_crop.get_dem_tif_ext_polygons(tifs)

    data = {'poly':polygons}
    pddata = pd.DataFrame(data)
    wkt_str = map_projection.get_raster_or_vector_srs_info_wkt(tifs[0])
    save_path = 'tif_extent.shp'
    vector_gpd.save_polygons_to_files(pddata,'poly',wkt_str,save_path)




if __name__ == '__main__':
    pass
#!/usr/bin/env python
# Filename: get_greenland_icesheet_bounary.py 
"""
introduction:  try to get the boundary of GreenLand Ice Sheet from a data product:
"Surface Elevation Change" from https://climate.esa.int/en/projects/ice-sheets-greenland/data/
after downloading: cci_sec_2021.zip, then unpacking

# gdal_translate NETCDF:CCI_GrIS_RA_SEC_5km_Vers3.0_2021-08-09.nc:SEC sec.tif : the width and high is not correct.
# QGIS also cannot ready the file correctly

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 05 August, 2022
"""

import os,sys
import netCDF4 as nc
import numpy as np

codes_dir2 = os.path.expanduser('~/codes/PycharmProjects/DeeplabforRS')
sys.path.insert(0, codes_dir2)
import raster_io
import vector_gpd

import rasterio
from rasterio.transform import from_origin

def main():
    nc_file = os.path.expanduser('~/Downloads/Release/CCI_GrIS_IS2_SEC_5km_Vers1.0_2021-08-09.nc')
    ds = nc.Dataset(nc_file, 'r', format='NETCDF4')
    print('dataset: ')
    print(ds)
    # print('dimensions:')
    # for item in ds.dimensions.values():
    #     print(item)
    # print('variable:')
    # for var in ds.variables.values():
    #     print(var)
    #     print('\n')

    # lat = ds.variables['lat'][:]
    # print(lat.shape)
    # print(lat)
    # lon = ds.variables['lon'][:]
    # print(lon.shape)
    # print(lon)


    SEC = ds.variables['SEC'][:]
    print(SEC.shape)
    # print(SEC)
    SEC = np.squeeze(SEC)   # Remove dimensions of size 1 from ndarray
    print(SEC.shape)
    print(SEC)

    height = ds.grid_ny
    width = ds.grid_nx
    minx = ds.grid_minx
    miny = ds.grid_miny
    res = float(ds.spatial_resolution.split()[0])
    maxy = miny + height * res

    prjection = ds.grid_projection



    # raster_io.save_numpy_array_to_rasterfile()
    # transform = from_origin(minx, miny, res, -res) # (west, north, xsize, ysize):
    # new_dataset = rasterio.open('mask.tif', 'w', driver='GTiff',
    #                             height=SEC.shape[0], width=SEC.shape[1],
    #                             count=1, dtype=str(SEC.dtype),
    #                             crs=prjection,
    #                             transform=transform)
    # new_dataset.write(SEC, 1)
    # new_dataset.close()

    mask = np.zeros(SEC.shape,dtype=np.uint8)
    mask[np.where(np.isnan(SEC) == False)] = 1
    transform = from_origin(minx, miny, res, -res) # (west, north, xsize, ysize):
    new_dataset = rasterio.open('mask.tif', 'w', driver='GTiff',
                                height=mask.shape[0], width=mask.shape[1],
                                count=1, dtype=str(mask.dtype),
                                crs=prjection,
                                nodata=0,
                                transform=transform)
    new_dataset.write(mask, 1)
    new_dataset.close()


    # to shapefile
    vector_gpd.raster2shapefile('mask.tif',out_shp='greenland_icesheet_mask.shp')
    vector_gpd.fill_holes_in_polygons_shp('greenland_icesheet_mask.shp','greenland_icesheet_mask.shp')

    # then buffer and remove small polygons in QGIS




if __name__ == '__main__':
    main()
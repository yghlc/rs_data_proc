#!/bin/bash

## Introduction:  convert grid polygons to raster


# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace


shp_dir=~/Data/Arctic/ArcticDEM/grid_shp
layer=ArcticDEM_grid_20km
grid_20km=${shp_dir}/${layer}.shp

res=20000
save_raster=grid_20km_bin.tif

rm ${save_raster} | true

# it turns out, there is a half pixel offset between grid and pixels,
# so dont use this gdal_rasterize, but using rasterio, see convert_grid_polygons_to_raster.py

# at all touch:  -at
gdal_rasterize -burn 1 -l ${layer} -tr ${res} ${res}  -ot Byte ${grid_20km} ${save_raster}



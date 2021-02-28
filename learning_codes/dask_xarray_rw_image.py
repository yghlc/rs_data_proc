#!/usr/bin/env python
# Filename: dask_xarray_rw_image.py 
"""
introduction: https://examples.dask.org/applications/satellite-imagery-geotiff.html

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 28 February, 2021
"""


import os
import json
import rasterio

import xarray as xr

work_dir=os.path.expanduser('~/Data/Arctic/canada_arctic/DEM/WR_dem_diff')
img_path = os.path.join(work_dir,'WR_extent_DEM_diff_sub_1_reg.tif')


def main():
    img = rasterio.open(img_path)
    print(img.is_tiled)
    print(img.block_shapes)

    # have not finished.


if __name__ == '__main__':
    main()
    pass


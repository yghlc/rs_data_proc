#!/bin/bash

## Introduction:  create a virtual mosaic for many rasters in different grids, then visualize it in QGIS.
## in QGIS, open many rasters at once are really slow

# if there are many files, it's relatively slow, but still manageable
# build overview files may can help

# ref: https://gis.stackexchange.com/questions/52367/handling-many-raster-files-in-qgis

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 10 September, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

dir=~/Data/dem_processing/dem_hillshade_newest_HWLine_grid
#ls ${dir}/*.tif | grep -v count > infile_list.txt
#save_vrt=~/Data/dem_processing/dem_hillshade_newest_HWLine_grid.vrt
#gdalbuildvrt -resolution average -r nearest -input_file_list infile_list.txt ${save_vrt}
#rm infile_list.txt

ls ${dir}/*count.tif > infile_list.txt
save_vrt=~/Data/dem_processing/dem_hillshade_newest_HWLine_count_grid.vrt
gdalbuildvrt -resolution average -r nearest -input_file_list infile_list.txt ${save_vrt}
rm infile_list.txt


# merge shapefile (GPKG format can support > 2GB)
#ogrmerge.py -o headwall_line_shps.gpkg -f GPKG -single -progress  headwall_shps_grid*/*.shp
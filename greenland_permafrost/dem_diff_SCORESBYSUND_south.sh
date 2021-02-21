#!/usr/bin/env bash

## Introduction:  do dem difference for KANGERLUSSUAQ

## run in: ~/Data/Greenland_permafrost/2021_NNA_PROJECT/KANGERLUSSUAQ/sub_dem

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 21 February, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

ext_shp=~/Data/Greenland_permafrost/2021_NNA_PROJECT/SCORESBYSUND_LANDSLIDES/extent/south_ext.shp

# Strip DEM shapefile
#dem_list=../dem_list.txt
dem_list=dem_south_list.txt

dem_dir=./
save_dir=./
format=GTiff


py=~/codes/PycharmProjects/rs_data_proc/tools/ArcticDEM_proc_grid.py

#--remove_inter_data
${py} ${ext_shp} -a ${dem_dir} -d ${save_dir} -f ${format} -l ${dem_list} \
        --create_mosaic_id --create_mosaic_date --create_dem_diff --remove_inter_data
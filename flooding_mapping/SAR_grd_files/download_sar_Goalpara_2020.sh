#!/bin/bash

## Introduction:  reprojection  and building overview images


# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

py=~/codes/PycharmProjects/rs_data_proc/flooding_mapping/SAR_grd_files/sar_grd_files_gee.py

res=10
start_date=2020-07-01
end_date=2020-07-06
ext_shp=~/Data/flooding_area/Goalpara/extent/Goalpara_ext.shp
region_name=Goalpara
save_dir=./

${py} -s ${start_date} -e ${end_date} --extent_shp=${ext_shp} -n ${region_name} -r ${res} ${save_dir}



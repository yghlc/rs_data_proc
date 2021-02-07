#!/usr/bin/env bash

## Introduction: create daily mosaic of Planet images

# run this script in ~/Data/Arctic/canada_arctic/rsImages

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 6 February, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace


out_res=3

# create the four bands, with 16 bit
#sr_max=2000
#sr_min=0

cloud_cover=0.3
# due to memory issues, we should not set this too high
process_num=4

py=~/codes/PycharmProjects/rs_data_proc/tools/mosaic_images_crop_grid.py
img_dir=~/Data/Arctic/canada_arctic/rsImages/planet_sr_images/2020_July_August

function daily_mosaic () {
    grid_shp=$1
    region=$2

    work_dir=${region}_daily_mosaic
    mkdir -p ${work_dir}
    cd ${work_dir}

    # --rgb  -u ${sr_max} -l ${sr_min}
    ${py} ${img_dir} ${grid_shp} -r ${out_res} -c ${cloud_cover} \
    -p ${process_num} --group_date

    # remove tmp folders
#    rm -r RGB_images_* || true
    rm -r planet_images_reproj_* || true

    cd ..
}

region=WR
grid_shp=~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent.shp
daily_mosaic  ${grid_shp} ${region}

region=Banks_east
grid_shp=~/Data/Arctic/canada_arctic/Banks_east/extent/Banks_east_extent.shp
daily_mosaic  ${grid_shp} ${region}


region=Ellesmere_Island
grid_shp=~/Data/Arctic/canada_arctic/Ellesmere_Island/extent/HotWC_extent.shp
daily_mosaic  ${grid_shp} ${region}




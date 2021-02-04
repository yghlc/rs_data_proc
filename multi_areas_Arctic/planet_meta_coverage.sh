#!/usr/bin/env bash

## Introduction: get the metadata and coverage of  downloaded Planet images

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 3 February, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

work_dir=~/Data/Arctic/canada_arctic/rsImages
cd ${work_dir}

#work_dir=${PWD}
code_dir=~/codes/PycharmProjects/rs_data_proc



cloud_cover_thr=0.1

py=${code_dir}/planetScripts/get_planet_image_list.py

function meta () {
   extent_file=$1
   region=$2
   meta_dir=${region}_Planet_meta
   mkdir -p ${meta_dir}
   save_shp=${meta_dir}/${region}_Planet_meta.shp

    # --group_date
   ${py} ${work_dir}/planet_sr_images -e ${extent_file} -m ${save_shp}
}

region=WR
ext_shp=~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent_latlon.shp
meta  ${ext_shp} ${region}


region=Banks_east
ext_shp=~/Data/Arctic/canada_arctic/Banks_east/extent/Banks_east_extent_latlon.shp
meta  ${ext_shp} ${region}


region=Ellesmere_Island
ext_shp=~/Data/Arctic/canada_arctic/Ellesmere_Island/extent/HotWC_extent_latlon.shp
meta  ${ext_shp} ${region}


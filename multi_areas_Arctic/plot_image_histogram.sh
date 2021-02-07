#!/usr/bin/env bash

## Introduction: Plot the histogram of images

# run this script in ~/Data/Arctic/canada_arctic/rsImages

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 6 February, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace


py=~/codes/PycharmProjects/rs_data_proc/tools/plot_Images_histogram.py

img_dir=~/Data/Arctic/canada_arctic/rsImages/planet_sr_images/2020_July_August

function histogram(){
    extent_file=$1
    ${py} ${img_dir} --planet_geojson -e ${extent_file}
}

#ext_shp=~/Data/Arctic/canada_arctic/Willow_River/extent/WR_extent_latlon.shp
#histogram ${ext_shp}
#
#ext_shp=~/Data/Arctic/canada_arctic/Banks_east/extent/Banks_east_extent_latlon.shp
#histogram ${ext_shp}
#
#ext_shp=~/Data/Arctic/canada_arctic/Ellesmere_Island/extent/HotWC_extent_latlon.shp
#histogram ${ext_shp}


# plot histogram of mosaic
function histogram_mosaic(){
    tif_dir=$1
    up_per=$2
    low_per=$3
    ${py} ${tif_dir} -u ${up_per} -l ${low_per}
}
dir=~/Data/Arctic/canada_arctic/rsImages
histogram_mosaic ${dir}/WR_daily_mosaic
histogram_mosaic ${dir}/Banks_east_daily_mosaic
histogram_mosaic ${dir}/Ellesmere_Island_daily_mosaic

# histogram for all
histogram_mosaic ${dir}




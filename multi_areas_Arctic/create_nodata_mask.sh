#!/bin/bash

## Introduction:  based on the images for generating trainig polygons, to create nodata mask and will be used to mask other daily mosaic

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 29 Oct, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

function nodata_mask() {
  # the input images are 16 bit
  img=$1
  nodata=$2
  output=$3
  # create a mask
  gdal_calc.py --calc="A!=${nodata}" --outfile=tmp.tif -A ${img} --NoDataValue 0 --type=Byte

  # remove nodata
  gdal_edit.py -unsetnodata tmp.tif

  # compress
  gdal_translate -co "compress=lzw" tmp.tif $output

  rm tmp.tif

}

# Willow Riveer
#img=~/Data/Arctic/canada_arctic/Willow_River/Planet2020/20200818_mosaic.tif
#output=~/Data/Arctic/canada_arctic/Willow_River/Planet2020/20200818_mosaic_nodataMask.tif
#nodata_mask ${img} 65535 ${output}

# Banks_east
#img=~/Data/Arctic/canada_arctic/Banks_east/Banks_Island_mosaic.tif
#output=~/Data/Arctic/canada_arctic/Banks_east/Banks_Island_mosaic_nodataMask.tif
#nodata_mask ${img} 0 ${output}


# Ellesmere_Island
img=~/Data/Arctic/canada_arctic/Ellesmere_Island/HotWeatherCreek_Mosaic.tif
output=~/Data/Arctic/canada_arctic/Ellesmere_Island/HotWeatherCreek_Mosaic_nodataMask.tif
nodata_mask ${img} 0 ${output}

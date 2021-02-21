#!/usr/bin/env bash

## Introduction:  convert to 8bit

## run in: ~/Data/Greenland_permafrost/2021_NNA_PROJECT/ITTOQQOTROORMIIT
# or
# ~/Data/Greenland_permafrost/2021_NNA_PROJECT/KANGERLUSSUAQ


#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 19 February, 2021

out_dir=sub_region_8bit
rm -r ${out_dir}
mkdir ${out_dir}
for tif in sub_region/*_ortho.tif; do

    echo $tif
    filename=$(basename "$tif")
#	  extension="${filename##*.}"
	  filename_noext="${filename%.*}"

    out8bit=${out_dir}/${filename_noext}_8bit.tif
    gdal_contrast_stretch -percentile-range 0.02 0.98 ${tif} ${out8bit}
#    gdal_translate -ot Byte -scale 100 500 1 255 ${tif} ${out8bit}

done
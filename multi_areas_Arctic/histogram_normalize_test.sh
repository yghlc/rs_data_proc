#!/usr/bin/env bash

## Introduction: Apply histogram normalization

# run this script in ~/Data/Arctic/canada_arctic/rsImages

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 16 June, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

## set LD_LIBRARY_PATH for gdal_contrast_stretch on server
#export LD_LIBRARY_PATH=~/programs/miniconda3/lib:$LD_LIBRARY_PATH
#export PATH=~/programs/dans-gdal-scripts_install/bin:$PATH

#gdal_contrast_stretch -histeq <target_stddev>  src.tif dst.tif   Histogram normalize to a target bell curve

plot_py=~/codes/PycharmProjects/rs_data_proc/tools/plot_Images_histogram.py

dir=~/Data/Arctic/canada_arctic/Willow_River/Planet2020


# 16 bit src
src=${dir}/20200818_mosaic.tif  # convert from 16 bit to 8bit
${plot_py} ${src} --value_range_min=1 --value_range_max=6000 -b 500

# 8bit src
#src=${dir}/20200818_mosaic_8bit_rgb.tif  # convert from 8bit to 8bit
#${plot_py} ${src} --value_range_min=1 --value_range_max=255 -b 254


#histnormalize,
#for std in 20 40 60 80 100 120 140 180 200; do
for std in 20 30 40 50 60 70 80 100; do
  echo $std
  dst=$(echo ${src}| cut -d '.' -f 1)
  dst=${dst}_histeq_${std}.tif
  gdal_contrast_stretch -histeq $std ${src} $dst

  # plot histogram
  ${plot_py} ${dst} --value_range_min=1 --value_range_max=255 -b 254
done








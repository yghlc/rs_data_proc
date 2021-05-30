#!/bin/bash

## Introduction:  reprojection  and building overview images


# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

py8=~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py

dir=~/Bhaltos2/lingcaoHuang/flooding_area/Houston/Houston_SAR_GRD_gee/S1_Houston
input=$1
if [ ! -z "$input" ]; then
  dir=$input
fi
echo "Input dir is" $dir

src_min=-20
src_max=1
dst_min=1
dst_max=255
nodata=0

for tif in ${dir}/*.tif; do
  echo $tif
  filename=$(basename $tif)
  filename_noext=${filename%.*}

  # to 8bit
  out_8bit=${filename_noext}_8bit.tif
  $py8  -s ${src_min} ${src_max} ${dst_min} ${dst_max} -n ${nodata} ${tif} ${out_8bit}

done








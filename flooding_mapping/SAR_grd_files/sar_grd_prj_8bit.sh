#!/bin/bash

## Introduction:  reprojection  and building overview images


# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

py8=~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py

res=10
# US National Atlas Equal Area, do we need to reproject?
prj=EPSG:2163

dir=~/Bhaltos2/lingcaoHuang/flooding_area/Houston/Houston_pre_processed_power_transform
input=$1
if [ ! -z "$input" ]; then
  dir=$input
fi
echo "Input dir is" $dir

src_min=0.60
src_max=0.75
dst_min=1
dst_max=255
nodata=0

tmp=/tmp

for tif in ${dir}/*/*.tif; do
  echo $tif
  filename=$(basename $tif)
  filename_noext=${filename%.*}

  # reporject
  prj_out=${tmp}/${filename_noext}_prj.tif
  gdalwarp -tr ${res} ${res} -t_srs ${prj} \
          -multi -wo NUM_THREADS=8  -r cubic -dstnodata ${nodata} $tif  ${prj_out}

  # to 8bit
  out_8bit=${filename_noext}_prj_8bit.tif
  $py8  -s ${src_min} ${src_max} ${dst_min} ${dst_max} -n ${nodata} ${prj_out} ${out_8bit}

  rm ${prj_out}

done








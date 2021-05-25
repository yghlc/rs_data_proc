#!/bin/bash

## Introduction:  reprojection  and building overview images


# Exit immediately if a command exits with a non-zero status. E: error trace
#set -eE -o functrace

# unpack to this folder first
dst=/home/lhuang/Data/tmp

res=0.5
# US National Atlas Equal Area
prj=EPSG:2163

for dd in $(ls -d ${dst}/*_RGB); do
  echo $dd
  folder=$(basename $dd)
  mkdir ${folder}
  for tif in ${dd}/*.tif; do
    echo $tif
    name=$(basename $tif)
    out=${folder}/${name}
    gdalwarp -co compress=lzw -co tiled=yes -co bigtiff=if_safer -tr ${res} ${res} -t_srs ${prj} \
          -multi -wo NUM_THREADS=8  -r cubic $tif  ${out}

    # add overview
    gdaladdo -ro -r cubic $out 4 8 16 32 64

  done

done

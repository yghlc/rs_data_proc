#!/bin/bash

# need gdal, on may laptop, run: conda activate base

dem_dir=~/Data/PDO/PDO_statistics_swatchs/ArcticDEM_eachSwatch
re_prj_dir=~/Data/PDO/PDO_statistics_swatchs/ArcticDEM_eachSwatch_latlon

mkdir -p ${re_prj_dir}

prj=EPSG:4326

for tif in ${dem_dir}/*crop.tif; do
  echo $tif
  filename=$(basename ${tif})
  filename_noext=${filename%.*}

  # reporject
  prj_out=${re_prj_dir}/${filename_noext}_prj.tif
  gdalwarp -t_srs ${prj} \
          -multi -wo NUM_THREADS=8 -r near -co "compress=lzw"  $tif  ${prj_out}

done



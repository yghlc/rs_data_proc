#!/bin/bash

# introduction: select a few tiles and reproject

# run this script in Bhaltos2/lingcaoHuang/global_surface_water

#authors: Huang Lingcao
#email:huanglingcao@gmail.com
#add time: 3 Sep, 2021

# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

# reproject to Polar Stereographic North
t_prj=EPSG:3413
out_dir=extent_epsg3413

mkdir -p ${out_dir}

# 170W_80N, 170W_70N, 160W_80N,160W_70N,150W_80N,150W_70N

for lat in 60N 70N 80N ;do
  for tif in extent/*${lat}*; do
    echo $tif

    # re-project
    basename=$(basename $tif)
    output=${out_dir}/${basename}

    if [ -f ${output} ]; then
      echo ${output} already exist, skip
      continue
    fi

    echo $output
    gdalwarp -r near -t_srs ${t_prj} -tr 30 30 -co compress=lzw $tif ${output}

  done
done






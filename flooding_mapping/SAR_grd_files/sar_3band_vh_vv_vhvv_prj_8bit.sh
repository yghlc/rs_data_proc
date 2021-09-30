#!/bin/bash

## Introduction:  reprojection  and building overview images


# Exit immediately if a command exits with a non-zero status. E: error trace
set -eE -o functrace

py8=~/codes/PycharmProjects/rs_data_proc/tools/convertTo8bit.py

res=10
# US National Atlas Equal Area, do we need to reproject?
prj=EPSG:2163

dir=~/Bhaltos2/lingcaoHuang/flooding_area/Houston/compose_3bands_VH_VV_VHVV
input=$1
if [ ! -z "$input" ]; then
  dir=$input
fi
echo "Input dir is" $dir

src_min_b1=0.70
src_max_b1=1.0

src_min_b2=0.74
src_max_b2=1.0

src_min_b3=1.0
src_max_b3=2.0

dst_min=1
dst_max=255

# set nodata  in the source images
# trouble: all the values in 3rd is 1, nodata region (outside image) is also 1, so sad.
# so we try to set pixel with band1+band2+band3 >=3 as 0, but may also mask some valid pixel.

src_nodata=0
nodata=0

# 98% cut histogram, min, max information
# 20170829_RGB_composite.tif
#0.692159, 1
#0.740218, 1
#1, 1

# 20170829_RGB_composite_north.tif
#0.73266, 1
#0.788721, 1
#1,1


tmp=/tmp

for tif in ${dir}/*.tif; do
  echo $tif
  filename=$(basename $tif)
  filename_noext=${filename%.*}

  out_8bit=${filename_noext}_prj_8bit.tif
  if [ -f ${out_8bit} ]; then
    echo ${out_8bit} exist, skip
    continue
  fi

  # mask nodata pixel
  mask_tif=${filename_noext}_mask.tif
  if [ -f ${mask_tif} ]; then
    echo ${mask_tif} exists, skip
  else
    gdal_calc.py -A $tif --A_band=1 -B $tif --B_band=2 -C $tif --C_band=3  -D $tif --allBands=D \
        --outfile=${mask_tif} --NoDataValue=${src_nodata}  --calc="((A+B+C)<3)*D"
  fi

  # reporject
  prj_out=${tmp}/${filename_noext}_prj.tif
  if [ -f ${prj_out} ]; then
    echo ${prj_out} exist, skip
  else
    gdalwarp -tr ${res} ${res} -t_srs ${prj} \
          -multi -wo NUM_THREADS=8  -r cubic -dstnodata ${src_nodata} $mask_tif  ${prj_out}
  fi

  # to 8bit

  $py8  -s ${src_min_b1} ${src_max_b1} ${dst_min} ${dst_max}  \
          -s ${src_min_b2} ${src_max_b2} ${dst_min} ${dst_max} \
          -s ${src_min_b3} ${src_max_b3} ${dst_min} ${dst_max} \
          -N ${src_nodata} -n ${nodata} ${prj_out} ${out_8bit}

  rm ${prj_out}

done




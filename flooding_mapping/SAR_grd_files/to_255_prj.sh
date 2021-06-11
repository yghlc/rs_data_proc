#!/bin/bash

res=10
# US National Atlas Equal Area, do we need to reproject?
prj=EPSG:2163
tmp=/tmp
nodata=0

for tif in $(ls ../Harvey_8-29/S1*/*.tif); do
  echo $tif
	

  filename=$(basename $tif)
  filename_noext=${filename%.*}

  # reporject
  prj_out=${tmp}/${filename_noext}_prj.tif
  gdalwarp -tr ${res} ${res} -t_srs ${prj} \
          -multi -wo NUM_THREADS=8  -r near -dstnodata ${nodata} $tif  ${prj_out}

  out=${filename_noext}_prj_255.tif
  gdal_calc.py -A ${prj_out} --calc="A*255" --outfile=${out} --NoDataValue=${nodata}
	#exit

done

#!/bin/bash

pdo_v3_dir=~/Data/PDO/pdo_v3_data/ABoVE_ReSALT_InSAR_PolSAR_V3/data

# convert netcdf to tif
# https://help.marine.copernicus.eu/en/articles/5029956-how-to-convert-netcdf-to-geotiff

# need gdal, on may laptop, run: conda activate base

#
for nc4 in ${pdo_v3_dir}/*.nc4; do
  echo $nc4
  filename=$(basename ${nc4})
  filename_noext=${filename%.*}
  # too many nodata inside the frame for alt
  #gdal_translate NETCDF:${nc4}:alt ${pdo_v3_dir}/${filename_noext}_alt.tif
  gdal_translate NETCDF:${nc4}:qa ${pdo_v3_dir}/${filename_noext}_qa.tif
done


# get shapefile for each tif
py=~/codes/PycharmProjects/rs_data_proc/tools/get_image_valid_extent.py
extent_dir=~/Data/PDO/extent_each_swatch
for tif in ${pdo_v3_dir}/*_qa.tif; do
  echo $tif
  filename=$(basename ${tif})
  filename_noext=${filename%.*}
  ${py} ${tif}  -o ${extent_dir}/${filename_noext}_exent.shp
done
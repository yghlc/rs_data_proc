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
  gdal_translate NETCDF:${nc4}:alt ${pdo_v3_dir}/${filename_noext}_alt.tif
  gdal_translate NETCDF:${nc4}:alt_unc ${pdo_v3_dir}/${filename_noext}_alt_unc.tif
done

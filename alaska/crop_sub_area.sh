#!/bin/bash

shpdir=/scratch/summit/lihu9680/Arctic/canada_arctic/Alaska_west

py=~/codes/PycharmProjects/Landuse_DL/sentinelScripts/get_subImages.py

shp1=${shpdir}/Area1outline/Area1outline_prj.shp
img_dir=/scratch/summit/lihu9680/Arctic/alaska/northern_alaska_rgb_2020_Jul_Aug/northern_alaska_rgb_2020_Jul_Aug_mosaic_3.0

buffersize=300
file_pattern=*.tif
out_dir=./


#${py} --rectangle --no_label_image -b ${buffersize} -e ${file_pattern} ${shp1} ${img_dir} 
for i in 2 3 4 ; do
	echo $i
	shp=${shpdir}/Area${i}outline/Area${i}outline_prj.shp
${py} --rectangle --no_label_image -b ${buffersize} -e ${file_pattern} ${shp} ${img_dir} 

	mv subImages subImages_${i} 

done



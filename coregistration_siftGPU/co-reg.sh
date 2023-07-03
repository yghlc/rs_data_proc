#!/bin/bash
set -eE -o functrace

# regions: BanksEast hotWC Tuk WillowRiver


py=~/codes/PycharmProjects/rs_data_proc/coregistration_siftGPU/coregistration_siftGPU.py
dataDir=${PWD}
#ref=${dataDir}/BanksEast_PlanetScope_20220823_psscene_analytic_sr_udm2_composite.tif
#sec=${dataDir}/BanksEast_RapidEye_20120801_reorthotile_analytic_sr_composite.tif
#python ${py} ${ref} ${sec}


for reg in BanksEast hotWC Tuk; do
	echo $reg

	ls ${reg}*/*composite.tif > ${reg}_img_list.txt
	#ls ${reg}*/*mosaic.tif >> ${reg}_img_list.txt
	ref=${dataDir}/${reg}*/*2022*composite.tif
	python ${py} ${ref} ${reg}_img_list.txt

done

reg=WillowRiver
ls ${reg}*/*composite.tif > ${reg}_img_list.txt
ls ${reg}*/*mosaic.tif >> ${reg}_img_list.txt
ref=${dataDir}/${reg}*/*2022*composite.tif
python ${py} ${ref} ${reg}_img_list.txt
